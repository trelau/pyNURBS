from __future__ import division

from numpy import array, cross, diff, float64
from numpy.linalg import norm

from pynurbs.config import Settings
from pynurbs.geometry.methods.divide import (extract_bezier_curve,
                                             extract_nurbs_curve,
                                             split_bezier_curve,
                                             split_bezier_surface)
from pynurbs.geometry.methods.geom_utils import (dehomogenize_array1d,
                                                 dehomogenize_array2d,
                                                 tessellate_cp_net)
from pynurbs.geometry.methods.modify import decompose_curve


def polygon_area2D(xy):
    """
    Calculate the signed area of a 2D (planar) polygon using the shoelace
    formula.

    :param array_like xy: Ordered array of points defining polygon boundary.

    :return: Area of polygon.
    :rtype: float
    """
    xy = array(xy, dtype=float64)
    nm = xy.shape

    x0, y0 = xy[0]
    xn, yn = xy[-1]
    a = xn * y0 - x0 * yn
    for i in range(0, nm[0] - 1):
        a += xy[i, 0] * xy[i + 1, 1]
        a -= xy[i + 1, 0] * xy[i, 1]
    a *= 0.5

    return a


def triangle_area(t):
    """
    Calculate the area of a triangle given by 3 points.

    :param array_like t: Array containing triangle vertices.
    """
    t = array(t, dtype=float64)
    return 0.5 * norm(cross(t[1] - t[0], t[2] - t[0]))


def cp_lengths(n, cp):
    """
    Calculate the linear length of the control polygon and its chord.

    :param int n: Number of control points - 1.
    :param ndarray cp: Dehomogenized control points.

    :return: Linear length of control polygon and its chord (lp, lc)
    :rtype: tuple
    """
    lp = sum(norm(diff(cp, axis=0), axis=1))
    lc = norm(cp[-1] - cp[0])
    return lp, lc


def cp_net_areas(n, m, cp):
    """
    Calculate the surface areas of the control point net by tessellating the
    control point net into triangles.

    :param int n: Control points - 1 in u-direction.
    :param int m: Control points - 1 in v-direction.
    :param ndarray cp: Dehomogenized control points.

    :return: Surface area of tessellated control point net and the area of the
        base plane (anet, ap).
    :rtype: tuple

    **Note:** The *base area* if the area of the two triangles formed by the
    four corner points of the control net.
    """
    # Tessellate the surface.
    verts, triangles = tessellate_cp_net(n, m, cp)
    ntri = triangles.shape[0]
    # Area of control net.
    anet = 0.
    for i in range(0, ntri):
        v1 = verts[triangles[i, 0]]
        v2 = verts[triangles[i, 1]]
        v3 = verts[triangles[i, 2]]
        anet += 0.5 * norm(cross(v2 - v1, v3 - v1))
    # Area of base plane.
    ap = (0.5 * norm(cross(cp[n, 0] - cp[0, 0], cp[n, m] - cp[0, 0])) +
          0.5 * norm(cross(cp[n, m] - cp[0, 0], cp[0, m] - cp[0, 0])))
    return anet, ap


def arc_length_bezier(cpw, n, u0=0., u1=1.):
    """
    Estimate the arc length of a Bezier curve.

    :param ndarray cpw: Control points.
    :param int n: Number of control points - 1 (degree).
    :param float u0: Starting parameter.
    :param float u1: Ending parameter.

    :return: Approximate arc length.
    :rtype: float
    """
    gtol = Settings.gtol
    arc_length = array([0.], dtype=float64)

    # Extract curve between u0 and u1 if needed.
    if u0 > 0. or u1 < 1.:
        cpw, _, _ = extract_bezier_curve(cpw, n, u0, u1)

    def _arc_length(cpw_i):
        cp, _ = dehomogenize_array1d(n, cpw_i)
        lp, lc = cp_lengths(n, cp)
        if abs(lp - lc) <= gtol:
            arc_length[0] += (2. * lc + (n - 1.) * lp) / (n + 1.)
        else:
            aw, _, _, bw, _, _ = split_bezier_curve(cpw_i, n, 0.5)
            _arc_length(aw)
            _arc_length(bw)

    _arc_length(cpw)
    return arc_length[0]


def surface_area_bezier(n, m, cpw):
    """
    Estimate the surface area of Bezier surface.

    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param ndarray cpw: Control points.

    :return: Approximate area of surface.
    :rtype: float
    """
    gtol = Settings.gtol
    surface_area = array([0.], dtype=float64)

    def _area(cpw_i):
        cp, _ = dehomogenize_array2d(n, m, cpw_i)
        anet, ap = cp_net_areas(n, m, cp)
        if abs(anet - ap) <= gtol:
            surface_area[0] += (4. * ap + (n + m + n * m - 3) * anet) / (
                    (n + 1) * (m + 1))
        else:
            qw1, qw2, _, _ = split_bezier_surface(cpw_i, n, m, u=0.5)
            qw3, qw4, _, _ = split_bezier_surface(qw1, n, m, v=0.5)
            qw5, qw6, _, _ = split_bezier_surface(qw2, n, m, v=0.5)
            _area(qw3)
            _area(qw4)
            _area(qw5)
            _area(qw6)

    _area(cpw)
    return surface_area[0]


def arc_length_nurbs(n, p, uk, cpw, u0=0., u1=1., a=0., b=1.):
    """
    Estimate the arc length of a NURBS curve between parameters *u0* and *u1*.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.
    :param float u0: Starting parameter.
    :param float u1: Ending parameter.
    :param float a: Lower domain of curve.
    :param float b: Upper domain of curve.

    :return: Estimated arc length.
    :rtype: float
    """
    # Extract curve between u0 and u1 if needed.
    if u0 > a or u1 < b:
        uk, cpw = extract_nurbs_curve(n, p, uk, cpw, u0, u1)
        n = cpw.shape[0] - 1
    # Decompose into Bezier segments
    nb, qw, _ = decompose_curve(n, p, uk, cpw)
    # Find arc length of each segment.
    d = 0.
    for i in range(nb):
        d += arc_length_bezier(qw[i], p)
    return d
