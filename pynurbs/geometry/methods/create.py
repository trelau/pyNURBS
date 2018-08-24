from __future__ import division

from numpy import argsort, array, cross, dot, float64, mean
from numpy.linalg import svd

from pynurbs.geometry.checker import CheckGeom
from pynurbs.geometry.intersector import IntersectGeom
from pynurbs.geometry.line import Line
from pynurbs.geometry.methods.geom_utils import angle_between_vecs
from pynurbs.geometry.methods.misc import is_local_domain
from pynurbs.geometry.plane import Plane
from pynurbs.geometry.point import Point
from pynurbs.geometry.projector import ProjectGeom
from pynurbs.geometry.system import System
from pynurbs.geometry.vector import Vector


def vector_by_points(p0, p1, p2=None):
    """
    Create a vector defined by two or three points.

    :param p0: Origin of plane.
    :type p0: :class:`.Point` or array_like
    :param p1: Point defining vector *p1* - *p0*.
    :type p1: :class:`.Point` or array_like
    :param p2: Point defining vector *p2* - *p0*.
    :type p2: :class:`.Point` or array_like

    :return: A vector from *p0* to *p1* if only two points are provided,
        or a vector defined by the cross product of *p10* x *p20* if all three
        points are provided.
    :rtype: :class:`.Vector`
    """
    if not isinstance(p0, Point):
        p0 = Point(p0)
    if not isinstance(p1, Point):
        p1 = Point(p1)
    # Cross product if three points are provided.
    if p2 is not None:
        if not isinstance(p2, Point):
            p2 = Point(p2)
        v10 = p1.xyz - p0.xyz
        v20 = p2.xyz - p0.xyz
        vn = cross(v10, v20)
        return Vector(vn, p0)
    # Straight vector if two points are provided.
    return Vector(p1.xyz - p0.xyz, p0)


def vector_by_axis(axis, origin):
    """
    Create a vector along the specified axis.

    :param str axis: Axis ('x', 'y', or 'z').
    :param array_like origin: Origin of vector.

    :return: Vector along specified axis.
    :rtype: :class:`.Vector`
    """
    if not isinstance(origin, Point):
        origin = Point(origin)

    if not isinstance(axis, str):
        return None
    if axis.lower() not in ['x', 'y', 'z']:
        return None

    if axis.lower() in ['x']:
        return Vector([1., 0., 0.], origin)
    if axis.lower() in ['y']:
        return Vector([0., 1., 0.], origin)
    if axis.lower() in ['z']:
        return Vector([0., 0., 1.], origin)


def line_by_points(p0, p1):
    """
    Create a line defined by two points.

    :param p0: Origin of plane.
    :type p0: :class:`.Point` or array_like
    :param p1: Point defining vector *p1* - *p0*.
    :type p1: :class:`.Point` or array_like

    :return: A line defined by an oriign at *p0* and a vector *p10*.
    :rtype: :class:`.Line`
    """
    if not isinstance(p0, Point):
        p0 = Point(p0)
    v = vector_by_points(p0, p1)
    return Line(p0, v)


def plane_by_points(p0, p1, p2):
    """
    Create a plane defined by three points.

    :param p0: Origin of plane.
    :type p0: :class:`.Point` or array_like
    :param p1: Point defining vector *p1* - *p0*.
    :type p1: :class:`.Point` or array_like
    :param p2: Point defining vector *p2* - *p0*.
    :type p2: :class:`.Point` or array_like

    :return: A plane with a normal vector defined by the cross product of
        *p10* x *p20* and an x-axis oriented towards *p1*
    :rtype: :class:`.Plane`
    """
    if not isinstance(p0, Point):
        p0 = Point(p0)
    if not isinstance(p1, Point):
        p1 = Point(p1)
    if not isinstance(p2, Point):
        p2 = Point(p2)
    v10 = p1.xyz - p0.xyz
    v20 = p2.xyz - p0.xyz
    vn = cross(v10, v20)
    vv = cross(vn, v10)
    vu = cross(vv, vn)
    vu = Vector(vu, p0)
    vv = Vector(vv, p0)
    vn = Vector(vn, p0)
    return Plane(p0, vn, vu, vv)


def fit_plane(pnts, tol=None):
    """
    Fit a plane to a scattered set of points.

    :param pnts: Points to fit (at least 3 points are required).
    :type pnts: list of :class:`.Point` instances or array_like
    :param float tol: Tolerance for checking the fit. If *None* is
            provided then the plane will be fit to the points. If a float is
            provided then the plane will not be created if the distance from
            any points is greater than *tol*.

    :return: Plane that best fits data.
    :rtype: :class:`.Plane`
    """
    # Convert points to array.
    pnts = array(pnts, dtype=float64)
    if pnts.shape[0] < 3:
        return None

    # Calculate average to use as the plane origin.
    p0 = mean(pnts, axis=0)

    # Move points to centroid.
    pc = pnts - p0

    # Use SVD.
    u, s, v = svd(pc, False)

    # Check that points are not on a line
    if abs(s[2] - s[1]) <= 1.0e-12:
        return None

    # Find min and max values that define normal vector and major and minor
    # axes of the plane.
    indx = argsort(s)
    vn = v[indx[0]]
    vv = v[indx[1]]
    vu = v[indx[2]]

    # Create plane.
    p0 = Point(p0)
    vu = Vector(vu, p0)
    vv = Vector(vv, p0)
    vn = Vector(vn, p0)
    plane = Plane(p0, vn, vu, vv)

    if tol is None:
        return plane

    # Check distance to each point.
    for pi in pnts:
        if abs(plane.dist2pnt(pi)) > tol:
            return None
    return plane


def plane_by_axes(p0, axes, sys=None):
    """
    Create a plane defined by an origin and standard axes.

    :param p0: Origin of plane.
    :type p0: :class:`.Point` or array_like
    :param axes: Standard axes, one of 'xy', 'xz', or 'yz'.
    :param sys: Reference system for axes.
    :type sys: :class:`.System`

    :return: Plane oriented by axes.
    :rtype: :class:`.Plane`
    """
    if not isinstance(axes, str):
        return None
    if axes.lower() not in ['xy', 'yx', 'xz', 'zx', 'yz', 'zy']:
        return None

    if not isinstance(p0, Point):
        p0 = Point(p0)

    vx = array([1., 0., 0.], dtype=float64)
    vy = array([0., 1., 0.], dtype=float64)
    vz = array([0., 0., 1.], dtype=float64)
    if CheckGeom.is_system(sys):
        vx = sys.vx.ijk
        vy = sys.vy.ijk
        vz = sys.vz.ijk

    if axes.lower() in ['xy', 'yx']:
        p1 = p0.xyz + vx
        p2 = p0.xyz + vy
        return plane_by_points(p0, p1, p2)
    if axes.lower() in ['xz', 'zx']:
        p1 = p0.xyz + vz
        p2 = p0.xyz + vx
        return plane_by_points(p0, p1, p2)
    if axes.lower() in ['yz', 'zy']:
        p1 = p0.xyz + vy
        p2 = p0.xyz + vz
        return plane_by_points(p0, p1, p2)


def planes_by_offset(plane, offset, n):
    """
    Create planes by offsetting an original plane.

    :param plane: Plane to offset.
    :type plane: :class:`.Plane`
    :param float offset: Distance to offset.
    :param int n: Number of planes to generate.

    :return: List of planes offset from original.
    :rtype: list
    """
    if not isinstance(plane, Plane):
        return None
    if n <= 0:
        n = 1
    planes = []
    for i in range(1, n + 1):
        planes.append(plane.offset(w=i * offset))
    return planes


def plane_by_normal(p0, vn):
    """
    Create a plane by an origin and normal vector.

    :param p0: Origin of plane.
    :type p0: :class:`.Point` or array_like
    :param vn: Normal vector of plane.
    :type vn: :class:`.Vector` or array_like

    :return: Plane with given origin and normal vector.
    :rtype: :class:`.Plane`
    """
    if not isinstance(p0, Point):
        p0 = Point(p0)
    if not isinstance(vn, Vector):
        vn = Vector(vn, p0)

    # Try x-axis.
    vu = cross([1., 0., 0.], vn.vxyz)
    if dot(vu, vu) > 0.:
        vv = Vector(cross(vu, vn.vxyz), p0)
        vu = Vector(vu, p0)
        return Plane(p0, vn, vu, vv)

    # Try y-axis.
    vu = cross([0., 1., 0.], vn.vxyz)
    if dot(vu, vu) > 0.:
        vv = Vector(cross(vu, vn.vxyz), p0)
        vu = Vector(vu, p0)
        return Plane(p0, vn, vu, vv)

    # Try z-axis.
    vu = cross([0., 0., 1.], vn.vxyz)
    if dot(vu, vu) > 0.:
        vv = Vector(cross(vu, vn.vxyz), p0)
        vu = Vector(vu, p0)
        return Plane(p0, vn, vu, vv)
    return None


def create_icurve_by_points(surface, p0, p1, isurf=None, trim=True):
    """
    Create an intersection curve between the two points on the surface.
    """
    # Get surface parameters.
    uv0 = ProjectGeom.invert(p0, surface)
    uv1 = ProjectGeom.invert(p1, surface)
    if None in uv0 or None in uv1:
        return None

    # Re-evaluate in case points were not on surface.
    p0 = surface.eval(uv0[0], uv0[1], domain='global')
    p1 = surface.eval(uv1[0], uv1[1], domain='global')

    # Generate intersection plane if no surface is given.
    if not CheckGeom.is_surface_like(isurf):
        vn = surface.norm(uv0[0], uv0[1], domain='global')
        p2 = Point(p0.xyz + vn.ijk)
        isurf = plane_by_points(p0, p1, p2)
        # If points are colinear return None for now.
        if isurf is None:
            return None

    # Intersect the surface.
    si = IntersectGeom.perform(surface, isurf)
    if not si.success:
        return None
    indx = si.curve_nearest_point(p0)
    icrv = si.get_icurve(indx)

    # Project points to curve.
    u0 = ProjectGeom.invert(p0, icrv)
    u1 = ProjectGeom.invert(p1, icrv)
    if None in [u0, u1]:
        return None

    # Reverse if necessary.
    if u0 > u1:
        icrv.reverse(True)
        u0, u1 = -u0 + icrv.a + icrv.b, -u1 + icrv.a + icrv.b

    # Trim if desired and return.
    if not trim:
        return icrv
    return icrv.extract(u0, u1, domain='global')


def points_at_kinks(curve, angle=30., u0=0., u1=1., domain='local'):
    """
    Create points at kinks of a curve.
    """
    if is_local_domain(domain):
        u0, u1 = curve.local_to_global_param(u0, u1)

    # Get knots and multiplicities.
    nu, um, uq = curve.get_mult_and_knots()

    kinks = []
    # If p > 1, create a point wherever the knot multiplicity equals the
    # degree, excluding the endpoints.
    if curve.p > 1:
        for mult, ui in zip(um, uq)[1:-1]:
            if ui <= u0 or ui >= u1:
                continue
            if mult == curve.p:
                pi = curve.eval(ui, domain='global')
                kinks.append((ui, pi))
        return kinks

    # If p = 1, then use the angle criteria and the control points.
    angle_tol = abs(angle)
    n, cp = curve.n, curve.cp
    for i in range(1, n):
        ui = uq[i]
        if ui <= u0 or ui >= u1:
            continue
        # Get previous, current, and next points.
        p0 = cp[i - 1]
        pi = cp[i]
        p1 = cp[i + 1]
        # Calculate angle.
        v1 = pi - p0
        v2 = p1 - pi
        angle = angle_between_vecs(v1, v2)
        if angle > angle_tol:
            pi = curve.eval(ui, domain='global')
            kinks.append((ui, pi))
    return kinks


def create_system_by_points(origin, x_axis, xz_plane):
    """
    Coordinate system by three points.
    """
    origin, x_axis, xz_plane = CheckGeom.to_points([origin, x_axis, xz_plane])
    for p in [origin, x_axis, xz_plane]:
        if not CheckGeom.is_point(p):
            return None
    vx = vector_by_points(x_axis, origin)
    vy = vector_by_points(origin, x_axis, xz_plane)
    vz = vx.cross(vy)
    return System(origin, vx, vy, vz)
