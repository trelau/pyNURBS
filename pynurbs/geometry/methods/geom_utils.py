from __future__ import division

from numpy import (zeros, float64, int32, cross, dot, hstack, array, arccos,
                   rad2deg)
from numpy.linalg import norm

from pynurbs.config import Settings
from pynurbs.geometry.point import Point


class GeomUtils(object):
    """
    Geometry utilities.
    """

    @staticmethod
    def points_to_array(pnts):
        """
         Convert the list of :class:`.Point` instances into a NumPy array.

        :param list pnts: List of :class:`.Point` instances.

        :return: NumPy array of points.
        :rtype: ndarray
        """
        return array(pnts, dtype=float64)

    @staticmethod
    def array_to_points1d(arr):
        """
        Convert the NumPy array to a list of :class:`.Point` instances.

        :param ndarray arr: NumPy array of points.

        :return: List of :class:`.Point` instances.
        :rtype: list

        If homogeneous points are provided, then the points should be in an
        *n* x 4 array and the Point coordinates and weight will be
        determined by:

        p = (x/w, y/w, z/w)
        """
        return array_to_points1d(arr)

    @staticmethod
    def array_to_points2d(arr):
        """
        Convert the NumPy array to a list of :class:`.Point` instances.

        :param ndarray arr: NumPy array of points.

        :return: n x m list of :class:`.Point` instances.
        :rtype: list

        If homogeneous points are provided, then the points should be in an
        *n* x *m* x 4 array and the Point instance coordinates will be
        determined by:

        p = (x/w, y/w, z/w)
        """
        return array_to_points2d(arr)

    @staticmethod
    def equivalence_points1d(pnts, tol=None):
        """
        Equivalence the 1D list of :class:`.Point` instances.

        :param list pnts: List of :class:`.Point` instances to equivalence.
        :param float tol: Tolerance. If *None* is provided then *gtol* from the
            default settings will be used.

        :return: List of unique :class:`.Point' instances.
        :rtype: list
        """
        return equivalence_points1d(pnts, tol)

    @staticmethod
    def equivalence_array1d(arr, tol=None):
        """
        Equivalence an 1D array.

        :param ndarray arr: Input array.
        :param float tol: Tolerance. If *None* is provided then *gtol* from the
            default settings will be used.

        :return: Array of unique points.
        :rtype: ndarray
        """
        return equivalence_array1d(arr, tol)

    @staticmethod
    def dehomogenize_array1d(n, cpw):
        """
        Convert the 1D homogeneous control points *cpw* to non-homogeneous
        coordinates (x/w, y/w, z/w).


        :param int n: Number of control points - 1.
        :param ndarray cpw: Control points.

        :return: Dehomogeneous control points *cp* and weights *w* (cp, w).
        :rtype: tuple
        """
        return dehomogenize_array1d(n, cpw)

    @staticmethod
    def dehomogenize_array2d(n, m, cpw):
        """
        Convert the 2D homogeneous control net *cpw* to non-homogeneous
        coordinates (x/w, y/w, z/w).


        :param int n: Number of control points - 1 in u-direction.
        :param int m: Number of control points - 1 in v-direction.
        :param ndarray cpw: Control net.

        :return: Dehomogeneous control points *cp* and weights *w* (cp, w).
        :rtype: tuple
        """
        return dehomogenize_array2d(n, m, cpw)

    @staticmethod
    def local_to_global_param(a, b, *args):
        """
        Convert parameter(s) from local domain 0 <= u <= 1 to global domain
        a <= u <= b.

        :param float a: Lower global parameter.
        :param float b: Upper global parameter.
        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        return local_to_global_param(a, b, args)

    @staticmethod
    def global_to_local_param(a, b, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0 <= u <= 1.

        :param float a: Lower global parameter.
        :param float b: Upper global parameter.
        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        return global_to_local_param(a, b, args)

    @staticmethod
    def barycentric_params(point, verts, constrain=True):
        """
        Determine the barycentric coordinates of the point with respect to the
        triangle defined by three vertices.

        :param ndarray point: Point.
        :param ndarray verts: Vertices of triangle as 3 x 3 array.
        :param bool constrain: Option to force the coordinates to be between
            0 and 1.

        :return: Barycentric coordinates of point (u, v)*where 0 <= u,v <= 1
            and a point on the triangle is,

        .. math::

            t(u, v) = (1 - u - v)v_{0} + uv_{1} + vv_{2}

        :rtype: tuple
        """
        return barycentric_params(point, verts, constrain)

    @staticmethod
    def is_tri_degen(t, tol=0.):
        """
        Check to see if the triangle is degenerate (area <= tol).

        :param array_like t: Array containing triangle vertices.
        :param float tol: Tolerance.

        :return: *True* if triangle is degenerate, *False* if not.
        :rtype: bool
        """
        t = array(t, dtype=float64)
        a = 0.5 * norm(cross(t[1] - t[0], t[2] - t[0]))
        if a <= tol:
            return True
        return False

    @staticmethod
    def point_in_triangle(p, t, tol=None):
        """
        Test to see if a point is inside a triangle. The point is first
        projected to the plane of the triangle for this test.

        :param ndarray p: Point inside triangle.
        :param ndarray t: Triangle vertices.
        :param float tol: Barycentric coordinate check.

        :return: *True* if point is inside triangle, *False* if not.
        :rtype: bool
        """
        return point_inside_triangle(p, t, tol)

    @staticmethod
    def point_in_polygon(p2d, poly2d, tol=None):
        """
        Check to see if a 2-D point is in a 2-D polygon using the winding
        number test.

        :param array_like p2d: Point in 2-D.
        :param array_like poly2d: Closed, ordered array defining the boundary
            of the polygon. The first and last points should be coincident.
        :param float tol: Tolerance for checking first/last point equivalence.

        :return: *True* if point is in polygon, *False* if not.
        :rtype: bool
        """
        wn = winding_number(p2d, poly2d, tol)
        return wn != 0


def array_to_points1d(arr):
    n = arr.shape
    # Detect if weighted points are provided.
    if n[1] == 4:
        # p = (x/w, y/w, z/w).
        return [Point(arr[i, :-1] / arr[i, -1]) for i in range(n[0])]
    elif n[1] == 3:
        # p = (x,y,z).
        return [Point(arr[i]) for i in range(n[0])]
    return []


def array_to_points2d(arr):
    nm = arr.shape
    # Detect if weighted points are provided.
    pnts = []
    for i in range(nm[0]):
        pi = []
        for j in range(nm[1]):
            if nm[2] == 4:
                # p = (x/w, y/w, z/w).
                pi.append(Point(arr[i, j, :-1] / arr[i, j, -1]))
            elif nm[2] == 3:
                # p = (x,y,z).
                pi.append(Point(arr[i, j]))
        pnts.append(pi)
    return pnts


def equivalence_points1d(pnts, tol=None):
    if len(pnts) <= 1:
        return pnts
    if tol is None:
        tol = Settings.gtol

    pout = [pnts[0]]
    for p in pnts:
        if any([p.is_equal(pi, tol) for pi in pout]):
            continue
        else:
            pout.append(p)
    return pout


def equivalence_array1d(arr, tol=None):
    if tol is None:
        tol = Settings.gtol
    arr = array(arr, dtype=float64)
    d0, d1 = arr.shape[0:2]

    arr_eqv = zeros(arr.shape, dtype=float64)
    arr_eqv[0] = arr[0]
    npts = 1
    for i in range(1, d0):
        unique = True
        for j in range(0, npts):
            d = norm(arr_eqv[j] - arr[i])
            if d <= tol:
                unique = False
                break
        if unique:
            arr_eqv[npts] = arr[i]
            npts += 1
    return arr_eqv[:npts]


def homogenize_array1d(cp, w):
    _w = w.reshape(-1, 1)
    return hstack((cp * _w, _w))


def homogenize_array2d(cp, w):
    n, m, _ = cp.shape
    cpw = zeros((n, m, 4), dtype=float64)
    for i in range(0, n):
        for j in range(0, m):
            cpw[i, j, :3] = cp[i, j] * w[i, j]
            cpw[i, j, 3] = w[i, j]
    return cpw


def dehomogenize_array1d(n, cpw):
    w = cpw[:, -1]
    cp = cpw[:, :-1] / cpw[:, -1].reshape(-1, 1)
    return cp, w


def dehomogenize_array2d(n, m, cpw):
    cp = zeros((n + 1, m + 1, 3), dtype=float64)
    w = zeros((n + 1, m + 1), dtype=float64)
    for i in range(0, n + 1):
        for j in range(0, m + 1):
            w[i, j] = cpw[i, j, -1]
            cp[i, j] = cpw[i, j, :-1] / w[i, j]
    return cp, w


def are_points_equal(p0, p1, tol=None):
    if tol is None:
        tol = Settings.gtol
    if norm(p0 - p1) <= tol:
        return True
    return False


def tessellate_cp_net(n, m, cp):
    """
    Tessellate the control net of the surface as triangles.

    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param ndarray cp: Control net.

    :return: Tessellated control net as an array of vertices and a
        connectivity matrix (verts, triangles).
    :rtype: tuple

    .. note::
        * The shape of *verts* is [(n + 1) x (m + 1), 3] where the first axis
          will be the index of the vertex, and the second axis is the
          x, y, and z coordinates.
        * The shape of *triangles* is [n x m x 2, 3], where the first axis is
          the index of the triangle, and the second axis are the three
          vertices that define the triangle.
    """
    nvert = (n + 1) * (m + 1)
    ntri = n * m * 2
    # Build the vertex matrix by stacking the control net in the u-direction
    # for each column. Build a temporary matrix to track the i, j indices and
    # the vertex index.
    verts = zeros((nvert, 3), dtype=float64)
    triangles = zeros((ntri, 3), dtype=int32)
    vert_id = zeros((n + 1, m + 1), dtype=int32)
    v = 0
    for j in range(0, m + 1):
        for i in range(0, n + 1):
            verts[v, :] = cp[i, j]
            vert_id[i, j] = v
            v += 1
    # Build triangle matrix in u-direction for each column.
    t = 0
    for j in range(0, m):
        for i in range(0, n):
            v1 = vert_id[i, j]
            v2 = vert_id[i + 1, j]
            v3 = vert_id[i + 1, j + 1]
            v4 = vert_id[i, j + 1]
            triangles[t, :] = [v1, v2, v3]
            triangles[t + 1, :] = [v1, v3, v4]
            t += 2
    return verts, triangles


def local_to_global_param(a, b, *args):
    if args[0] is None:
        return None
    local_u = args
    is_list = True
    if len(local_u) == 1:
        is_list = False
    global_u = []
    for ui in local_u:
        if ui < 0.:
            ui = 0.
        if ui > 1.:
            ui = 1.
        global_u.append(a + ui * (b - a))
    if is_list:
        return global_u
    else:
        return global_u[0]


def global_to_local_param(a, b, *args):
    if args[0] is None:
        return None
    global_u = args
    is_list = True
    if len(global_u) == 1:
        is_list = False
    local_u = []
    for ui in global_u:
        u = (ui - a) / (b - a)
        if u < 0.:
            u = 0.
        if u > 1.:
            u = 1.
        local_u.append(u)
    if is_list:
        return local_u
    else:
        return local_u[0]


def is_curve_flat(n, cp, tol=None):
    """
    Check flatness of curve by measuring the distance of each control point to
    the chord line formed between the first and last control points.

    :param int n: Number of control points - 1 in u-direction.
    :param ndarray cp: Dehomogeneous control points.
    :param float tol: Tolerance. If *None* is provided then *gtol* from the
        default settings will be used.

    :return: *True* if curve is flat, *False* if not.
    :rtype: bool
    """
    if tol is None:
        tol = Settings.ftol

    dline = norm(cp[n] - cp[0])
    if dline <= 0.:
        # Curve must be closed. Use the distance between points.
        for i in range(1, n):
            if norm(cp[0] - cp[i]) > tol:
                return False
        return True
    for i in range(1, n):
        v0 = cp[i] - cp[0]
        v1 = cp[i] - cp[n]
        d = norm(cross(v0, v1)) / dline
        if d > tol:
            return False
    return True


def is_surface_flat(n, m, cp, tol=None):
    """
    Check flatness of surface by measuring the distance of each control
    point to the chord line formed between the first and last control
    points. This is performed for each row/column of the control net and fails
    if any point fails the test.

    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param ndarray cp: Dehomogeneous control points.
    :param float tol: Tolerance. If *None* is provided then *gtol* from the
        default settings will be used.

    :return: *True* if surface is flat, *False* if not.
    :rtype: bool

    .. note::
        Method not yet valid for closed surfaces.
    """
    if tol is None:
        tol = Settings.ftol

    for j in range(0, m + 1):
        if not is_curve_flat(n, cp[:, j], tol):
            return False
    for i in range(0, n + 1):
        if not is_curve_flat(m, cp[i, :], tol):
            return False
    return True


def barycentric_params(p, t, constrain=True):
    # eps = Settings.eps
    v01 = t[1] - t[0]
    v02 = t[2] - t[0]
    vp = p - t[0]
    d01 = dot(v01, v01)
    d12 = dot(v01, v02)
    d02 = dot(v02, v02)
    dp1 = dot(vp, v01)
    dp2 = dot(vp, v02)
    denom = d01 * d02 - d12 * d12
    u = (d02 * dp1 - d12 * dp2) / denom
    v = (d01 * dp2 - d12 * dp1) / denom
    if not constrain:
        return u, v
    if u < 0.:
        u = 0
    elif u > 1.:
        u = 1.
    if v < 0.:
        v = 0
    elif v > 1.:
        v = 1.
    return u, v


def invert_point_in_triangle(p, tri3d, tri2d, constrain=True):
    # Barycentric parameters.
    ub, vb = barycentric_params(p, tri3d, constrain)
    # Convert to 2D.
    u2d, v2d = (1. - ub - vb) * tri2d[0] + ub * tri2d[1] + vb * tri2d[2]
    return u2d, v2d


def point_inside_triangle(p, t, tol=None):
    """
    Test to see if a point is inside a triangle. The point is first
    projected to the plane of the triangle for this test.

    :param ndarray p: Point inside triangle.
    :param ndarray t: Triangle vertices.
    :param float tol: Tolerance for barycentric coordinate check.

    :return: *True* if point is inside triangle, *False* if not.
    :rtype: bool
    """
    if tol is None:
        tol = Settings.ptol
    v01 = t[1] - t[0]
    v02 = t[2] - t[0]
    vp = p - t[0]
    d01 = dot(v01, v01)
    d12 = dot(v01, v02)
    d02 = dot(v02, v02)
    dp1 = dot(vp, v01)
    dp2 = dot(vp, v02)
    denom = d01 * d02 - d12 * d12
    if denom == 0.:
        return False
    u = (d02 * dp1 - d12 * dp2) / denom
    v = (d01 * dp2 - d12 * dp1) / denom
    if u >= -tol and v >= -tol and u + v <= 1. + tol:
        return True
    return False


def winding_number(p2d, poly2d, tol=None):
    """
    Calculate the winding number of the point and the polygon (2-D only).

    :param array_like p2d:
    :param poly2d:
    :param float tol:

    :return: Winding number.
    :rtype: int
    """
    if tol is None:
        tol = Settings.ptol

    p2d = array(p2d, dtype=float64)
    poly2d = array(poly2d, dtype=float64)

    # Make sure polygon is closed.
    if norm(poly2d[0] - poly2d[-1]) > tol:
        poly2d[-1] = poly2d[0]

    # Define is_left method needed for winding number.
    def _is_left(uv0, uv1):
        return ((uv1[0] - uv0[0]) * (p2d[1] - uv0[1]) -
                (p2d[0] - uv0[0]) * (uv1[1] - uv0[1]))

    # Calculate winding number.
    wn = 0
    n = poly2d.shape[0]
    for i in range(0, n - 1):
        if poly2d[i, 1] <= p2d[1]:
            if poly2d[i + 1, 1] > p2d[1]:
                if _is_left(poly2d[i], poly2d[i + 1]) > 0.:
                    wn += 1
        else:
            if poly2d[i + 1, 1] <= p2d[1]:
                if _is_left(poly2d[i], poly2d[i + 1]) < 0.:
                    wn -= 1
    return wn


def angle_between_vecs(v1, v2):
    """
    Calculate the angle between two vectors.
    """
    v1 = array(v1, dtype=float64)
    v2 = array(v2, dtype=float64)
    x = dot(v1, v2) / (norm(v1) * norm(v2))
    if x > 1.:
        x = 1.
    elif x < -1.:
        x = -1.
    return rad2deg(arccos(x))
