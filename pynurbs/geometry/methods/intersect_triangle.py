from __future__ import division

from numpy import dot, cross, zeros, float64
from numpy.linalg import norm

from pynurbs.config import Settings


def intersect_triangle_ray(v0, v1, v2, origin, ray, test_cull=False, eps=None):
    """
    Compute the intersection of a triangle with a ray.

    :param ndarray v0: Vertex 0 of triangle.
    :param ndarray v1: Vertex 1 of triangle.
    :param ndarray v2: Vertex 2 of triangle.
    :param ndarray origin: Ray origin.
    :param ndarray ray: Ray direction (unit vector).
    :param bool test_cull: Option to cull back facing triangles.
    :param float eps: Tolerance for ray-triangle intersection. If *None* is
        provided, then the *eps* value from the settings class is used.

    :return: *None* if the ray does not intersect the triangle. If an
        intersection is found, then the parameter *t* on the ray and and
        parameters (u,v) on the triangle are returned (t, u, v).
    :rtype: tuple

    *Reference:* Moller, T., Trumbore, B., "Fast Minimum Storage Ray Triangle
    Intersection."
    """
    if eps is None:
        eps = Settings.atol

    e1 = v1 - v0
    e2 = v2 - v0
    pvec = cross(ray, e2)
    det = dot(e1, pvec)
    if test_cull:
        if det <= eps:
            return None
        tvec = origin - v0
        u = dot(tvec, pvec)
        if u <= -eps or u >= det + eps:
            return None
        qvec = cross(tvec, e1)
        v = dot(ray, qvec)
        if v <= -eps or u + v >= det:
            return None
        t = dot(e2, qvec)
        inv_det = 1. / det
        t *= inv_det
        u *= inv_det
        v *= inv_det
    else:
        if abs(det) < eps:
            return None
        inv_det = 1. / det
        tvec = origin - v0
        u = dot(tvec, pvec) * inv_det
        if u <= -eps or u >= 1. + eps:
            return None
        qvec = cross(tvec, e1)
        v = dot(ray, qvec) * inv_det
        if v <= -eps or u + v >= 1. + eps:
            return None
        t = dot(e2, qvec) * inv_det
    return t, u, v


def intersect_triangles(t1, t2, tol=None, eps=None):
    """
    Find the intersection segment of two triangles.

    :param ndarray t1: Array containing vertices of triangle 1 (3 x 3).
    :param ndarray t2: Array containing vertices of triangle 2 (3 x 3).
    :param float tol: Optional tolerance to use for length of intersection
        segment. If none is provided then the default *gtol* will be used.
    :param float eps: Tolerance for ray-triangle intersection. If *None* is
        provided, then the *eps* value from the settings class is used.

    :return: The number of intersection points and the intersection points that
        form the intersection segment between two triangles (if any) as a
        (2 x 3) array (ni, pi).
    :rtype: tuple
    """
    if tol is None:
        tol = Settings.gtol
    if eps is None:
        eps = Settings.atol

    # Find the normal vectors of the planes of the triangles.
    vn1 = cross(t1[1] - t1[0], t1[2] - t1[0])
    vn2 = cross(t2[1] - t2[0], t2[2] - t2[0])

    # TODO Test for parallel triangles.
    # if abs(dot(vn1, vn2)) - norm(vn1) * norm(vn2) <= Settings.atol:
    #     return 0, []

    # Test the two triangles for possible intersect by checking the signed
    # distance of their vertices to the other's plane.
    test = False
    d0 = dot(vn1, t2[0] - t1[0])
    for i in range(1, 3):
        di = dot(vn1, t2[i] - t1[0])
        if d0 * di <= 0.:
            test = True
            break
    if not test:
        return 0, []

    test = False
    d0 = dot(vn2, t1[0] - t2[0])
    for i in range(1, 3):
        di = dot(vn2, t1[i] - t2[0])
        if d0 * di <= 0.:
            test = True
            break
    if not test:
        return 0, []

    # Intersect edges of triangle 1 with triangle 2.
    verts = zeros((3, 3), dtype=float64)
    edges = zeros((3, 3), dtype=float64)
    edges[0] = t1[1] - t1[0]
    verts[0] = t1[0]
    edges[1] = t1[2] - t1[0]
    verts[1] = t1[0]
    edges[2] = t1[2] - t1[1]
    verts[2] = t1[1]
    pnts = zeros((2, 3), dtype=float64)
    npts = 0
    for i in range(3):
        tuv = intersect_triangle_ray(t2[0], t2[1], t2[2], verts[i], edges[i])
        if tuv is None or tuv[0] < -eps or tuv[0] > 1. + eps:
            continue
        t, u, v = tuv
        # pi = (1. - u - v) * t2[0] + u * t2[1] + v * t2[2]
        pi = verts[i] + t * edges[i]
        pnts[npts] = pi
        npts += 1
        if npts == 2:
            # Keep going if there are vertices still to check.
            di = norm(pnts[1] - pnts[0])
            if di <= tol and i < 2:
                npts -= 1
            elif di > tol:
                return npts, pnts
            elif di <= tol and i == 2:
                npts -= 1

    # Intersect edges of triangle 2 with triangle 1.
    edges[0] = t2[1] - t2[0]
    verts[0] = t2[0]
    edges[1] = t2[2] - t2[0]
    verts[1] = t2[0]
    edges[2] = t2[2] - t2[1]
    verts[2] = t2[1]
    for i in range(3):
        tuv = intersect_triangle_ray(t1[0], t1[1], t1[2], verts[i], edges[i])
        if tuv is None or tuv[0] < -eps or tuv[0] > 1. + eps:
            continue
        t, u, v = tuv
        # pi = (1. - u - v) * t1[0] + u * t1[1] + v * t1[2]
        pi = verts[i] + t * edges[i]
        pnts[npts] = pi
        npts += 1
        if npts == 2:
            # Keep going if there are vertices still to check.
            di = norm(pnts[1] - pnts[0])
            if di <= tol and i < 2:
                npts -= 1
            elif di > tol:
                return npts, pnts
            elif di <= tol and i == 2:
                npts -= 1

    return 0, []


def intersect_triangle_plane(t, p0, pnorm, tol=None):
    """
    Find the intersection segment of a plane and a triangle.

    :param ndarray t: Array containing vertices of triangle (3 x 3).
    :param ndarray p0: Origin of plane.
    :param ndarray pnorm: Normal vector of plane.
    :param float tol: Optional tolerance to use for length of intersection
        segment. If none is provided then the default *gtol* will be used.

    :return: Intersection segment between triangle and plane (if any) as a
        (2 x 3) array.
    :rtype: ndarray
    """
    if tol is None:
        tol = Settings.gtol

    pnorm /= norm(pnorm)

    # Check for parallel triangle and plane. Return no intersection if
    # they're parallel.
    tri_norm = cross(t[1] - t[0], t[2] - t[0])
    tri_norm /= norm(tri_norm)
    if 1. - abs(dot(pnorm, tri_norm)) <= 1.0e-12:
        return 0, zeros((2, 3), dtype=float64)

    # Test for possible intersection by checking the signed distance of the
    # vertices to the plane.
    test = False
    d0 = dot(pnorm, t[0] - p0)
    for i in range(1, 3):
        di = dot(pnorm, t[i] - p0)
        if d0 * di <= 0.:
            test = True
            break
    if not test:
        return 0, zeros((2, 3), dtype=float64)

    # Intersect the edges of the triangle with the plane.
    verts = zeros((3, 3), dtype=float64)
    edges = zeros((3, 3), dtype=float64)
    edges[0] = t[1] - t[0]
    verts[0] = t[0]
    edges[1] = t[2] - t[0]
    verts[1] = t[0]
    edges[2] = t[2] - t[1]
    verts[2] = t[1]
    pnts = zeros((2, 3), dtype=float64)
    npts = 0
    for i in range(3):
        d0 = dot(pnorm, verts[i] - p0)
        d1 = dot(pnorm, verts[i] + edges[i] - p0)
        if d0 * d1 > 0.:
            continue
        if abs(d0) <= tol and abs(d1) <= tol:
            pnts[0] = verts[i]
            pnts[1] = verts[i] + edges[i]
            return 2, pnts
        denom = dot(pnorm, edges[i])
        ui = dot(pnorm, p0 - verts[i]) / denom
        if ui < 0.:
            ui = 0.
        if ui > 1.:
            ui = 1.
        pi = verts[i] + ui * edges[i]
        pnts[npts] = pi
        npts += 1
        if npts == 2:
            # Keep going if there are vertices still to check.
            di = norm(pnts[1] - pnts[0])
            if di <= tol and i < 2:
                npts -= 1
            elif di <= tol and i == 2:
                return 0, zeros((2, 3), dtype=float64)
            else:
                return npts, pnts
    return 0, zeros((2, 3), dtype=float64)
