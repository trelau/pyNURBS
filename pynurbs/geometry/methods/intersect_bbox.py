from __future__ import division

from numpy import dot, float64, array, divide, errstate, mean
from numpy.linalg import norm


def bboxes_intersect(bbox1, bbox2, tol=0.):
    """
    Test if the bounding boxes intersect.

    :param BoundingBox bbox1: Bounding box 1.
    :param BoundingBox bbox2: Bounding box 2.
    :param float tol: Tolerance for bounding boxes.

    :return: *True* if they intersect, *False* if not.
    :rtype: bool
    """
    bounds1 = bbox1.bounds
    bounds2 = bbox2.bounds
    # Adjust bounding box by the tolerance.
    bounds1[:, 0] -= tol
    bounds1[:, 1] += tol
    bounds2[:, 0] -= tol
    bounds2[:, 1] += tol
    for i in range(3):
        if bounds1[i, 0] > bounds2[i, 1] or bounds2[i, 0] > bounds1[i, 1]:
            return False
    return True


def bbox_intersects_plane(bbox, plane, tol=0.):
    """
    Test if the bounding box intersects a plane.

    :param BoundingBox bbox: Bounding box.
    :param Plane plane: Plane.
    :param float tol: Tolerance for bounding box.

    :return: *True* if they intersect, *False* if not.
    :rtype: bool
    """
    p0 = plane.p0.xyz
    pnorm = plane.vn.ijk
    bounds = bbox.bounds
    # Adjust bounding box by the tolerance.
    bounds[:, 0] -= tol
    bounds[:, 1] += tol
    # Find centroid and radius of bounding box.
    cg = mean(bounds, axis=1)
    r = norm(bounds[:, 1] - bounds[:, 0]) / 2.
    # Distance from cg to plane.
    d = abs(dot(pnorm, cg - p0))
    return d <= r


def bbox_ray_intersect(bbox, p0, v, t0=None, t1=None, tol=0.):
    """
    Test to see if a ray intersects a bounding box.

    :param bbox: Axis-aligned bounding box.
    :type bbox: :class:`.BoundingBox`
    :param array_like p0: Origin of ray.
    :param array_like v: Direction vector of ray.
    :param float t0: Lower domain of valid intersection interval.
    :param float t1: Upper domain of valid intersection interval.
    :param float tol: Tolerance for bounding box.

    :return: *True* if ray intersects bounding box, *False* if not.
    :rtype: bool
    """
    # Get arrays.
    bounds = bbox.bounds
    p0 = array(p0, float64)
    v = array(v, float64)

    # Adjust bounding box by the tolerance.
    bounds[:, 0] -= tol
    bounds[:, 1] += tol

    # Set t0 and t1 to +/- Inifity is they are not.
    with errstate(all='ignore'):
        if t0 is None:
            t0 = divide(1., -0.)
        if t1 is None:
            t1 = divide(1., 0.)

        # Calculate inverse direction vectors.
        divx = divide(1., v[0])
        divy = divide(1., v[1])
        divz = divide(1., v[2])

        # x-direction
        if divx >= 0.:
            tmin = (bounds[0, 0] - p0[0]) * divx
            tmax = (bounds[0, 1] - p0[0]) * divx
        else:
            tmin = (bounds[0, 1] - p0[0]) * divx
            tmax = (bounds[0, 0] - p0[0]) * divx

        # y-direction
        if divy >= 0.:
            tymin = (bounds[1, 0] - p0[1]) * divy
            tymax = (bounds[1, 1] - p0[1]) * divy
        else:
            tymin = (bounds[1, 1] - p0[1]) * divy
            tymax = (bounds[1, 0] - p0[1]) * divy

        if tmin > tymax or tymin > tmax:
            return False

        if tymin > tmin:
            tmin = tymin
        if tymax < tmax:
            tmax = tymax

        # z-direction
        if divz >= 0.:
            tzmin = (bounds[2, 0] - p0[2]) * divz
            tzmax = (bounds[2, 1] - p0[2]) * divz
        else:
            tzmin = (bounds[2, 1] - p0[2]) * divz
            tzmax = (bounds[2, 0] - p0[2]) * divz

        if tmin > tzmax or tzmin > tmax:
            return False

        if tzmin > tmin:
            tmin = tzmin
        if tzmax < tmax:
            tmax = tzmax

        return (tmin < t1) and (tmax > t0)
