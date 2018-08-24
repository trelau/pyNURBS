from __future__ import division

from numpy import cross, dot
from numpy.linalg import norm


def point_to_line(p0, p1, pi):
    """
    Calculate the distance from a point to a line.

    :param array_like p0: Starting point of line.
    :param array_like p1: Ending point of line.
    :param array_like pi: Point to calculate distance to.

    :return: Distance from point to line.
    :rtype: float
    """
    denom = norm(p1 - p0)
    if denom == 0.:
        return norm(pi - p0)
    numer = norm(cross(pi - p0, pi - p1))
    return numer / denom


def point_to_plane(p0, vn, pi, signed=False):
    """
    Calculate the distance from a point to a plane.

    :param array_like p0: Origin of plane.
    :param array_like vn: Normal vector of plane.
    :param array_like pi: Point to calculate distance to.
    :param bool signed: Option to returned signed distance (*True*), or the
        absolute value (*False*).

    :return: Distance from point to plane.
    :rtype: float
    """
    vnorm = vn / norm(vn)
    d = dot(vnorm, pi - p0)
    if signed:
        return d
    return abs(d)
