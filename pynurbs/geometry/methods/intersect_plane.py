from __future__ import division

from numpy import cross, dot
from numpy.linalg import norm

from pynurbs.config import Settings
from pynurbs.geometry.line import Line
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


def intersect_plane_plane(plane1, plane2):
    """
    Find the intersection of two planes.

    :param plane1: Plane 1 to intersect.
    :type plane1: :class:`.Plane`
    :param plane2: Plane 2 to intersect.
    :type plane2: :class:`.Plane`

    :return: Intersection line of the two planes. Returns *None* if planes
        are parallel.
    :rtype: :class:`.Line`
    """
    # Global parameters.
    n1 = plane1.vn.ijk
    n2 = plane2.vn.ijk
    p1 = plane1.p0.xyz
    p2 = plane2.p0.xyz

    # Cross product to find line vector.
    v = cross(n1, n2)
    vn = norm(v)
    # Return None if planes are parallel.
    if vn <= Settings.atol:
        return None
    v /= vn

    # Find point along intersection vector using intersection of 3 planes.
    d1 = dot(n1, p1)
    d2 = dot(n2, p2)
    n2v = cross(n2, v)
    p0 = Point((d1 * n2v + d2 * cross(v, n1)) / dot(n1, n2v))
    line = Line(p0, Vector(v, p0))
    return line
