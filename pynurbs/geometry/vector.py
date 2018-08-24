from __future__ import division

from numpy import add, array, cross, dot, float64, subtract
from numpy.linalg import norm

from pynurbs.geometry.geom import Geometry


class Vector(Geometry):
    """
    Vector.

    The Vector class defines direction and magnitude in Cartesian space.

    :param array_like v: The vector components.
    :param origin: Origin of vector.
    :type origin: :class:`.Point`

    :var origin: Origin of vector.
    :type origin: :class:`.Point`
    :var ndarray vxyz: Vector array.
    :var float mag: Vector magnitude.
    :var ndarray ijk: Unit vector array.
    """

    def __init__(self, v, origin):
        super(Vector, self).__init__('vector')
        self._vxyz = array(v, dtype=float64)
        self._p0 = origin
        self._mag = norm(self._vxyz)
        self._ijk = self._vxyz / self._mag

    def __array__(self, dtype=float64):
        return array(self._vxyz, dtype=dtype)

    def __iter__(self):
        for elm in self._vxyz:
            yield elm

    def __len__(self):
        return len(self._vxyz)

    def __getitem__(self, item):
        return self._vxyz[item]

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def vxyz(self):
        return self._vxyz

    @property
    def mag(self):
        return self._mag

    @property
    def ijk(self):
        return self._ijk

    @property
    def origin(self):
        return self._p0

    def add(self, v):
        """
        Add vector to current instance.

        :param v: Other vector.
        :type v: :class:`.Vector`

        :return: Vector formed after addition.
        :rtype: :class:`.Vector`
        """
        return Vector(self._vxyz + v.vxyz, self._p0)

    def sub(self, v):
        """
        Subtract vector from current instance.

        :param v: Other vector.
        :type v: :class:`.Vector`

        :return: Vector formed after subtraction.
        :rtype: :class:`.Vector`
        """
        return Vector(self._vxyz - v.vxyz, self._p0)

    def dot(self, v):
        """
        Find the dot product of the two vectors.

        :param v: Other vector.
        :type v: :class:`.Vector`

        :return: Dot product of two vectors.
        :rtype: float
        """
        return dot(self._vxyz, v.vxyz)

    def cross(self, v):
        """
        Find the cross product of the two vectors using right-hand rule.

        :param v: Other vector.
        :type v: :class:`.Vector`

        :return: Cross product of two vectors.
        :rtype: :class:`.Vector`
        """
        return Vector(cross(self._vxyz, v.vxyz), self._p0)
