from numpy import array, dot, float64, ndarray

from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.misc import is_array_type
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


class Plane(Geometry):
    """
    Plane.

    The Plane class defines an infinite plane in Cartesian space.

    :param p0: Origin of plane.
    :type p0: :class:`.Point` or array_like
    :param vn: Normal vector of plane.
    :type vn: :class:`.Vector`
    :param vu: Vector along plane x-axis.
    :type vu: :class:`.Vector`
    :param vv: Vector along plane y-axis.
    :type vv: :class:`.Vector`

    :var p0: Plane origin.
    :type p0: :class:`.Point`
    :var vn: Plane normal vector.
    :type vn: :class:`.Vector`
    :var vu: Vector along plane x-axis.
    :type vu: :class:`.Vector`
    :var vv: Vector along plane y-axis.
    :type vv: :class:`.Vector`

    ..warning::
        The vectors *vu*, *vv*, and *vn* should all be orthogonal to each
        other. This is currently not checked at object initialization.
    """

    def __init__(self, p0, vn, vu, vv):
        super(Plane, self).__init__('plane')
        self._p0 = p0
        self._vn = vn
        self._vu = vu
        self._vv = vv

    def __call__(self, u, v, w=0., rtype='Point'):
        return self.eval(u, v, w, rtype)

    @property
    def p0(self):
        return self._p0

    @property
    def vn(self):
        return self._vn

    @property
    def vu(self):
        return self._vu

    @property
    def vv(self):
        return self._vv

    def local_to_global_param(self, d, *args):
        """
        Convert parameter(s) from local domain 0 <= u,v <= 1 to global domain
        a <= u,v <= b.

        :param str d: Direction ('u' or 'v').
        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        if len(args) == 1:
            return args[0]
        return args

    def global_to_local_param(self, d, *args):
        """
        Convert surface parameter(s) from global domain a <= u,v <= b to local
        domain 0 <= u,v <= 1.

        :param str d: Direction ('u' or 'v').
        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        if len(args) == 1:
            return args[0]
        return args

    def eval(self, u=0., v=0., w=0., rtype='Point', *args, **kwargs):
        """
        Evaluate point on the plane at (u, v).

        :param u:
        :param v:
        :param w:
        :param rtype:

        :return:
        """
        pnt = (self._p0.xyz + u * self._vu.ijk + v * self._vv.ijk +
               w * self._vn.ijk)
        if is_array_type(rtype):
            return pnt
        else:
            return Point(pnt)

    def norm(self, u, v, rtype='Vector', **kwargs):
        """
        Compute the plane normal.

        :param float u: Parametric point.
        :param float v: Parametric point.
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray).


        :return:
        """
        vnorm = self._vn.ijk
        if is_array_type(rtype):
            return vnorm
        p0 = self.eval(u, v)
        return Vector(vnorm, p0)

    def dist2pnt(self, p):
        """
        Calculate the signed distance to a point.

        :param p: Point to calculate distance to.
        :type p: :class:`.Point` or array_like

        :return: Distance to point.
        :rtype: float

        .. note::
            A positive distance means the point is in the same direction as
            the plane normal vector. A negative distance means that the
            point is in the opposite direction of the plane normal vector.
        """
        if isinstance(p, Point):
            p = p.xyz
        elif isinstance(p, (tuple, list, ndarray)):
            p = array(p, dtype=float64)
        return dot(self._vn.ijk, p - self._p0.xyz)

    def offset(self, u=0., v=0., w=1.):
        """
        Offset the plane.

        :param u:
        :param v:
        :param w:

        :return:
        """
        p0 = self.eval(u, v, w)
        return Plane(p0, self._vn, self._vu, self._vv)
