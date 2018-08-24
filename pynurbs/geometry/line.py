from numpy import array, float64
from numpy.linalg import norm

from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.misc import is_array_type
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


class Line(Geometry):
    """
    Line.

    The Line class defines an infinite straight line in Cartesian space
    defined by an origin and direction vector.

    :param p0: Origin of line.
    :type p0: :class:`.Point`
    :param v: Line direction vector.
    :type v: :class:`.Vector`

    :var p0: Origin of line.
    :type p0: :class:`.Point`
    :var v: Line direction vector.
    :type v: :class:`.Vector`
    """

    def __init__(self, p0, v):
        super(Line, self).__init__('line')
        self._p0 = p0
        self._v = v

    def __call__(self, u, rtype='Point'):
        return self.eval(u, rtype)

    @property
    def p0(self):
        return self._p0

    @property
    def v(self):
        return self._v

    @property
    def a(self):
        return 0.

    @property
    def b(self):
        return self.length

    @property
    def knots(self):
        return array([self.a, self.b], dtype=float64)

    @property
    def length(self):
        return norm(self._v)

    @staticmethod
    def local_to_global_param(*args):
        """
        For a line this method returns the input arguments and is implemented
        only for consistency with curve_like objects.

        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        if len(args) == 1:
            return args[0]
        return args

    @staticmethod
    def global_to_local_param(*args):
        """
        For a line this method returns the input arguments and is implemented
        only for consistency with curve_like objects.

        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        if len(args) == 1:
            return args[0]
        return args

    def eval(self, u, rtype='Point', *args, **kwargs):
        """
        Evaluate point on line.

        :param float u: Parameter (-inf < u < inf).
        :param str rtype: Option to return a NumPy array or :class:`.Point`
            instance (rtype = 'ndarray' or 'Point').

        :return: Point on line.
        :rtype: :class:`.Point` or ndarray
        """
        pnt = self._p0.xyz + u * self._v.ijk
        if is_array_type(rtype):
            return pnt
        else:
            return Point(pnt)

    def deriv(self, u, rtype='Vector'):
        """
        Compute the line derivative. This only returns the line direction
        vector with magnitude, and is implemented only for consistency with
        Bezier and NURBS curve objects.

        :param float u: Parametric point (-inf < u < inf).
        :param str rtype: Option to return a NumPy array or :class:`.Vector`
            instance (rtype = 'ndarray' or 'Vector').

        :return: Line derivative.
        :rtype: :class:`.Vector` or ndarray
        """
        if is_array_type(rtype):
            return self._v.vxyz
        else:
            p0 = self.eval(u)
            return Vector(self._v.vxyz, p0)

    def extract(self, u0, u1, *args, **kwargs):
        """
        Extract a line.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.

        :return: Curve between *u0* and *u1*.
        :rtype: :class:`.NurbsCurve`
        """
        p0 = self.eval(u0)
        p1 = self.eval(u1)
        v = Vector(p1.xyz - p0.xyz, p0)
        return Line(p0, v)

    def arc_length(self, u0=0., u1=1., *args, **kwargs):
        """
        Estimate the arc length between line parameters.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.

        :return: Arc length of line between *u0* and *u1*.
        :rtype: float
        """
        p0 = self.eval(u0, rtype='ndarray')
        p1 = self.eval(u1, rtype='ndarray')
        return norm(p1 - p0)

    @staticmethod
    def param_at_arc_length(x, *args, **kwargs):
        """
        Estimate the line parameter at a given arc length.

        :param float x: Arc length.

        :return: Global line parameter at arc length.
        :rtype: float
        """
        return float64(x)

    def point_at_arc_length(self, x, rtype='Point', *args, **kwargs):
        """
        Estimate a line point at a given arc length.

        :param float x: Arc length.
        :param str rtype: Option to return a NumPy array or Point instance
            (rtype = 'Point' or 'ndarray').

        :return: Line point at arc length.
        :rtype: :class:`.Point` or ndarray
        """
        return self.eval(x, rtype)
