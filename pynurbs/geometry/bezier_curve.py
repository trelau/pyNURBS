from __future__ import division

from numpy import array, float64, ndarray, ones, zeros

from pynurbs.geometry.bounding_box import BoundingBox
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.calculate import arc_length_bezier
from pynurbs.geometry.methods.divide import (extract_bezier_curve,
                                             split_bezier_curve)
from pynurbs.geometry.methods.evaluate import (bezier_curve_points,
                                               check_param,
                                               deCasteljau1,
                                               find_mult_knots,
                                               rat_curve_derivs)
from pynurbs.geometry.methods.geom_utils import (are_points_equal,
                                                 dehomogenize_array1d,
                                                 global_to_local_param,
                                                 homogenize_array1d,
                                                 local_to_global_param)
from pynurbs.geometry.methods.misc import (is_array_like, is_array_type,
                                           is_local_domain)
from pynurbs.geometry.methods.modify import elevate_bezier_curve_degree
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


class BezierCurve(Geometry):
    """
    Bezier curve.

    The Beizer curve class is defined entirely by a set of control points and
    can be either rational or non-rational. The number of control points will
    determine the degree of the curve.

    :var float a: Lower bound of curve domain.
    :var float b: Upper bound of curve domain.
    :var int p: Degree.
    :var int n: Number of control points - 1.
    :var int m: Number of knots - 1 (m = n + p + 1).
    :var ndarray uk: Knot vector.
    :var ndarray cp: Control points.
    :var ndarray w: Weights.
    :var ndarray cpw: Homogenious control points.
    """

    def __init__(self):
        super(BezierCurve, self).__init__('bezier_curve')
        self._n = None
        self._cp = None
        self._w = None
        self._cpw = None
        self._uk = None
        self._p = None
        self._m = None
        self._u0 = 0.
        self._u1 = 1.

    def __call__(self, u, rtype='Point', domain='local'):
        return self.eval(u, rtype, domain)

    @property
    def a(self):
        return self._u0

    @property
    def b(self):
        return self._u1

    @property
    def p(self):
        return self._p

    @property
    def m(self):
        return self._m

    @property
    def n(self):
        return self._n

    @property
    def uk(self):
        return self._uk

    @property
    def knots(self):
        return self.get_mult_and_knots()[-1]

    @property
    def mult(self):
        return self.get_mult_and_knots()[1]

    @property
    def cp(self):
        return self._cp

    @property
    def w(self):
        return self._w

    @property
    def cpw(self):
        return self._cpw

    @property
    def data(self):
        return self._n, self._p, self._uk, self._cpw

    @property
    def length(self):
        return self.arc_length()

    @property
    def is_closed(self):
        if self._cpw is None:
            return False
        return are_points_equal(self._cpw[0], self._cpw[-1])

    @property
    def bbox(self):
        return self.get_bbox()

    def clear(self):
        """
        Clear all values assigned to curve.
        """
        self._u0 = 0.
        self._u1 = 1.
        self._n = None
        self._cp = None
        self._w = None
        self._cpw = None
        self._uk = None
        self._p = None
        self._m = None

    def set_domain(self, a=0., b=1.):
        """
        Set the domain of the Bezier curve between [a, b] where a < b.

        :param float a: Lower domain.
        :param float b: Upper domain.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if self._p is None or a >= b:
            return False
        self._u0 = a
        self._u1 = b
        self._uk = zeros(2 * (self._p + 1), dtype=float64)
        self._uk[:self._p + 1] = a
        self._uk[self._p + 1:] = b
        self._m = self._uk.size - 1
        return True

    def set_cp(self, cp):
        """
        Set control points (weights are set to 1).

        :param cp: Control points.
        :type: array_like

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        self._cp = array(cp, dtype=float64)
        self._n = self._cp.shape[0] - 1
        self._p = self._n
        self._w = ones(self._cp.shape[0], dtype=float64)
        self._cpw = homogenize_array1d(self._cp, self._w)
        self.set_domain()
        return True

    def set_w(self, w):
        """
        Set the weights of the control points.

        :param w: Weights of control points (w > 0.)
        :type w: array_like

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if self._cp is None or self._w is None:
            return False
        if isinstance(w, (tuple, list, ndarray)):
            w = array(w, dtype=float64)
            if w.shape != self._w.shape:
                return False
            indx = w <= 0.
            w[indx] = 1.
            self._w[:] = w
            self._cpw[:, :] = homogenize_array1d(self._cp, self._w)
            return True
        return False

    def set_cpw(self, cpw):
        """
        Set control points and weights.

        :param array_like cpw: Array of homogeneous control points. This array
            is dehomogenized and used to set *cp* and *w*.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        cpw = array(cpw, dtype=float64)
        cp, w = dehomogenize_array1d(cpw.shape[0] - 1, cpw)
        return self.set_cp(cp) and self.set_w(w)

    def modify_cp(self, indx, cpi):
        """
        Modify a control point of the curve at the given index.

        :param int indx: Index of control points ( 0 <= indx <= n).

        :param cpi: New control point location.
        :type cpi: :class:`.Point` or array_like

        :return: *True* if modified, *False* if not.
        :rtype: bool
        """
        if not isinstance(self._cp, ndarray):
            return False
        cpi = array(cpi, dtype=float64)
        if indx > self._n:
            indx = -1
        dim = cpi.shape[0]
        self._cp[indx, :dim] = cpi
        self._cpw[indx, :dim] = self._cp[indx, :dim] * self._w[indx]
        return True

    def local_to_global_param(self, *args):
        """
        Convert parameter(s) from local domain 0 <= u <= 1 to global domain
        a <= u <= b.

        :param args: Local parameter(s).

        :return: Global parameter(s).
        :rtype: float or list of floats
        """
        return local_to_global_param(self._u0, self._u1, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0 <= u <= 1.

        :param args: Global parameter(s).

        :return: Local parameter(s).
        :rtype: float or list of floats
        """
        return global_to_local_param(self._u0, self._u1, *args)

    def eval(self, u, rtype='Point', domain='local'):
        """
        Evaluate curve at parametric point.

        :param float u: Parametric point.
        :param str rtype: Option to return a NumPy array or Point instance
            (rtype = 'Point' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Point on curve.
        :rtype: :class:`.Point` or ndarray
        """
        if is_array_like(u):
            return self.eval_params(u, rtype, domain)

        if not is_local_domain(domain):
            u = self.global_to_local_param(u)
        pw = deCasteljau1(self._cpw, self._p, u)
        pnt = pw[:-1] / pw[-1]
        if is_array_type(rtype):
            return pnt
        else:
            return Point(pnt)

    def eval_params(self, ulist, rtype='Point', domain='local'):
        """
        Evaluate the curve at multiple parameters.

        :param array_like ulist: Curve parameters.
        :param str rtype: Option to return a NumPy array or list of Point
            instances (rtype = 'Point' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Points at parameters.
        :rtype: List of :class:`.Point` instances or NumPy array.
        """
        if not is_array_like(ulist):
            return self.eval(ulist, rtype, domain)

        if not is_local_domain(domain):
            ulist = map(self.global_to_local_param, ulist)
        pnts = bezier_curve_points(self._n, self._cpw, ulist)
        if is_array_type(rtype):
            return pnts
        else:
            return [Point(pi) for pi in pnts]

    def deriv(self, u, k, d=None, rtype='Vector', domain='local'):
        """
        Compute the curve derivative.

        :param float u: Parametric point.
        :param int k: Compute and return the *k*-th derivative.
        :param int d: Highest derivative to compute. If *d* = *None*, then only
            the *k*-th derivative will be computed.
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Curve derivative.
        :rtype: :class:`.Vector` or ndarray
        """
        if self._cp is None:
            return None
        if d is None:
            d = k
        # Use NURBS rational curve derivative to compute Bezier derivative.
        # Convert to global since using a NURBS method.
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        der = rat_curve_derivs(self._n, self._p, self._uk, self._cpw, u, d)
        if is_array_type(rtype):
            return der[k]
        else:
            p0 = Point(der[0])
            return Vector(der[k], p0)

    def split(self, u=0.5, domain='local'):
        """
        Split the curve.

        :param float u: Parametric point to split curve at.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Two new Bezier curves *(c1, c2)*.
        :rtype: tuple
        """
        if not is_local_domain(domain):
            u = self.global_to_local_param(u)
        qw1, a1, b1, qw2, a2, b2 = split_bezier_curve(self._cpw, self._n, u,
                                                      self._u0, self._u1)
        c1, c2 = BezierCurve(), BezierCurve()
        c1.set_cpw(qw1)
        c2.set_cpw(qw2)
        c1.set_domain(a1, b1)
        c2.set_domain(a2, b2)
        return c1, c2

    def extract(self, u0=0., u1=1., domain='local'):
        """
        Extract a curve.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Curve between *u0* and *u1*.
        :rtype: :class:`.BezierCurve`
        """
        if not is_local_domain(domain):
            u0, u1 = self.global_to_local_param(u0, u1)
        qw, a, b = extract_bezier_curve(self._cpw, self._n, u0, u1, self._u0,
                                        self._u1)
        c = BezierCurve()
        c.set_cpw(qw)
        c.set_domain(a, b)
        return c

    def elevate_degree(self, deg=None, inplace=False):
        """
        Elevate degree of curve.

        :param int deg: Degree to elevate curve to. If *None* is provided, the
            curve degree will be elevated once.
        :param bool inplace: Option to return a new curve or modify the curve
            in place.

        :return: New curve with elevated degree.
        :rtype: :class:`.BezierCurve`
        """
        if self._cp is None:
            return None
        if deg is None:
            deg = self._p + 1
        qw = elevate_bezier_curve_degree(self._n, deg, self._cpw)
        if inplace:
            self.clear()
            self.set_cpw(qw)
        else:
            c = BezierCurve()
            c.set_cpw(qw)
            return c

    def decompose(self):
        """
        Decompose curve into Bezier segments. For a Bezier curve this simply
        returns a copy of itself.

        :return: Number and list of Bezier curve segments (nb, curves).
        :rtype: tuple
        """
        return 1, [self.copy()]

    def get_mult_and_knots(self):
        """
        Get an array of unique knots values for the curve and their
        multiplicities.

        :return: Number of unique knots, multiplicities, and the knots
            (nu, um, uq).
        :rtype: tuple
        """
        return find_mult_knots(self._n, self._p, self._uk)

    def get_bbox(self):
        """
        Get curve bounding box.

        :return: Curve bounding box.
        :rtype: :class:`.BoundingBox`
        """
        bbox = BoundingBox()
        bbox.add_curve(self)
        return bbox

    def check_param(self, u):
        """
        Check that the parameter is within the global domain of the knot vector
        or is within tolerance of a unique knot value. Use
        :class:`.CompareFloats` to compare floats.

        :param float u: Global parameter.

        :return: Parameter within global domain or near interior knot value.
        :rtype: float
        """
        return check_param(self._n, self._p, self._uk, u)

    def arc_length(self, u0=0., u1=1., domain='local'):
        """
        Estimate the arc length between curve parameters.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.
        :param str domain: Option to use local or global domain.

        :return: Arc length of curve between *u0* and *u1*.
        :rtype: float
        """
        if not is_local_domain(domain):
            u0, u1 = self.global_to_local_param(u0, u1)

        return arc_length_bezier(self.cpw, self.n, u0, u1)
