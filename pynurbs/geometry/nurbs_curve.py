from numpy import array, float64, ndarray, ones

from pynurbs.geometry.bezier_curve import BezierCurve
from pynurbs.geometry.bounding_box import BoundingBox
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.calculate import arc_length_nurbs
from pynurbs.geometry.methods.divide import (extract_nurbs_curve,
                                             split_nurbs_curve)
from pynurbs.geometry.methods.evaluate import (check_param, curve_point,
                                               curve_points,
                                               find_mult_knots,
                                               rat_curve_derivs)
from pynurbs.geometry.methods.geom_utils import (are_points_equal,
                                                 dehomogenize_array1d,
                                                 global_to_local_param,
                                                 homogenize_array1d,
                                                 local_to_global_param)
from pynurbs.geometry.methods.misc import (is_array_like, is_array_type,
                                           is_local_domain)
from pynurbs.geometry.methods.modify import (adjust_knot_vector,
                                             curve_knot_ins,
                                             decompose_curve,
                                             elevate_nurbs_curve_degree,
                                             refine_knot_vect_curve,
                                             reverse_nurbs_curve)
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


class NurbsCurve(Geometry):
    """
    NURBS curve.

    The NURBS (Non-Uniform Rational B-Spline) curve is defined by control
    points, a degree, and a knot vector. No parameters are set at
    initialization.

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
        super(NurbsCurve, self).__init__('nurbs_curve')
        self._n = None
        self._cp = None
        self._w = None
        self._cpw = None
        self._uk = None
        self._p = None
        self._m = None
        self._a = 0.
        self._b = 1.

    def __call__(self, u, rtype='Point', domain='local'):
        return self.eval(u, rtype, domain)

    @property
    def a(self):
        return self._a

    @property
    def b(self):
        return self._b

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
        self._a = 0.
        self._b = 1.
        self._n = None
        self._cp = None
        self._w = None
        self._uk = None
        self._p = None
        self._m = None

    def set_domain(self, a=0., b=1.):
        """
        Adjust the knot vector *uk* to be between [a, b], where a < b. The
        interior knots are adjusted to maintain original parametrization.

        :param float a: Lower domain.
        :param float b: Upper domain.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        qk = adjust_knot_vector(self._n, self._p, self._uk, a, b)
        self.set_knots(qk)
        self._a = self._uk[self._p]
        self._b = self._uk[self._n + 1]
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
        self._w = ones(self._cp.shape[0], dtype=float64)
        self._cpw = homogenize_array1d(self._cp, self._w)
        self.set_deg()
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

    def set_knots(self, uk):
        """
        Set knot vector of NURBS curve.

        :param array_like uk: Knot vector.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        self._uk = array(uk, dtype=float64)
        self._m = self._uk.size - 1
        if self._p is not None and self._n is not None:
            self._a = self._uk[self._p]
            self._b = self._uk[self._n + 1]
        self.set_deg()
        return True

    def set_deg(self, p=None):
        """
        Set the degree of the curve. If the control points and knots are
        already set, then p = m - n - 1.

        :param int p: Degree.
        """
        if self._uk is not None:
            self._m = self._uk.size - 1
        if self._cp is not None and self._uk is not None:
            self._p = self._m - self._n - 1
            self._a = self._uk[self._p]
            self._b = self._uk[self._n + 1]
            return True
        elif p is not None:
            self._p = int(p)
            return True
        return False

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
        """
        return local_to_global_param(self._a, self._b, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0 <= u <= 1.

        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        return global_to_local_param(self._a, self._b, *args)

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

        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        pnt = curve_point(self._n, self._p, self._uk, self._cpw, u)
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
        :rtype: List of :class:`.Point` instances or NumPy array
        """
        if not is_array_like(ulist):
            return self.eval(ulist, rtype, domain)

        if is_local_domain(domain):
            ulist = map(self.global_to_local_param, ulist)
        pnts = curve_points(self._n, self._p, self._uk, self._cpw, ulist)
        if is_array_type(rtype):
            return pnts
        else:
            return [Point(pi) for pi in pnts]

    def deriv(self, u, k, d=None, rtype='Vector', domain='local'):
        """
        Compute the *k* -th derivative at point *u*.

        :param float u: Parametric point.
        :param int k: Derivative to return (0 <= k <= d).
        :param int d: Highest derivative to compute. If *d* = *None*, then only
            the *k* th derivative will be computed.
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Curve *k* -th derivative.
        :rtype: :class:`.Vector` or ndarray
        """
        if self._cp is None:
            return None
        if d is None:
            d = k
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        der = rat_curve_derivs(self._n, self._p, self._uk,
                               self._cpw, u, d)
        if is_array_type(rtype):
            return der[k]
        else:
            p0 = Point(der[0])
            return Vector(der[k], p0)

    def split(self, u=0.5, domain='local'):
        """
        Split the curve into two segments at *u*.

        :param float u: Parametric point to split curve at.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Two new NURBS curves *(c1, c2)*.
        :rtype: tuple
        """
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        uk1, qw1, uk2, qw2 = split_nurbs_curve(self._n, self._p, self._uk,
                                               self._cpw, u)
        c1, c2 = NurbsCurve(), NurbsCurve()
        c1.set_cpw(qw1)
        c1.set_knots(uk1)
        c1.set_deg(self._p)
        c2.set_cpw(qw2)
        c2.set_knots(uk2)
        c2.set_deg(self._p)
        return c1, c2

    def extract(self, u0, u1, domain='local'):
        """
        Extract a curve.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Curve between *u0* and *u1*.
        :rtype: :class:`.NurbsCurve`
        """
        if is_local_domain(domain):
            u0, u1 = self.local_to_global_param(u0, u1)
        uq, qw = extract_nurbs_curve(self._n, self._p, self._uk, self._cpw, u0,
                                     u1)
        c = NurbsCurve()
        c.set_cpw(qw)
        c.set_knots(uq)
        c.set_deg(self._p)
        return c

    def elevate_degree(self, deg=None, inplace=False):
        """
        Elevate degree of curve.

        :param int deg: Degree to elevate curve to. If *None* is provided, the
            curve degree will be elevated once.
        :param bool inplace: Option to return a new curve or modify the curve
         in place.

        :return: New curve with elevated degree.
        :rtype: :class:`.NurbsCurve`
        """
        if self._cp is None:
            return None
        if deg is None:
            t = 1
        else:
            t = deg - self._p
        nq, uq, qw = elevate_nurbs_curve_degree(self._n, self._p, self._uk,
                                                self._cpw, t)
        if inplace:
            self.clear()
            self.set_cpw(qw)
            self.set_knots(uq)
            self.set_deg(deg)
        else:
            c = NurbsCurve()
            c.set_cpw(qw)
            c.set_knots(uq)
            c.set_deg(deg)
            return c

    def insert_knot(self, u, r=1, inplace=False, domain='local'):
        """
        Insert knot.

        :param float u: Knot value to insert.
        :param int r: Number of times to insert knot (r + s <= p, where *s*
            is knot multiplicity).
        :param bool inplace: Option to return a new curve or modify the curve
            in place.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: New curve after knot insertion.
        :rtype: :class:`.NurbsCurve`
        """
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        uk, cpw = curve_knot_ins(self._n, self._p, self._uk, self._cpw, u, r)
        if inplace:
            self.set_cpw(cpw)
            self.set_knots(uk)
            return True
        else:
            c = NurbsCurve()
            c.set_cpw(cpw)
            c.set_knots(uk)
            c.set_deg(self._p)
            return c

    def refine_knots(self, x, inplace=False, domain='local'):
        """
        Refine the knot vector by inserting the elements of *x*.

        :param array_like x: Knot vector to insert.
        :param bool inplace: Option to return new curve or modify
            existing curve.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: New curve after knot refinement.
        :rtype: :class:`.NurbsCurve`
        """
        if is_local_domain(domain):
            x = self.local_to_global_param(*x)
        if isinstance(x, (tuple, list, ndarray)):
            x = array(x, dtype=float64)
        elif isinstance(x, (int, float)):
            x = array([x], dtype=float64)
        uq, qw = refine_knot_vect_curve(self._n, self._p, self._uk,
                                        self._cpw, x)
        if inplace:
            self.set_cpw(qw)
            self.set_knots(uq)
            return True
        else:
            c = NurbsCurve()
            c.set_cpw(qw)
            c.set_knots(uq)
            c.set_deg(self._p)
            return c

    def decompose(self):
        """
        Decompose curve into Bezier segments. The domain of the Bezier segments
        will be derived from the NURBS curve knot vector.

        :return: Number and list of Bezier curve segments (nb, curves).
        :rtype: tuple
        """
        nb, qw, ab = decompose_curve(self._n, self._p, self._uk, self._cpw)
        curves = []
        for i in range(0, nb):
            c = BezierCurve()
            c.set_cpw(qw[i])
            c.set_domain(ab[i, 0], ab[i, 1])
            curves.append(c)
        return nb, curves

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
        Return a bounding box instance for the curve.

        :return: Curve bounding box.
        :rtype: :class:`.BoundingBox`
        """
        bbox = BoundingBox()
        bbox.add_curve(self)
        return bbox

    def reverse(self, inplace=True):
        """
        Reverse direction of curve.

        :param bool inplace: Option to return new curve or modify
            existing curve.

        :return: Reversed NURBS curve.
        :rtype: :class:`.NurbsCurve`
        """
        uq, qw = reverse_nurbs_curve(self._n, self._p, self._uk, self._cpw)
        if inplace:
            self.set_knots(uq)
            self.set_cpw(qw)
            return True
        c = NurbsCurve()
        c.set_cpw(qw)
        c.set_knots(uq)
        c.set_deg(self._p)
        return c

    def reverse_param(self, u, domain='local'):
        """
        Reverse the parameter.

        :param float u: Original parameter.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Reversed parameter.
        :rtype: float
        """
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        u = -u + self._a + self._b
        if not is_local_domain(domain):
            return u
        return self.global_to_local_param(u)

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
        if is_local_domain(domain):
            u0, u1 = self.local_to_global_param(u0, u1)
        return arc_length_nurbs(self.n, self.p, self.uk, self.cpw, u0, u1,
                                self.a, self.b)
