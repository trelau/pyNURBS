from __future__ import division

from numpy import array, cross, float64, ndarray, ones, zeros

from pynurbs.geometry.bezier_surface import BezierSurface
from pynurbs.geometry.bounding_box import BoundingBox
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.calculate import surface_area_bezier
from pynurbs.geometry.methods.divide import (extract_nurbs_isocurve,
                                             extract_nurbs_surface,
                                             split_nurbs_surface)
from pynurbs.geometry.methods.evaluate import (check_param, find_mult,
                                               find_mult_knots,
                                               rat_surface_derivs,
                                               surface_point, surface_points)
from pynurbs.geometry.methods.geom_utils import (are_points_equal,
                                                 dehomogenize_array2d,
                                                 global_to_local_param,
                                                 homogenize_array2d,
                                                 local_to_global_param)
from pynurbs.geometry.methods.misc import (is_array_like, is_array_type,
                                           is_local_domain)
from pynurbs.geometry.methods.modify import (adjust_knot_vector,
                                             decompose_surface,
                                             move_nurbs_surface_seam,
                                             refine_knot_vect_surface,
                                             surface_knot_ins)
from pynurbs.geometry.nurbs_curve import NurbsCurve
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


class NurbsSurface(Geometry):
    """
    NURBS surface.

    The NURBS (Non-Uniform Rational B-Spline) surface is defined by control
    points, a degree, and a knot vector in both the u- and v-directions. No
    parameters are set at initialization.

    :var float au: Lower bound of surface domain in u-direction.
    :var float bu: Upper bound of surface domain in u-direction.
    :var float av: Lower bound of surface domain in v-direction.
    :var float bv: Upper bound of surface domain in v-direction.
    :var int p: Degree in u-direction.
    :var int q: Degree in v-direction.
    :var int n: Number of control points - 1 in u-direction.
    :var int m: Number of control points - 1 in v-direction.
    :var int r: Number of knots - 1 in u-direction.
    :var int s: Number of knots - 1 in v-direction.
    :var ndarray uk: Knot vector in u-direction.
    :var ndarray vk: Knot vector in v-direction.
    :var ndarray cp: Control points.
    :var ndarray w: Weights.
    :var ndarray cpw: Homogeneous control points.
    """

    def __init__(self):
        super(NurbsSurface, self).__init__('nurbs_surface')
        self._n = None
        self._m = None
        self._cp = None
        self._w = None
        self._cpw = None
        self._p = None
        self._q = None
        self._uk = None
        self._vk = None
        self._r = None
        self._s = None
        self._au = 0.
        self._bu = 1.
        self._av = 0.
        self._bv = 1.

    def __call__(self, u, v, rtype='Point', domain='local'):
        return self.eval(u, v, rtype, domain)

    @property
    def au(self):
        return self._au

    @property
    def bu(self):
        return self._bu

    @property
    def av(self):
        return self._av

    @property
    def bv(self):
        return self._bv

    @property
    def p(self):
        return self._p

    @property
    def q(self):
        return self._q

    @property
    def n(self):
        return self._n

    @property
    def m(self):
        return self._m

    @property
    def uk(self):
        return self._uk

    @property
    def vk(self):
        return self._vk

    @property
    def uknots(self):
        return self.get_mult_and_knots('u')[-1]

    @property
    def umult(self):
        return self.get_mult_and_knots('u')[1]

    @property
    def vknots(self):
        return self.get_mult_and_knots('v')[-1]

    @property
    def vmult(self):
        return self.get_mult_and_knots('v')[1]

    @property
    def r(self):
        return self._r

    @property
    def s(self):
        return self._s

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
        return (self._n, self._p, self._uk, self._m, self._q, self._vk,
                self._cpw)

    @property
    def bbox(self):
        return self.get_bbox()

    def clear(self):
        """
        Clear all values assigned to surface.
        """
        self._n = None
        self._m = None
        self._cp = None
        self._w = None
        self._cpw = None
        self._p = None
        self._q = None
        self._uk = None
        self._vk = None
        self._r = None
        self._s = None
        self._au = 0.
        self._bu = 1.
        self._av = 0.
        self._bv = 1.

    def set_domain(self, au=0., bu=1., av=0., bv=1.):
        """
        Adjust the knot vectors to be between [a, b], where a < b. The
        interior knots are adjusted to maintain original parametrization.

        :param float au: Lower bound in u-direction.
        :param float bu: Upper bound in u-direction.
        :param float av: Lower bound in v-direction.
        :param float bv: Upper bound in v-direction.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        # u-direction
        uk = adjust_knot_vector(self._n, self._p, self._uk, au, bu)
        vk = adjust_knot_vector(self._m, self._q, self._vk, av, bv)
        self.set_knots(uk, vk)
        self._au = self._uk[self._p]
        self._bu = self._uk[self._n + 1]
        self._av = self._vk[self._q]
        self._bv = self._vk[self._m + 1]
        return True

    def set_cp(self, cp):
        """
        Set control points (weights are set to 1).

        :param cp: Control point net.
        :type: cp: array_like

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        self._cp = array(cp, dtype=float64)
        nm = self._cp.shape
        self._n = nm[0] - 1
        self._m = nm[1] - 1
        self._w = ones((nm[0], nm[1]), dtype=float64)
        self._cpw = homogenize_array2d(self._cp, self._w)
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
            self._w[:, :] = w
            self._cpw[:, :, :] = homogenize_array2d(self._cp, self._w)
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
        nm = cpw.shape
        cp, w = dehomogenize_array2d(nm[0] - 1, nm[1] - 1, cpw)
        return self.set_cp(cp) and self.set_w(w)

    def set_knots(self, uk=None, vk=None):
        """
        Set the knots of NURBS surface.

        :param uk: Knots to set for u-direction.
        :param vk: Knots to set for v-direction.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if uk is not None:
            self._uk = array(uk, dtype=float64)
            self._r = self._uk.size - 1
        if vk is not None:
            self._vk = array(vk, dtype=float64)
            self._s = self._vk.size - 1
        if self._p is not None and self._n is not None:
            self._au = self._uk[self._p]
            self._bu = self._uk[self._n + 1]
        if self._q is not None and self._m is not None:
            self._av = self._vk[self._q]
            self._bv = self._vk[self._m + 1]
        return self.set_deg()

    def set_deg(self, p=None, q=None):
        """
        Set the degree of the surface. If the control points and knots are
        already set, then p = h - n - 1 and q = k - m - 1.

        :param int p: Degree in u-direction.
        :param int q: Degree in v-direction.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if self._cp is not None:
            self._n = len(self._cp) - 1
            self._m = len(self.cp[0]) - 1
        if self._uk is not None and self._vk is not None:
            self._r = self._uk.size - 1
            self._s = self._vk.size - 1
        if (self._cp is not None and self._uk is not None and
                self._vk is not None):
            self._p = self._r - self._n - 1
            self._q = self._s - self._m - 1
            self._au = self._uk[self._p]
            self._bu = self._uk[self._n + 1]
            self._av = self._vk[self._q]
            self._bv = self._vk[self._m + 1]
            return True
        elif p is not None:
            self._p = int(p)
            return True
        elif q is not None:
            self._q = int(q)
            return True
        return False

    def local_to_global_param(self, d, *args):
        """
        Convert parameter(s) from local domain 0 <= u,v <= 1 to global domain
        a <= u,v <= b.

        :param str d: Direction ('u' or 'v').
        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        if d.lower() in ['u']:
            return local_to_global_param(self._au, self._bu, *args)
        else:
            return local_to_global_param(self._av, self._bv, *args)

    def global_to_local_param(self, d, *args):
        """
        Convert surface parameter(s) from global domain a <= u,v <= b to local
        domain 0 <= u,v <= 1.

        :param str d: Direction ('u' or 'v').
        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        if d.lower() in ['u']:
            return global_to_local_param(self._au, self._bu, *args)
        else:
            return global_to_local_param(self._av, self._bv, *args)

    def eval(self, u, v, rtype='Point', domain='local'):
        """
        Evaluate surface at parametric points.

        :param float u: Parametric point in u-direction.
        :param float v: Parametric point in v-direction.
        :param str rtype: Option to return a NumPy array or Point instance
            (rtype = 'ndarray' or 'Point').
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Point on surface.
        :rtype: :class:`.Point` or ndarray
        """
        if is_array_like(u) and is_array_like(v):
            return self.eval_params(u, v, rtype, domain)

        if is_local_domain(domain):
            u = self.local_to_global_param('u', u)
            v = self.local_to_global_param('v', v)
        pnt = surface_point(self._n, self._p, self._uk,
                            self._m, self._q, self._vk,
                            self._cpw, u, v)
        if is_array_type(rtype):
            return pnt
        else:
            return Point(pnt)

    def eval_params(self, ulist, vlist, rtype='Point', domain='local'):
        """
        Evaluate the surface at multiple parameters.

        :param array_like ulist: Parameters in u-direction.
        :param array_like vlist: Parameters in v-direction.
        :param str rtype: Option to return a NumPy array or Point instance
            (rtype = 'ndarray' or 'Point').
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Points on surface.
        :rtype: List :class:`.Point` instances or ndarray
        """
        if not is_array_like(ulist) or not is_array_like(vlist):
            return self.eval(ulist, vlist, rtype, domain)

        if is_local_domain(domain):
            ulist = [self.local_to_global_param('u', ui) for ui in ulist]
            vlist = [self.local_to_global_param('v', vi) for vi in vlist]
        pnts = surface_points(self._n, self._p, self._uk,
                              self._m, self._q, self._vk,
                              self._cpw, ulist, vlist)
        if is_array_type(rtype):
            return pnts
        else:
            return [Point(pi) for pi in pnts]

    def deriv(self, u, v, k, l, rtype='Vector', domain='local'):
        """
        Compute the surface derivative.

        :param float u: Parametric point.
        :param float v: Parametric point.
        :param int k: Derivative to return in u-direction (0 <= k <= d).
        :param int l: Derivative to return in v-direction (0 <= l <= d).
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray').
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Surface derivative.
        :rtype: :class:`.Vector` or ndarray
        """
        if self._cp is None:
            return None
        if is_local_domain(domain):
            u = self.local_to_global_param('u', u)
            v = self.local_to_global_param('v', v)
        d = k + l
        der = rat_surface_derivs(self._n, self._p, self._uk, self._m, self._q,
                                 self._vk, self._cpw, u, v, d)
        if is_array_type(rtype):
            return der[k, l]
        else:
            p0 = Point(der[0, 0])
            return Vector(der[k, l], p0)

    def norm(self, u, v, rtype='Vector', domain='local'):
        """
        Compute the surface normal.

        :param float u: Parametric point.
        :param float v: Parametric point.
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray).
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Surface normal.
        :rtype: :class:`.Vector` or ndarray
        """
        du = self.deriv(u, v, 1, 0, rtype='ndarray', domain=domain)
        dv = self.deriv(u, v, 0, 1, rtype='ndarray', domain=domain)
        norm = cross(du, dv)
        if is_array_type(rtype):
            return norm
        else:
            p0 = self.eval(u, v, rtype='Point', domain=domain)
            return Vector(norm, p0)

    def tangent(self, u, v, rtype='Vector', domain='local'):
        """
        Surface tangent at (u, v).

        :param float u: Parametric point.
        :param float v: Parametric point.
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray).
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Surface tangent vector.
        :rtype: :class:`.Vector` or ndarray
        """
        du = self.deriv(u, v, 1, 0, rtype='ndarray', domain=domain)
        dv = self.deriv(u, v, 0, 1, rtype='ndarray', domain=domain)
        if is_array_type(rtype):
            return du + dv
        else:
            p0 = self.eval(u, v, 'Point', domain=domain)
            return Vector(du + dv, p0)

    def is_closed(self, d):
        """
        Check if surface is closed along u- or v-direction.

        :param str d: Direction to check for a closed surface ('u' or 'v').

        :return: *True* if surface is closed, *False* if not.
        :rtype: bool
        """
        # Check control points at v=0 and v=1 in the u-direction.
        cpw = self._cpw
        if d.lower() in ['v']:
            for i in range(0, self._n + 1):
                if not are_points_equal(cpw[i, 0], cpw[i, -1]):
                    return False
            return True
        # Check control points at u=0 and u=1 in the v-direction.
        if d.lower() in ['u']:
            for j in range(0, self._m + 1):
                if not are_points_equal(cpw[0, j], cpw[-1, j]):
                    return False
            return True
        return False

    def surface_area(self, u0=0., u1=1., v0=0., v1=1., domain='local'):
        """
        Estimate the surface area. Use parameters to restrict area
        calculation if desired.

        :param float u0: Starting parameter in u-direction.
        :param float u1: Ending parameter in u-direction.
        :param float v0: Starting parameter in v-direction.
        :param float v1: Ending parameter in v-direction.
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Surface area.
        :rtype: float
        """
        s = self.extract(u0, u1, v0, v1, domain)
        area = 0.
        for si in s.decompose():
            area += surface_area_bezier(si.n, si.m, si.cpw)
        return area

    def split(self, u=None, v=None, domain='local'):
        """
        Split NURBS surface.

        :param float u: Location of split (in v-direction).
        :param float v: Location of split (in u-direction).
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Tuple of NURBS surfaces after splitting.
        :rtype: tuple

        .. note::
            If only one parameter is provided, then two surfaces are
            created and the length of the tuple is two (s1, s2). If both
            parameters are provided, four surfaces are created and the length
            of tuple is four (s1, s2, s3, s4).
        """
        if is_local_domain(domain):
            u = self.local_to_global_param('u', u)
            v = self.local_to_global_param('v', v)
        if u is None or v is None:
            # Split in one direction only.
            uk1, vk1, qw1, uk2, vk2, qw2 = split_nurbs_surface(self._n,
                                                               self._p,
                                                               self._uk,
                                                               self._m,
                                                               self._q,
                                                               self._vk,
                                                               self._cpw, u, v)
            s1, s2 = NurbsSurface(), NurbsSurface()
            s1.set_cpw(qw1)
            s1.set_knots(uk1, vk1)
            s2.set_cpw(qw2)
            s2.set_knots(uk2, vk2)
            return s1, s2
        elif u is not None and v is not None:
            # Split in both directions, starting at u (in v-direction).
            uk1, vk1, qw1, uk2, vk2, qw2 = split_nurbs_surface(self._n,
                                                               self._p,
                                                               self._uk,
                                                               self._m,
                                                               self._q,
                                                               self._vk,
                                                               self._cpw, u,
                                                               None)
            # Split first surface at v (in u-direction).
            nm1, nm2 = qw1.shape[:-1], qw2.shape[:-1]
            n1, m1 = nm1[0] - 1, nm1[1] - 1
            uk3, vk3, qw3, uk4, vk4, qw4 = split_nurbs_surface(n1, self._p,
                                                               uk1, m1,
                                                               self._q, vk1,
                                                               qw1, None, v)
            # Split second surface at v (in u-direction).
            n2, m2 = nm2[0] - 1, nm2[1] - 1
            uk5, vk5, qw5, uk6, vk6, qw6 = split_nurbs_surface(n2, self._p,
                                                               uk2, m2,
                                                               self._q, vk2,
                                                               qw2, None, v)
            # Create surfaces.
            s3, s4, s5, s6 = (NurbsSurface(), NurbsSurface(),
                              NurbsSurface(), NurbsSurface())
            s3.set_cpw(qw3)
            s3.set_knots(uk3, vk3)
            s4.set_cpw(qw4)
            s4.set_knots(uk4, vk4)
            s5.set_cpw(qw5)
            s5.set_knots(uk5, vk5)
            s6.set_cpw(qw6)
            s6.set_knots(uk6, vk6)
            return s3, s4, s5, s6

    def isocurve(self, u=None, v=None, domain='local'):
        """
        Extract iso-curve from NURBS surface.

        :param float u: Parameter *u* to extract curve at (in v-direction).
        :param float v: Parameter *v* to extract curve at (in u-direction).
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: NURBS curve along parameter.
        :rtype: :class:`.NurbsCurve`
        """
        if is_local_domain(domain):
            u = self.local_to_global_param('u', u)
            v = self.local_to_global_param('v', v)
        uq, qw = extract_nurbs_isocurve(self._n, self._p, self._uk,
                                        self._m, self._q, self._vk,
                                        self._cpw, u, v)
        c = NurbsCurve()
        c.set_cpw(qw)
        c.set_knots(uq)
        return c

    def extract(self, u0=0., u1=1., v0=0., v1=1., domain='local'):
        """
        Extract NURBS surface bounded by parameters.

        :param float u0: Starting parameter in u-direction.
        :param float u1: Ending parameter in u-direction.
        :param float v0: Starting parameter in v-direction.
        :param float v1: Ending parameter in v-direction.
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Extracted NURBS surface.
        :rtype: :class:`.NurbsSurface`
        """
        if is_local_domain(domain):
            u0 = self.local_to_global_param('u', u0)
            u1 = self.local_to_global_param('u', u1)
            v0 = self.local_to_global_param('v', v0)
            v1 = self.local_to_global_param('v', v1)
        _, _, uq, vq, qw = extract_nurbs_surface(self._n, self._p, self._uk,
                                                 self._m, self._q, self._vk,
                                                 self._cpw, u0, u1, v0, v1)
        s = NurbsSurface()
        s.set_cpw(qw)
        s.set_knots(uq, vq)
        return s

    def insert_knot(self, uv, d, r=1, inplace=False, domain='local'):
        """
        Insert knot in direction.

        :param float uv: Knot value to insert.
        :param str d: Direction to insert knot ('u1 or 'v').
        :param int r: Number of times to insert knot.
        :param bool inplace: Option to return new surface or modify
            existing surface.
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: New surface after knot insertion.
        :rtype: :class:`.NurbsSurface`
        """
        if is_local_domain(domain):
            uv = self.local_to_global_param(d, uv)
        uk, vk, cpw = surface_knot_ins(self._n, self._p, self._uk,
                                       self._m, self._q, self._vk,
                                       self._cpw, d, uv, r)
        if inplace:
            self.set_cpw(cpw)
            self.set_knots(uk, vk)
            return True
        else:
            s = NurbsSurface()
            s.set_cpw(cpw)
            s.set_knots(uk, vk)
            return s

    def refine_knots(self, x, d, inplace=False, domain='local'):
        """
        Refine the knot vector by inserting the elements of *x* in specified
        direction.

        :param array_like x: Knot vector to insert.
        :param str d: Direction to refine knots ('u' or 'v').
        :param bool inplace: Option to return new surface or modify
            existing surface.
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: New surface after knot refinement.
        :rtype: :class:`.NurbsSurface`
        """
        if is_local_domain(domain):
            x = self.local_to_global_param(d, *x)
        if isinstance(x, (tuple, list)):
            x = array(x, dtype=float64)
        elif isinstance(x, (int, float)):
            x = array([x], dtype=float64)
        uk, vk, cpw = refine_knot_vect_surface(self._n, self._p, self._uk,
                                               self._m, self._q, self._vk,
                                               self._cpw, x, d)
        if inplace:
            self.set_cpw(cpw)
            self.set_knots(uk, vk)
            return True
        else:
            s = NurbsSurface()
            s.set_cpw(cpw)
            s.set_knots(uk, vk)
            return s

    def decompose(self):
        """
        Decompose surface into Bezier patches. The domain of the Bezier patches
        will be derived from the NURBS surface knot vectors.

        :return: Number and list of Bezier patches (nbu, nbv, patches).
        :rtype: tuple
        """
        nbu, nbv = 0, 0
        # Decompose in u-direction to get Bezier strips.
        nbu, qwu, abu = decompose_surface(self._n, self._p, self._uk,
                                          self._m, self._q, self._vk,
                                          self._cpw, 'u')

        # Decompose the strips in v-direction to get patches.
        patches = []
        uk = zeros(2 * (self._p + 1), dtype=float64)
        for ii in range(0, nbu):
            cpw = qwu[ii]
            uk[:self._p + 1] = abu[ii, 0]
            uk[self._p + 1:] = abu[ii, 1]
            nbv, qwv, abv = decompose_surface(self._p, self._p, uk,
                                              self._m, self._q, self._vk,
                                              cpw, 'v')
            # Create Bezier surface for each v-direction patch.
            row = []
            for jj in range(0, nbv):
                s = BezierSurface()
                s.set_cpw(qwv[jj])
                s.set_domain(abu[ii, 0], abu[ii, 1], abv[jj, 0], abv[jj, 1])
                row.append(s)
            patches.append(row)
        return nbu, nbv, patches

    def get_mult_and_knots(self, d):
        """
        Get an array of unique knots values for the curve and their
        multiplicities in the specified direction.

        :param str d: Direction to get knot multiplicities and values ('u'
            or 'v').

        :return: Number of unique knots, multiplicities, and the knots
            (nu, um, uq) in specified direction.
        :rtype: tuple
        """
        if d.lower() in ['u']:
            return find_mult_knots(self._n, self._p, self._uk)
        else:
            return find_mult_knots(self._m, self._q, self._vk)

    def get_mult(self, uv, d):
        """
        Get multiplicity of parameter.

        :param float uv: Parameter.
        :param str d: Direction to check ('u' or 'v')

        :return: Multiplicity of parameter.
        :rtype: int
        """
        if d.lower() in ['u']:
            return find_mult(self._n, self._p, self._uk, uv)
        return find_mult(self._m, self._q, self._vk, uv)

    def get_bbox(self):
        """
        Return a bounding box instance for the surface.

        :return: Surface bounding box.
        :rtype: :class:`.BoundingBox`
        """
        bbox = BoundingBox()
        bbox.add_surface(self)
        return bbox

    def check_params(self, u=None, v=None, is_uclosed=False, is_vclosed=False):
        """
        Check that the parameter is within the global domain of the knot vector
        or is within tolerance of a unique knot value. Use
        :class:`.CompareFloats` to compare floats.

        :param float u: Global parameter in u-direction.
        :param float v: Global parameter in v-direction.
        :param bool is_uclosed: Is surface closed in u-direction.
        :param bool is_vclosed: Is surface closed in v-direction.

        :return: Parameters within global domain. Returns same value if
            already within domain. Returns *None* if parameter isn't provided.
        :rtype: tuple
        """
        if u is not None:
            u = check_param(self._n, self._p, self._uk, u, is_uclosed)
        if v is not None:
            v = check_param(self._m, self._q, self._vk, v, is_vclosed)
        return u, v

    def move_seam(self, uv, d, inplace=False, domain='local'):
        """
        Move the seam of a closed surface to a new parameter value.

        :param float uv: New seam parameter.
        :param str d: Direction for seam movement (surface must be closed in
            that direction).
        :param bool inplace: Option to return new surface or modify inpalce.
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: New surface or *True* if successful.
        :rtype: :class:`.NurbsSurface` or bool
        """
        if not self.is_closed(d):
            return False
        if is_local_domain(domain):
            uv = self.local_to_global_param(d, uv)
        uq, vq, qw = move_nurbs_surface_seam(self._n, self._p, self._uk,
                                             self._m, self._q, self._vk,
                                             self._cpw, uv, d)

        if inplace:
            return self.set_cpw(qw) and self.set_knots(uq, vq)

        s = NurbsSurface()
        s.set_cpw(qw)
        s.set_knots(uq, vq)
        s.set_deg(self._p, self._q)
        return s
