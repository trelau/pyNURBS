from __future__ import division

from numpy import array, cross, float64, ndarray, ones, zeros

from pynurbs.geometry.bezier_curve import BezierCurve
from pynurbs.geometry.bounding_box import BoundingBox
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.calculate import surface_area_bezier
from pynurbs.geometry.methods.divide import (extract_bezier_isocurve,
                                             extract_bezier_surface,
                                             split_bezier_surface)
from pynurbs.geometry.methods.evaluate import (bezier_surface_points,
                                               check_param, deCasteljau2,
                                               find_mult_knots,
                                               rat_surface_derivs)
from pynurbs.geometry.methods.geom_utils import (are_points_equal,
                                                 dehomogenize_array2d,
                                                 global_to_local_param,
                                                 homogenize_array2d,
                                                 local_to_global_param)
from pynurbs.geometry.methods.misc import (is_array_like, is_array_type,
                                           is_local_domain)
from pynurbs.geometry.point import Point
from pynurbs.geometry.vector import Vector


class BezierSurface(Geometry):
    """
    Bezier surface.

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
        super(BezierSurface, self).__init__('bezier_surface')
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
        self._u0 = 0.
        self._u1 = 1.
        self._v0 = 0.
        self._v1 = 1.

    def __call__(self, u, v, rtype='Point', domain='local'):
        return self.eval(u, v, rtype, domain)

    @property
    def au(self):
        return self._u0

    @property
    def bu(self):
        return self._u1

    @property
    def av(self):
        return self._v0

    @property
    def bv(self):
        return self._v1

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
        self._u0 = 0.
        self._u1 = 1.
        self._v0 = 0.
        self._v1 = 1.

    def set_domain(self, au=0., bu=1., av=0., bv=1.):
        """
        Set the domain of the Bezier surface between [a, b] where a < b.

        :param float au: Lower bound in u-direction.
        :param float bu: Upper bound in u-direction.
        :param float av: Lower bound in v-direction.
        :param float bv: Upper bound in v-direction.

        :return: *True* if successful, *False* if not.
        :rtype: bool
        """
        if self._cp is not None and au <= bu and av <= bv:
            self._u0 = au
            self._u1 = bu
            self._v0 = av
            self._v1 = bv
            self._uk = zeros(2 * (self._p + 1), dtype=float64)
            self._uk[:self._p + 1] = au
            self._uk[self._p + 1:] = bu
            self._r = self._uk.size - 1
            self._vk = zeros(2 * (self._q + 1), dtype=float64)
            self._vk[:self._q + 1] = av
            self._vk[self._q + 1:] = bv
            self._s = self._vk.size - 1
            return True
        return False

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
        self._p = self._n
        self._m = nm[1] - 1
        self._q = self._m
        self._w = ones((nm[0], nm[1]), dtype=float64)
        self._cpw = homogenize_array2d(self._cp, self._w)
        return self.set_domain()

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

    def local_to_global_param(self, d, *args):
        """
        Convert parameter(s) from local domain 0 <= u,v <= 1 to global domain
        a <= u,v <= b.

        :param str d: Direction ('u' or 'v').
        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        if d.lower() in ['u']:
            return local_to_global_param(self._u0, self._u1, *args)
        else:
            return local_to_global_param(self._v0, self._v1, *args)

    def global_to_local_param(self, d, *args):
        """
        Convert surface parameter(s) from global domain a <= u,v <= b to local
        domain 0 <= u,v <= 1.

        :param str d: Direction ('u' or 'v').
        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        if d.lower() in ['u']:
            return global_to_local_param(self._u0, self._u1, *args)
        else:
            return global_to_local_param(self._v0, self._v1, *args)

    def eval(self, u, v, rtype='Point', domain='local'):
        """
        Evaluate surface at parametric points.

        :param float u: Parametric point in u-direction.
        :param float v: Parametric point in v-direction.
        :param str rtype: Option to return a NumPy array or Point instance
            (rtype = 'Point' or 'ndarray').
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Point on surface.
        :rtype: :class:`.Point` or ndarray
        """
        if not is_local_domain(domain):
            u = self.global_to_local_param('u', u)
            v = self.global_to_local_param('v', v)
        pw = deCasteljau2(self._cpw, self._n, self._m, u, v)
        pnt = pw[:-1] / pw[-1]
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

        if not is_local_domain(domain):
            ulist = [self.global_to_local_param('u', ui) for ui in ulist]
            vlist = [self.global_to_local_param('v', vi) for vi in vlist]
        pnts = bezier_surface_points(self._n, self._m, self._cpw,
                                     ulist, vlist)
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
        # Use NURBS rational curve derivative to compute Bezier derivative.
        # Convert to global since using a NURBS method.
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
            p0 = self.eval(u, v, domain=domain)
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
            p0 = self.eval(u, v, domain=domain)
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
        return surface_area_bezier(s.n, s.m, s.cpw)

    def split(self, u=None, v=None, domain='local'):
        """
        Split Bezier surface.

        :param float u: Location of split (in v-direction).
        :param float v: Location of split (in u-direction).
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Tuple of Bezier surfaces after splitting.
        :rtype: tuple

        .. note::
            If only one parameter is provided, then two surfaces are
            created and the length of the tuple is two (s1, s2). If both
            parameters are provided, four surfaces are created and the length
            of tuple is four (s1, s2, s3, s4).
        """
        if not is_local_domain(domain):
            u = self.global_to_local_param('u', u)
            v = self.global_to_local_param('v', v)
        if u is None or v is None:
            # Split in one direction only.
            qw1, qw2, ab1, ab2 = split_bezier_surface(self._cpw, self._n,
                                                      self._m, u, v,
                                                      self._u0, self._u1,
                                                      self._v0, self._v1)
            s1, s2 = BezierSurface(), BezierSurface()
            s1.set_cpw(qw1)
            s2.set_cpw(qw2)
            s1.set_domain(ab1[0, 0], ab1[0, 1], ab1[1, 0], ab1[1, 1])
            s2.set_domain(ab2[0, 0], ab2[0, 1], ab2[1, 0], ab2[1, 1])
            return s1, s2
        elif u is not None and v is not None:
            # Split in both directions, starting  at u (in v-direction).
            qw1, qw2, ab1, ab2 = split_bezier_surface(self._cpw, self._n,
                                                      self._m, u, None,
                                                      self._u0, self._u1,
                                                      self._v0, self._v1)
            # Split first surface at v (in u-direction).
            qw3, qw4, ab3, ab4 = split_bezier_surface(qw1, self._n, self._m,
                                                      None, v, ab1[0, 0],
                                                      ab1[0, 1], ab1[1, 0],
                                                      ab1[1, 1])
            # Split second surface at v (in u-direction).
            qw5, qw6, ab5, ab6 = split_bezier_surface(qw2, self._n, self._m,
                                                      None, v, ab2[0, 0],
                                                      ab2[0, 1], ab2[1, 0],
                                                      ab2[1, 1])
            # Create surfaces
            s3, s4, s5, s6 = (BezierSurface(), BezierSurface(),
                              BezierSurface(), BezierSurface())
            s3.set_cpw(qw3)
            s4.set_cpw(qw4)
            s5.set_cpw(qw5)
            s6.set_cpw(qw6)
            s3.set_domain(ab3[0, 0], ab3[0, 1], ab3[1, 0], ab3[1, 1])
            s4.set_domain(ab4[0, 0], ab4[0, 1], ab4[1, 0], ab4[1, 1])
            s5.set_domain(ab5[0, 0], ab5[0, 1], ab5[1, 0], ab5[1, 1])
            s6.set_domain(ab6[0, 0], ab6[0, 1], ab6[1, 0], ab6[1, 1])
            return s3, s4, s5, s6

    def isocurve(self, u=None, v=None, domain='local'):
        """
        Extract iso-curve from Bezier surface.

        :param float u: Parameter *u* to extract curve at (in v-direction).
        :param float v: Parameter *v* to extract curve at (in u-direction).
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Bezier curve along parameter.
        :rtype: :class:`.BezierCurve`
        """
        if not is_local_domain(domain):
            u = self.global_to_local_param('u', u)
            v = self.global_to_local_param('v', v)
        qw = extract_bezier_isocurve(self._n, self._m, self._cpw, u, v)
        c = BezierCurve()
        c.set_cpw(qw)
        if u is not None:
            c.set_domain(self._v0, self._v1)
        if v is not None:
            c.set_domain(self._u0, self._u1)
        return c

    def extract(self, u0=0., u1=1., v0=0., v1=1., domain='local'):
        """
        Extract Bezier surface bounded by parameters.

        :param float u0: Starting parameter in u-direction.
        :param float u1: Ending parameter in u-direction.
        :param float v0: Starting parameter in v-direction.
        :param float v1: Ending parameter in v-direction.
        :param str domain: Option to use local (0 <= u, v <= 1) or global
            (a <= u, v <= b) domain ('local', 'l', 'global', 'g').

        :return: Extracted Bezier surface.
        :rtype: :class:`.BezierSurface`
        """
        if not is_local_domain(domain):
            u0 = self.global_to_local_param('u', u0)
            u1 = self.global_to_local_param('u', u1)
            v0 = self.global_to_local_param('v', v0)
            v1 = self.global_to_local_param('v', v1)
        qw, ab = extract_bezier_surface(self._n, self._m, self._cpw, u0, u1,
                                        v0, v1, self._u0, self._u1, self._v0,
                                        self._v1)
        s = BezierSurface()
        s.set_cpw(qw)
        s.set_domain(ab[0, 0], ab[0, 1], ab[1, 0], ab[1, 1])
        return s

    def decompose(self):
        """
        Decompose surface into Bezier patches. For a Bezier surface this simply
        returns a copy of itself.

        :return: Number and list of Bezier patches (nbu, nbv, patches).
        :rtype: tuple
        """
        return 1, 1, [[self.copy()]]

    def get_mult(self, uv, d):
        """
        Get multiplicity of parameter.

        :param float uv: Parameter.
        :param str d: Direction to check ('u' or 'v')

        :return: Multiplicity of parameter.
        :rtype: int
        """
        if d.lower() in ['u']:
            if uv <= self._u0 or uv >= self._u1:
                return self._p + 1
            return 0
        if uv <= self._v0 or uv >= self._v1:
            return self._q + 1
        return 0

    def get_mult_and_knots(self, d):
        """
        Get an array of unique knot values for the curve and their
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
