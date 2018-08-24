from numpy import array, cross, dot, float64

from pynurbs.config import Settings
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.evaluate import check_param
from pynurbs.geometry.methods.geom_utils import (global_to_local_param,
                                                 local_to_global_param)
from pynurbs.geometry.methods.intersect_surface import (refine_spi_point,
                                                        refine_ssi_point)
from pynurbs.geometry.methods.misc import is_array_type, is_local_domain
from pynurbs.geometry.plane import Plane
from pynurbs.geometry.point import Point, Point2D
from pynurbs.geometry.vector import Vector


class ICurve(Geometry):
    """
    Intersection curve.
    """

    def __init__(self, s1, s2, crv3d, crv2d_s1, crv2d_s2, ftol):
        super(ICurve, self).__init__('icurve')
        self._s1 = s1
        self._s2 = s2
        self._crv3d = crv3d
        self._crv2d_s1 = crv2d_s1
        self._crv2d_s2 = crv2d_s2
        self._ftol = ftol
        self._crv3d.set_color(*self.color)
        # self._check_order()

    def __call__(self, u, rtype='Point', domain='local', tol=None):
        return self.eval(u, rtype, domain, tol)

    @property
    def s1(self):
        return self._s1

    @property
    def s2(self):
        return self._s2

    @property
    def crv3d(self):
        return self._crv3d

    @property
    def a(self):
        return self._crv3d.a

    @property
    def b(self):
        return self._crv3d.b

    @property
    def p(self):
        return self._crv3d.p

    @property
    def m(self):
        return self._crv3d.m

    @property
    def n(self):
        return self._crv3d.n

    @property
    def uk(self):
        return self._crv3d.uk

    @property
    def knots(self):
        return self.get_mult_and_knots()[-1]

    @property
    def mult(self):
        return self.get_mult_and_knots()[1]

    @property
    def cp(self):
        return self._crv3d.cp

    @property
    def w(self):
        return self._crv3d.w

    @property
    def cpw(self):
        return self._crv3d.cpw

    @property
    def data(self):
        return self.n, self.p, self.uk, self.cpw

    @property
    def ftol(self):
        return self._ftol

    @property
    def is_spi(self):
        return isinstance(self._s2, Plane)

    @property
    def is_ssi(self):
        return not isinstance(self._s2, Plane)

    @property
    def is_closed(self):
        return self._crv3d.is_closed

    @property
    def length(self):
        return self.arc_length()

    @property
    def cp2d_s1(self):
        return self._crv2d_s1.cp

    @property
    def cp2d_s2(self):
        return self._crv2d_s2.cp

    @property
    def occ(self):
        return self.get_occ()

    @property
    def handle(self):
        return self._crv3d.handle

    @property
    def adaptor(self):
        return self._crv3d.adaptor

    def get_occ(self):
        """
        Get OCC object.

        :return:
        """
        return self._crv3d.get_occ()

    def _check_order(self):
        """
        Make sure that the derivative of the intersection curve is in the same
        direction as the 3-D curve. If not, reverse the curves.
        """
        v3d = self._crv3d.deriv(0.5, 1, rtype='ndarray')
        du = self.deriv(0.5, 1, rtype='ndarray')
        if dot(v3d, du) < 0.:
            self._crv3d.reverse(True)
            self._crv2d_s1.reverse(True)
            self._crv2d_s2.reverse(True)

    def local_to_global_param(self, *args):
        """
        Convert parameter(s) from local domain 0 <= u <= 1 to global domain
        a <= u <= b.

        :param args: Local parameter(s).

        :return: Global parameter(s).
        """
        return local_to_global_param(self.a, self.b, *args)

    def global_to_local_param(self, *args):
        """
        Convert parameter(s) from global domain a <= u <= b to local domain
        0 <= u <= 1.

        :param args: Global parameter(s).

        :return: Local parameter(s).
        """
        return global_to_local_param(self.a, self.b, *args)

    def has_surf(self, sref):
        """
        Check to see if the surface is used in the intersection curve.

        :param surface_like sref: Reference surface.

        :return: *True* if used, *False* if not.
        :rtype: bool
        """
        if sref is self._s1 or sref is self._s2:
            return True
        return False

    def get_curve2d(self, sref):
        """
        Get the 2-D curve associated with the reference surface.

        :param sref: Reference surface (must have been used in the
            intersection).

        :return: 2-D intersection curve.
        :rtype: :class:`.NurbsCurve`
        """
        if sref is self._s1:
            return self._crv2d_s1
        if sref is self._s2:
            return self._crv2d_s2
        return None

    def eval2d(self, u, rtype='Point', domain='local', tol=None, sref=None):
        """
        Evaluate the intersection curve in 2-D space for each surface.

        :param float u: Parametric point.
        :param str rtype: Option to return a NumPy array or Point2D instance
            (rtype = 'Point' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').
        :param float tol: Tolerance for point refinement.
        :param surface_like sref: Option to provide one of the two surfaces
            used in the intersection curve. If present, the method will
            only return the 2-D parameters associated to that surface.

        :return: Point on curve. Will return *None* if *sref* is not in the
            intersection.
        :rtype: :class:`.Point2D` or ndarray
        """
        if is_local_domain(domain):
            u = self.local_to_global_param(u)

        # Evaluate 2-D curves.
        s1, s2 = self._s1, self._s2
        u1, v1 = self._crv2d_s1.eval(u, rtype='ndarray', domain='global')[:-1]
        u2, v2 = self._crv2d_s2.eval(u, rtype='ndarray', domain='global')[:-1]

        # Project 3-D point to surfaces to get initial parameters.
        # s1, s2 = self._s1, self._s2
        # p3d = self._crv3d.eval(u, domain='global')
        # u1, v1 = invert_point_on_surface(p3d, s1)
        # if self.is_spi:
        #     u2, v2 = invert_point_on_plane(p3d, s2)
        # else:
        #     u2, v2 = invert_point_on_surface(p3d, s2)

        # Refine parameters.
        if tol is None:
            tol = Settings.gtol / 100.
        if self.is_spi:
            u1, v1, u2, v2 = refine_spi_point(s1, s2, u1, v1, tol)[:-1]
        else:
            u1, v1, u2, v2 = refine_ssi_point(s1, s2, u1, v1, u2, v2, tol)[:-1]

        if is_array_type(rtype):
            if sref is s1:
                return array([u1, v1], dtype=float64)
            if sref is s2:
                return array([u2, v2], dtype=float64)
            return (array([u1, v1], dtype=float64),
                    array([u2, v2], dtype=float64))
        else:
            if sref is s1:
                return Point2D((u1, v1))
            if sref is s2:
                return Point2D((u2, v2))
            return Point2D((u1, v1)), Point2D((u2, v2))

    def eval(self, u, rtype='Point', domain='local', tol=None):
        """
        Evaluate curve at parametric point.

        :param float u: Parametric point.
        :param str rtype: Option to return a NumPy array or Point instance
            (rtype = 'Point' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').
        :param float tol: Tolerance for point refinement.

        :return: Point on curve.
        :rtype: :class:`.Point` or ndarray
        """
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        uv1, uv2 = self.eval2d(u, rtype='ndarray', domain='global', tol=tol)
        p3d = self._s1.eval(uv1[0], uv1[1], rtype='ndarray', domain='global')
        if is_array_type(rtype):
            return p3d
        else:
            return Point(p3d)

    def deriv(self, u, k=1, d=1, rtype='Vector', domain='local', tol=None):
        """
        Compute the derivative of the intersection curve. Only supports
        first derivates.

        :param float u: Parametric point.
        :param int k: Derivative to return (0 <= k <= 1).
        :param int d: Highest derivative to compute. Currently only supports
            first derivative.
        :param str rtype: Option to return a NumPy array or a Vector instance
            (rtype = 'Vector' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').
        :param float tol: Tolerance for point refinement.

        :return: First derivative or intersection curve.
        :rtype: :class:`.Vector` or ndarray
        """
        if is_local_domain(domain):
            u = self.local_to_global_param(u)
        # Evaluate 2-D points.
        s1, s2 = self._s1, self._s2
        uv1, uv2 = self.eval2d(u, rtype='ndarray', domain='global', tol=tol)
        # Evaluate surface normals.
        vn1 = s1.norm(uv1[0], uv1[1], rtype='ndarray', domain='global')
        if self.is_spi:
            vn2 = s2.vn.ijk
        else:
            vn2 = s2.norm(uv2[0], uv2[1], rtype='ndarray', domain='global')
        # First derivative is cross product
        du = cross(vn1, vn2)
        if is_array_type(rtype):
            return du
        else:
            p0 = self.eval(u, domain='global')
            return Vector(du, p0)

    def extract(self, u0, u1, domain='local'):
        """
        Extract a curve.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').

        :return: Curve between *u0* and *u1*.
        :rtype: :class:`.ICurve`
        """
        if is_local_domain(domain):
            u0, u1 = self.local_to_global_param(u0, u1)
        if u0 > u1:
            u0, u1 = u1, u0

        # Evaluate points on intersection curve.
        p3d_0 = self.eval(u0, rtype='ndarray', domain='global')
        p3d_1 = self.eval(u1, rtype='ndarray', domain='global')
        uv0_s1_0, uv0_s2_0 = self.eval2d(u0, rtype='ndarray', domain='global')
        uv0_s1_1, uv0_s2_1 = self.eval2d(u1, rtype='ndarray', domain='global')
        # Extract curves
        crv3d = self._crv3d.extract(u0, u1, domain='global')
        crv2d_s1 = self._crv2d_s1.extract(u0, u1, domain='global')
        crv2d_s2 = self._crv2d_s2.extract(u0, u1, domain='global')
        # Force endpoints of curves to match original points by adjusting
        # first and last control points.
        crv3d.modify_cp(0, p3d_0)
        crv3d.modify_cp(-1, p3d_1)
        crv2d_s1.modify_cp(0, uv0_s1_0)
        crv2d_s1.modify_cp(-1, uv0_s1_1)
        crv2d_s2.modify_cp(0, uv0_s2_0)
        crv2d_s2.modify_cp(-1, uv0_s2_1)
        return ICurve(self._s1, self._s2, crv3d, crv2d_s1, crv2d_s2,
                      self._ftol)

    def check_param(self, u):
        """
        Check that the parameter is within the global domain of the knot vector
        or is within tolerance of a unique knot value. Use
        :class:`.CompareFloats` to compare floats.

        :param float u: Global parameter.

        :return: Parameter within global domain or near interior knot value.
        :rtype: float
        """
        c3d = self._crv3d
        return check_param(c3d.n, c3d.p, c3d.uk, u)

    def reverse(self, inplace=True):
        """
        Reverse direction of curve.

        :param bool inplace: Option to return new curve or modify
            existing curve.

        :return: Reversed NURBS curve.
        :rtype: :class:`.ICurve`
        """
        if inplace:
            self._crv3d.reverse()
            self._crv2d_s1.reverse()
            self._crv2d_s2.reverse()
            return True
        crv3d = self._crv3d.reverse(False)
        crv2d_s1 = self._crv2d_s1.reverse(False)
        crv2d_s2 = self._crv2d_s1.reverse(False)
        return ICurve(self._s1, self._s2, crv3d, crv2d_s1, crv2d_s2,
                      self._ftol)

    def arc_length(self, u0=0., u1=1., domain='local'):
        """
        Estimate the arc length between curve parameters.

        :param float u0: Starting parameter.
        :param float u1: Ending parameter.
        :param str domain: Option to use local or global domain.

        :return: Arc length of curve between *u0* and *u1*.
        :rtype: float

        ..note:
        Arc length estimates will be a function of the tolerance used in the
        original intersection method.
        """
        return self._crv3d.arc_length(u0, u1, domain)

    def get_mult_and_knots(self):
        """
        Get an array of unique knots values for the curve and their
        multiplicities.

        :return: Number of unique knots, multiplicities, and the knots
            (nu, um, uq).
        :rtype: tuple
        """
        return self._crv3d.get_mult_and_knots()
