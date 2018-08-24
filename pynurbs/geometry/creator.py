from numpy import array, float64

from pynurbs.config import Settings
from pynurbs.geometry.bezier_curve import BezierCurve
from pynurbs.geometry.bezier_surface import BezierSurface
from pynurbs.geometry.checker import CheckGeom
from pynurbs.geometry.icurve import ICurve
from pynurbs.geometry.line import Line
from pynurbs.geometry.methods.create import (create_icurve_by_points,
                                             create_system_by_points,
                                             fit_plane,
                                             line_by_points, plane_by_axes,
                                             plane_by_normal, plane_by_points,
                                             planes_by_offset,
                                             points_at_kinks, vector_by_axis,
                                             vector_by_points)
from pynurbs.geometry.methods.interpolate import (global_curve_interpolation,
                                                  interpolate_curves)
from pynurbs.geometry.methods.transform import translate_curve
from pynurbs.geometry.nurbs_curve import NurbsCurve
from pynurbs.geometry.nurbs_surface import NurbsSurface
from pynurbs.geometry.point import Point, Point2D
from pynurbs.geometry.vector import Vector


class CreateGeom(object):
    """
    Geometry creator.
    """

    @staticmethod
    def point2d(uv=(0., 0.)):
        """
        Create a 2-D point.

        :param array_like uv: Location of point.

        :return: Point at *uv*.
        :rtype: :class:`.Point2D`
        """
        if CheckGeom.is_point2d(uv):
            return uv.copy()
        if len(uv) == 2:
            return Point2D(uv)
        return None

    @staticmethod
    def point(xyz=(0., 0., 0.)):
        """
        Create a Point.

        :param array_like xyz: Location of point.

        :return: Point at *xyz*.
        :rtype: :class:`.Point`
        """
        if CheckGeom.is_point(xyz):
            return xyz.copy()
        if len(xyz) == 3:
            return Point(xyz)
        return None

    @staticmethod
    def points_by_array(pnts):
        """
        Create a list of Point instances.

        :param array_like pnts: Array of point locations.

        :return: List of Point instances.
        :rtype: list
        """
        return [CreateGeom.point(pi) for pi in pnts]

    @staticmethod
    def point_from_other(curve, d, u0=None, p0=None, domain='local'):
        """
        Create a point along a curve at a given distance from another point.

        :param curve: Curve to generate point on.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param float d: Distance from given point (positive or negative).
        :param float u0: Parameter of other point.
        :param p0: Other point instead of parameter (will be inverted).
        :type p0: :class:`.Point` or array_like
        :param str domain: Option to use local or global domain.

        :return: Point on curve at given distance from reference point and
            parameter (pnt, param).
        :rtype: tuple
        """
        raise NotImplementedError('No Python implementation.')
        # return curve_point_from_other(curve, d, u0, p0, domain)

    @staticmethod
    def points_along_curve(curve, ds=None, npts=None, u0=0., u1=1.,
                           domain='local', s0=None, s1=None):
        """
        Generate points along a curve.

        :param curve: Curve used to generate points.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param float ds: Point spacing based on arc-length of curve.
        :param int npts: Number of desired points along the curve. This will be
            treated as a minimum number if both *ds and *npts* are provided.
        :param float u0: Starting parameter.
        :param float u1: Ending parameter.
        :param str domain: Domain for *u0* and *u1* ('local' or 'global').
        :param float s0: Option to create first point at specified distance
            from beginning of curve (s0 >= 0).
        :param float s1: Option to create last point at specified distance
            from end of curve (s1 >= 0).

        :return: Class containing point along curve results.
        :rtype: :class:`.PointsAlongCurve`
        """
        raise NotImplementedError('No Python implementation.')
        # results = points_along_curve(curve, ds, npts, u0, u1, 'Point',
        #                              domain, s0, s1)
        #
        # return CreatedPoints(results[1], results[2])

    @staticmethod
    def points_at_kinks(curve, angle=30., u0=0., u1=1., domain='local'):
        """
        Create points at kinks of a curve.

        :param curve:
        :param angle:
        :param u0:
        :param u1:
        :param domain:

        :return:
        """
        if not isinstance(angle, (int, float)):
            return None
        if CheckGeom.is_line(curve):
            results = [(u0, curve.eval(u0)),
                       (u1, curve.eval(u1))]
        else:
            results = points_at_kinks(curve, angle, u0, u1, domain)
        params = [row[0] for row in results]
        pnts = [row[1] for row in results]
        return CreatedPoints(pnts, params)

    @staticmethod
    def vector_by_points(p0, p1, p2=None):
        """
        Create a vector defined by two or three points.

        :param p0: Origin of plane.
        :type p0: :class:`.Point` or array_like
        :param p1: Point defining vector *p1* - *p0*.
        :type p1: :class:`.Point` or array_like
        :param p2: Point defining vector *p2* - *p0*.
        :type p2: :class:`.Point` or array_like

        :return: A vector from *p0* to *p1* if only two points are provided,
            or a vector defined by the cross product of *p10* x *p20* if all
            three points are provided.
        :rtype: :class:`.Vector`
        """
        return vector_by_points(p0, p1, p2)

    @staticmethod
    def vector_by_axis(axis='x', origin=(0., 0., 0.)):
        """
        Create a vector along the specified axis.

        :param str axis: Axis ('x', 'y', or 'z').
        :param array_like origin: Origin of vector.

        :return: Vector along specified axis.
        :rtype: :class:`.Vector`
        """
        return vector_by_axis(axis, origin)

    @staticmethod
    def line_by_vector(p0, v):
        """
        Create a line by an origin and direction vector.

        :param p0:
        :param v:

        :return:
        """
        p0 = CheckGeom.to_point(p0)
        if not CheckGeom.is_vector(v):
            v = array(v, dtype=float64)
            v = Vector(v, p0)
        return Line(p0, v)

    @staticmethod
    def line_by_points(p0, p1):
        """
        Create a line defined by two points.

        :param p0: Origin of plane.
        :type p0: :class:`.Point` or array_like
        :param p1: Point defining vector *p1* - *p0*.
        :type p1: :class:`.Point` or array_like

        :return: A line defined by an oriign at *p0* and a vector *p10*.
        :rtype: :class:`.Line`
        """
        return line_by_points(p0, p1)

    @staticmethod
    def plane_by_normal(p0, vn):
        """
        Create a plane by an origin and normal vector.

        :param p0: Origin of plane.
        :type p0: :class:`.Point` or array_like
        :param vn: Normal vector of plane.
        :type vn: :class:`.Vector` or array_like

        :return: Plane with given origin and normal vector.
        :rtype: :class:`.Plane`
        """
        return plane_by_normal(p0, vn)

    @staticmethod
    def plane_by_points(p0, p1, p2):
        """
        Create a plane defined by three points.

        :param p0: Origin of plane.
        :type p0: :class:`.Point` or array_like
        :param p1: Point defining vector *p1* - *p0*.
        :type p1: :class:`.Point` or array_like
        :param p2: Point defining vector *p2* - *p0*.
        :type p2: :class:`.Point` or array_like

        :return: A plane with a normal vector defined by the cross product of
            *p10* x *p20* and an x-axis oriented towards *p1*
        :rtype: :class:`.Plane`
        """
        return plane_by_points(p0, p1, p2)

    @staticmethod
    def fit_plane(pnts, tol=None):
        """
        Fit a plane to a scattered set of points.

        :param pnts: Points to fit.
        :type pnts: list of :class:`.Point` instances or array_like
        :param float tol: Tolerance for checking the fit. If *None* is
            provided then the plane will be fit to the points. If a float is
            provided then the plane will not be created if the distance from
            any points is greater than *tol*.

        :return: Plane that best fits data.
        :rtype: :class:`.Plane`
        """
        return fit_plane(pnts, tol)

    @staticmethod
    def plane_by_axes(p0=(0., 0., 0.), axes='xz', sys=None):
        """
        Create a plane defined by an origin and standard axes.

        :param p0: Origin of plane.
        :type p0: :class:`.Point` or array_like
        :param axes: Standard axes, one of 'xy', 'xz', or 'yz'.
        :param sys: Reference system for axes.
        :type sys: :class:`.System`

        :return: Plane oriented by axes.
        :rtype: :class:`.Plane`
        """
        return plane_by_axes(p0, axes, sys)

    @staticmethod
    def planes_by_offset(plane, offset=1., n=1):
        """
        Create planes by offsetting an original plane.

        :param plane: Plane to offset.
        :type plane: :class:`.Plane`
        :param float offset: Distance to offset.
        :param int n: Number of planes to generate.

        :return: List of planes offset from original.
        :rtype: list
        """
        return planes_by_offset(plane, offset, n)

    @staticmethod
    def planes_along_curve(curve, ds=None, npts=None, pref=None, u0=0., u1=1.,
                           domain='local', s0=None, s1=None):
        """
        Create planes along a curve.

        :param curve:
        :param ds:
        :param int npts:
        :param pref:
        :param u0:
        :param u1:
        :param domain:
        :param s0:
        :param s1:

        :return: List of :class:`.Plane` instances along the curve.
        :rtype: list
        """
        raise NotImplementedError('No Python implementation.')
        # return planes_along_curve(curve, ds, npts, pref, u0, u1, domain, s0,
        #                           s1)

    @staticmethod
    def planes_between_planes(plane0, plane1, ds):
        """
        Create planes between two planes.

        :param plane0:
        :param plane1:
        :param ds:

        :return:
        """
        raise NotImplementedError('No Python implementation.')
        # if not CheckGeom.is_plane(plane0) or not CheckGeom.is_plane(plane1):
        #     return []
        # return planes_between_planes(plane0, plane1, ds)

    @staticmethod
    def bezier_curve(cp=None, cpw=None):
        """
        Create a Bezier curve.

        :param array_like cp: Control points of curve.
        :param array_like cpw: Homogeneous control points of curve.

        :return: Bezier curve.
        :rtype: :class:`.BezierCurve`
        """
        c = BezierCurve()
        if cp is not None:
            c.set_cp(cp)
            return c
        if cpw is not None:
            c.set_cpw(cpw)
            return c
        return c

    @staticmethod
    def bezier_surface(cp=None, cpw=None):
        """
        Create a Bezier surface.

        :param array_like cp: Control points of surface.
        :param array_like cpw: Homogeneous control points of surface.

        :return: Bezier surface.
        :rtype: :class:`.BezierSurface`
        """
        s = BezierSurface()
        if cp is not None:
            s.set_cp(cp)
            return s
        if cpw is not None:
            s.set_cpw(cpw)
            return s
        return s

    @staticmethod
    def nurbs_curve():
        """
        Create a NURBS curve.

        :return: NURBS curve. Properties can be set after object creation.
        :rtype: :class:`.NurbsCurve`
        """
        return NurbsCurve

    @staticmethod
    def nurbs_surface():
        """
        Create a NURBS surface.

        :return: NURBS surface. Properties can be set after object creation.
        :rtype: :class:`.NurbsSurface`
        """
        return NurbsSurface()

    @staticmethod
    def interpolate_points(pnts, p=3, method='chord'):
        """
        Create a NURBS curve by interpolating points.

        :param array_like pnts: Array of points.
        :param int p: Degree for interpolated curve.
        :param str method: Method for determining curve parameters
            ('uniform', 'chord', or 'centripetal').

        :return: Curve interpolating points.
        :rtype: :class:`.NurbsCurve`

        .. warning::
            The method does not automatically check for coincident points
            which will cause a fatal error.
        """
        qp = CreateGeom.points_by_array(pnts)
        cp, uk, p = global_curve_interpolation(qp, p, method)
        c = NurbsCurve()
        c.set_cp(cp)
        c.set_knots(uk)
        c.set_deg(p)
        return c

    @staticmethod
    def interpolate_curves(crvs, q=3, method='chord', auto_reverse=True):
        """
        Create a NURBS surface by interpolating curves.

        :param list crvs: List of curves to interpolate.
        :param int q: Degree for interpolated surface in direction of curves
            (v-direction).
        :param str method: Method for determining curve parameters
            ('uniform', 'chord', or 'centripetal').
        :param bool auto_reverse: Option to check direction of curves and
            reverse if necessary.

        :return: Surface interpolating curves. The  degree of the surface in
            the v-direction will be specified by *p*, while the degree in
            the u-direction will be the highest degree of any curve (i.e.,
            the curve degrees are elevated if necessary).
        :rtype: :class:`.NurbsSurface`
        """
        # Convert to NURBS if necessary.
        _crvs = []
        for c in crvs:
            if isinstance(c, NurbsCurve):
                _crvs.append(c)
            elif isinstance(c, BezierCurve):
                ci = NurbsCurve()
                ci.set_cpw(c.cpw)
                ci.set_knots(c.uk)
                ci.set_deg(c.p)
                _crvs.append(ci)
            elif isinstance(c, ICurve):
                _crvs.append(c.crv3d)

        cpw, uk, vk, p, q = interpolate_curves(_crvs, q, method, False,
                                               auto_reverse)
        s = NurbsSurface()
        s.set_cpw(cpw)
        s.set_knots(uk, vk)
        s.set_deg(p, q)
        return s

    @staticmethod
    def surface_by_curve_drag(curve, vector, vmax=1., vmin=-1.):
        """
        Create a surface by dragging a curve along a vector in both direction.

        :param curve: Curve to drag.
        :type curve: :class`.BezierCurve`, :class`.NurbsCurve`,
            or :class`.ICurve`
        :param array_like vector: Drag direction vector.
        :param float vmax: Distance to drag in positive v-direction.
        :param float vmin: Distance to drag in negative v-direction.

        :return: Surface created from dragging curves.
        :rtype: :class:`.NurbsSurface`
        """
        tol = Settings.gtol
        if not CheckGeom.is_curve(curve):
            return None
        if CheckGeom.is_icurve(curve):
            curve = curve.crv3d
        if vmax < 0.:
            vmax = 0.
        if vmin > 0.:
            vmin = 0.
        if abs(vmax) <= tol and abs(vmin) <= tol:
            return None
        c0 = translate_curve(curve, vector, vmin, False)
        c1 = translate_curve(curve, vector, vmax, False)
        return CreateGeom.interpolate_curves([c0, c1])

    @staticmethod
    def surface_from_plane(plane, w=1., h=1.):
        """
        Create a surface from a plane by specifying surface height and width.

        :param plane:
        :param w:
        :param h:

        :return:
        """
        if not CheckGeom.is_plane(plane):
            return None
        if w <= 0. or h <= 0.:
            return None
        w2 = w / 2.
        h2 = h / 2.
        p0 = plane.eval(-w2, 0.)
        p1 = plane.eval(w2, 0.)
        c = CreateGeom.interpolate_points([p0, p1], p=1)
        return CreateGeom.surface_by_curve_drag(c, plane.vv.ijk, h2, -h2)

    @staticmethod
    def icurve_by_points(surface, p0, p1, isurf=None, trim=True):
        """
        Create an intersection curve between two points on the surface.

        :param surface:
        :type surface: :class:`.BezierSurface` or :class`.NurbsSurface`
        :param point_like p0: Starting point.
        :param point_like p1: Ending point.
        :param surface_like isurf: Intersection surface. If *None* if
            provided then a plane will be created between the two points using
            the surface normal at *p0*.
        :param bool trim: Option to trim the curve or not.

        :return: Intersection curve between *p0* and *p1*. Returns *None* if
            method fails.
        :rtype: :class:`.ICurve`
        """
        if not CheckGeom.is_surface(surface):
            return None
        if not CheckGeom.is_point_like(p0) or not CheckGeom.is_point_like(p1):
            return None
        return create_icurve_by_points(surface, p0, p1, isurf, trim)

    @staticmethod
    def system_by_points(origin=(0., 0., 0.), x_axis=(1., 0., 0.),
                         xz_plane=(0., 0., 1.)):
        """
        Create a coordinate system by three points.

        :param origin:
        :param x_axis:
        :param xz_plane:

        :return:
        :rtype:
        """
        return create_system_by_points(origin, x_axis, xz_plane)


class CreatedPoints(object):
    """
    Class to handle results of point creation methods.
    """

    def __init__(self, pnts, params):
        self._pnts = pnts
        self._params = params

    @property
    def npts(self):
        return len(self._pnts)

    @property
    def pnts(self):
        return self._pnts

    @property
    def pnt(self):
        try:
            return self._pnts[0]
        except IndexError:
            return None

    @property
    def params(self):
        return self._params

    @property
    def param(self):
        try:
            return self._params[0]
        except IndexError:
            return None

    @property
    def zip_results(self):
        return zip(self._pnts, self._params)
