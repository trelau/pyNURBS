from itertools import combinations

from numpy import float64, zeros

from pynurbs.config import Settings
from pynurbs.geometry.checker import CheckGeom
from pynurbs.geometry.icurve import ICurve
from pynurbs.geometry.methods.check import si_curve_nearest_point
from pynurbs.geometry.methods.interpolate import global_curve_interpolation
from pynurbs.geometry.methods.intersect_curve import (intersect_curve_curve,
                                                      intersect_curve_icurve,
                                                      intersect_curve_line,
                                                      intersect_curve_plane,
                                                      intersect_curve_surface,
                                                      intersect_icurve_icurve,
                                                      intersect_icurve_line,
                                                      intersect_icurve_plane,
                                                      intersect_icurve_surface)
from pynurbs.geometry.methods.intersect_line import (intersect_line_line,
                                                     intersect_line_plane,
                                                     intersect_line_surface)
from pynurbs.geometry.methods.intersect_plane import intersect_plane_plane
from pynurbs.geometry.methods.intersect_surface import (
    intersect_surface_plane,
    intersect_surface_surface, refine_spi_point, refine_ssi_point)
from pynurbs.geometry.methods.invert import invert_points_on_plane
from pynurbs.geometry.nurbs_curve import NurbsCurve
from pynurbs.geometry.point import Point, Point2D


class IntersectGeom(object):
    """
    Intersect geometry.
    """

    @staticmethod
    def perform(geom1, geom2, itol=None, ftol=None, t0=None, t1=None):
        """
        Perform the intersection of two geometries.

        :param geom1: Geometry entity 1.
        :param geom2: Geometry entity 2.
        :param float itol: Intersection tolerance.
        :param float ftol: Surface flatness tolerance.
        :param float t0: Lower parameter on line to exclude results.
        :param float t1: Upper parameter on line to exclude results.

        :return: Intersector object depending on type of *geom1* and *geom2*.
            *None* if returned if entities are not supported.
        """
        # Curve intersection with...
        if CheckGeom.is_curve(geom1):
            # Curve
            if CheckGeom.is_curve(geom2):
                return IntersectCurveCurve(geom1, geom2, itol)
            # Plane
            if CheckGeom.is_plane(geom2):
                return IntersectCurvePlane(geom1, geom2, itol)
            # Surface
            if CheckGeom.is_surface(geom2):
                return IntersectCurveSurface(geom1, geom2, itol)
            # Line
            if CheckGeom.is_line(geom2):
                return IntersectCurveLine(geom1, geom2, itol, t0, t1)

        # Surface intersection with...
        if CheckGeom.is_surface(geom1):
            # Curve
            if CheckGeom.is_curve_like(geom2):
                return IntersectCurveSurface(geom2, geom1, itol, t0, t1)
            # Plane
            if CheckGeom.is_plane(geom2):
                return IntersectSurfacePlane(geom1, geom2, ftol)
            # Surface
            if CheckGeom.is_surface(geom2):
                return IntersectSurfaceSurface(geom1, geom2, ftol)

        # Plane intersection with...
        if CheckGeom.is_plane(geom1):
            # Curve
            if CheckGeom.is_curve_like(geom2):
                return IntersectCurvePlane(geom2, geom1, itol, t0, t1)
            # Surface
            if CheckGeom.is_surface(geom2):
                return IntersectSurfacePlane(geom2, geom1, ftol)
            # Plane
            if CheckGeom.is_plane(geom2):
                return IntersectPlanePlane(geom1, geom2)

        # Line intersection with...
        if CheckGeom.is_line(geom1):
            # Curve
            if CheckGeom.is_curve(geom2):
                return IntersectCurveLine(geom2, geom1, itol, t0, t1)
            # Line
            if CheckGeom.is_line(geom2):
                return IntersectLineLine(geom1, geom2, itol)
            # Plane
            if CheckGeom.is_plane(geom2):
                return IntersectCurvePlane(geom1, geom2, itol, t0, t1)
            # Surface
            if CheckGeom.is_surface(geom2):
                return IntersectCurveSurface(geom1, geom2, itol, t0, t1)

        # Return error if combination not supported.
        return IntersectError()

    @staticmethod
    def curves(curves, itol=None):
        """
        Intersect a list of curves.

        :param curves:
        :param float itol: Intersection tolerance.

        :return: Curve intersection results.
        :type: :class:`.IntersectCurves`
        """
        return IntersectCurves(curves, itol)

    @staticmethod
    def refine_spi(pnt, u0, v0, surface, plane, tol=None):
        """
        Refine a point on a surface-plane intersection.

        :param pnt:
        :param u0:
        :param v0:
        :param surface:
        :param plane:
        :param float tol:

        :return:
        """
        if not CheckGeom.is_point(pnt) or not CheckGeom.is_surface(surface) \
                or not CheckGeom.is_plane(plane):
            return None, None, None, None

        # Refine and update point.
        if tol is None:
            tol = Settings.gtol / 100.
        u0, v0, up, vp, pi = refine_spi_point(surface, plane, u0, v0, tol)
        pnt.set_xyz(pi)
        return u0, v0, up, vp

    @staticmethod
    def refine_ssi(pnt, u1, v1, surface1, u2, v2, surface2, tol=None):
        """
        Refine a point on a surface-surface intersection.

        :param pnt:
        :param u1:
        :param v1:
        :param surface1:
        :param u2:
        :param v2:
        :param surface2:
        :param float tol:

        :return:
        """
        if not CheckGeom.is_point(pnt) or not CheckGeom.is_surface(surface1) \
                or not CheckGeom.is_surface(surface2):
            return None, None, None, None

        # Refine and update point.
        if tol is None:
            tol = Settings.gtol / 100.
        u1, v1, u2, v2, pi = refine_ssi_point(surface1, surface2, u1, v1,
                                              u2, v2, tol)
        pnt.set_xyz(pi)
        return u1, v1, u2, v2


class IntersectError(object):
    """
    Class for handling intersection errors.
    """

    @property
    def success(self):
        return False

    @property
    def npts(self):
        return 0

    @property
    def ncrvs(self):
        return 0


class CurveIntersector(object):
    """
    Base class for handling curve intersection methods and results.
    """

    def __init__(self, c1, c2):
        self._c1 = c1
        self._c2 = c2
        self._npts = 0
        self._results = []
        self._kdt = None

    def _set_results(self, npts, results):
        """
        Set curve intersection results.
        """
        if npts > 0:
            self._npts = npts
            # Replace ndarrays with point instances and build kd-tree.
            data = zeros((npts, 3), dtype=float64)
            for i in range(npts):
                data[i] = results[i][1]
                results[i][1] = Point(results[i][1])
            self._results = results

    @property
    def npts(self):
        return self._npts

    @property
    def success(self):
        if self._npts > 0:
            return True
        return False

    @property
    def points(self):
        if self._npts <= 0:
            return []
        return [results[1] for results in self._results]

    @property
    def parameters(self):
        return [results[0] for results in self._results]

    def point(self, indx=0):
        """
        Return the point result by index.

        :param int indx: Index for point selection.

        :return: Intersection point at index.
        :rtype: :class:`.Point`
        """
        if indx > self._npts - 1:
            return self._results[-1][1]
        return self._results[indx][1]

    def params_by_cref(self, cref):
        """
        Get the parameters of the intersection results by curve reference.

        :param cref: Reference curve (must have been using in the
            intersection method).

        :return: List of parameters for the reference curve. Returns empty
            list if  reference curve is not in the intersection.
        :rtype: list
        """
        if self._c1 is cref:
            return [results[0][0] for results in self._results]
        elif self._c2 is cref:
            if CheckGeom.is_surface_like(cref):
                return [results[0][1:] for results in self._results]
            return [results[0][1] for results in self._results]
        return []


class IntersectLineLine(CurveIntersector):
    """
    Line-line intersection.
    """

    def __init__(self, line1, line2, itol=None):
        super(IntersectLineLine, self).__init__(line1, line2)
        if CheckGeom.is_line(line1) and CheckGeom.is_line(line2):
            self._perform(line1, line2, itol)

    def _perform(self, line1, line2, itol):
        """
        Perform the line-curve intersection.
        """
        _, npts, results = intersect_line_line(line1, line2, itol)
        self._set_results(npts, results)


class IntersectCurveLine(CurveIntersector):
    """
    Curve-line intersection.

    :param line: Line to intersect.
    :type line: :class:`.Line`
    :param curve: Curve to intersect.
    :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param float t0: Lower parameter on line to exclude results.
    :param float t1: Upper parameter on line to exclude results.
    :param float itol: Intersection tolerance. If *None* is provided then
        *itol* from the default settings will be used.
    :param float t0: Lower parameter on line to exclude results.
    :param float t1: Upper parameter on line to exclude results.

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all intersection points in the solution set as
        :class:`.Point` instances.
    :var list parameters: List of parameters on the line and curve after
        intersection as a list of tuples [(t1, u1), (t2, u2), ...].
    """

    def __init__(self, curve, line, itol=None, t0=None, t1=None):
        super(IntersectCurveLine, self).__init__(curve, line)
        if CheckGeom.is_line(line) and CheckGeom.is_curve(curve):
            self._perform(curve, line, itol, t0, t1)

    def _perform(self, curve, line, itol, t0, t1):
        """
        Perform the line-curve intersection.
        """
        if CheckGeom.is_icurve(curve):
            _, npts, results = intersect_icurve_line(curve, line, itol, t0, t1)
        else:
            _, npts, results = intersect_curve_line(curve, line, itol, t0, t1)
        self._set_results(npts, results)


class IntersectCurveCurve(CurveIntersector):
    """
    Curve-curve intersection.

    :param curve1: Curve 1 to intersect.
    :type curve1: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param curve2: Curve 2 to intersect.
    :type curve2: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param float itol: Intersection tolerance. If *None* is provided then
        *itol* from the default settings will be used.

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all intersection points in the solution set as
        :class:`.Point` instances.
    :var list parameters: List of parameters on each curve after intersection
        as a list of tuples [(u11, u21), (u12, u22), ...].
    """

    def __init__(self, curve1, curve2, itol=None):
        super(IntersectCurveCurve, self).__init__(curve1, curve2)
        if CheckGeom.is_curve(curve1) and CheckGeom.is_curve(curve2):
            self._perform(curve1, curve2, itol)

    def _perform(self, curve1, curve2, itol):
        """
        Perform the curve-curve intersection.
        """
        if CheckGeom.is_icurve(curve1) and CheckGeom.is_icurve(curve2):
            _, npts, results = intersect_icurve_icurve(curve1, curve2, itol)
        elif CheckGeom.is_icurve(curve2):
            _, npts, results = intersect_curve_icurve(curve1, curve2, itol)
        elif CheckGeom.is_icurve(curve1):
            _, npts, results = intersect_curve_icurve(curve2, curve1, itol)
            self._c1, self._c2 = self._c2, self._c1
        else:
            _, npts, results = intersect_curve_curve(curve1, curve2, itol)
        self._set_results(npts, results)


class IntersectCurvePlane(CurveIntersector):
    """
    Curve-plane intersection.

    :param curve: Curve to intersect.
    :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param plane: Intersection plane.
    :type plane: :class:`.Plane`
    :param float itol: Intersection tolerance. If *None* is provided then
        *itol* from the default settings will be used.
    :param float t0: Lower parameter on line to exclude results.
    :param float t1: Upper parameter on line to exclude results.

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all intersection points in the solution set as
        :class:`.Point` instances.
    :var list parameters: List of parameters on curve at intersection points.
        [u1, u2, ...]
    """

    def __init__(self, curve, plane, itol=None, t0=None, t1=None):
        super(IntersectCurvePlane, self).__init__(curve, plane)
        if CheckGeom.is_curve_like(curve) and CheckGeom.is_plane(plane):
            self._perform(curve, plane, itol, t0, t1)

    def _perform(self, curve, plane, itol, t0, t1):
        """
        Perform the curve-plane intersection.
        """
        if CheckGeom.is_icurve(curve):
            _, npts, _results = intersect_icurve_plane(curve, plane, itol)
        elif CheckGeom.is_line(curve):
            _, npts, _results = intersect_line_plane(curve, plane, itol, t0,
                                                     t1)
        else:
            _, npts, _results = intersect_curve_plane(curve, plane, itol)
        # Generate plane parameters for each point.
        results = []
        for i in range(0, npts):
            t = _results[i][0]
            pi = _results[i][1]
            u, v = invert_points_on_plane([pi], plane)[0]
            results.append([(t, u, v), pi])
        self._set_results(npts, results)

    @property
    def curve_parameters(self):
        return [results[0][0] for results in self._results]

    @property
    def plane_parameters(self):
        return [results[0][1:] for results in self._results]


class IntersectCurveSurface(CurveIntersector):
    """
    Curve-surface intersection.

    :param curve: Curve to intersect.
    :type curve: :class:`.BezierCurve`, :class:`.NurbsCurve`, or :class:`.Line`
    :param surface: Surface to intersect.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float itol: Intersection tolerance. If *None* is provided then
        *itol* from the default settings will be used.
    :param float t0: Lower parameter on line to exclude results.
    :param float t1: Upper parameter on line to exclude results.

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all intersection points in the solution set as
        :class:`.Point` instances.
    :var list curve_parameters: List of parameters on curve at intersection
        points [u1, u2, ...].
    :var list surface_parameters: List of parameters on surface at intersection
        points [(u1, v1), (u2, v2), ...].
    """

    def __init__(self, curve, surface, itol=None, t0=None, t1=None):
        super(IntersectCurveSurface, self).__init__(curve, surface)
        if CheckGeom.is_curve_like(curve) and CheckGeom.is_surface(surface):
            self._perform(curve, surface, itol, t0, t1)

    def _perform(self, curve, surface, itol, t0, t1):
        """
        Perform the curve-surface intersection.
        """
        if CheckGeom.is_icurve(curve):
            _, npts, results = intersect_icurve_surface(curve, surface, itol)
        elif CheckGeom.is_line(curve):
            _, npts, results = intersect_line_surface(curve, surface, itol,
                                                      t0, t1)
        else:
            _, npts, results = intersect_curve_surface(curve, surface, itol)
        self._set_results(npts, results)

    @property
    def curve_parameters(self):
        return [results[0][0] for results in self._results]

    @property
    def surface_parameters(self):
        return [results[0][1:] for results in self._results]


class IntersectCurves(object):
    """
    Intersect a list of curves.
    """

    def __init__(self, curves, itol=None):
        self._crvs = []
        self._ncrvs = 0
        self._pnts = {}
        self._params = {}
        for c in curves:
            if CheckGeom.is_curve_like(c):
                self._crvs.append(c)
                self._ncrvs += 1
        self._perform(itol)

    def _perform(self, itol):
        """
        Perform intersection for each curve.
        """
        # Initialize results.
        for c in self._crvs:
            self._pnts[c] = []
            self._params[c] = []

        for c1, c2 in combinations(self._crvs, 2):
            cci = IntersectGeom.perform(c1, c2, itol)
            if cci.success:
                # Add points to curve.
                for pi in cci.points:
                    self._pnts[c1].append(pi)
                    self._pnts[c2].append(pi)
                # Add parameters to curve.
                for ui in cci.params_by_cref(c1):
                    self._params[c1].append(ui)
                for ui in cci.params_by_cref(c2):
                    self._params[c2].append(ui)

    @property
    def npts(self):
        npts = 0
        for c in self._crvs:
            npts += self.npts_by_curve(c)
        return npts

    @property
    def success(self):
        if self.npts > 0:
            return True
        return False

    def npts_by_curve(self, curve):
        """
        Get number of intersection points for a given curve.

        :param curve: Reference curve.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`

        :return: Number of intersection points.
        :rtype: int
        """
        if curve not in self._pnts:
            return 0
        return len(self._pnts[curve])

    def points_by_curve(self, curve, sort=False):
        """
        Get all intersection points for a given curve.

        :param curve: Reference curve.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param bool sort: Option to sort points by increasing curve parameters.

        :return: List of intersection points along curve.
        :rtype: list
        """
        if curve not in self._pnts:
            return []
        if not sort:
            return self._pnts[curve]
        data = self.data_by_curve(curve, True)
        return [p[1] for p in data]

    def params_by_curve(self, curve, sort=False):
        """
        Get all intersection parameters for a given curve.

        :param curve: Reference curve.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param bool sort: Option to sort by increasing curve parameters.

        :return: List of intersection parameters along curve.
        :rtype: list
        """
        if curve not in self._pnts:
            return []
        if not sort:
            return self._params[curve]
        data = self.data_by_curve(curve, True)
        return [u[0] for u in data]

    def data_by_curve(self, curve, sort=False):
        """
        Return a list of sorted data (parameters and points) by for a given
        curve.

        :param curve: Reference curve.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param bool sort: Option to sort data by increasing curve parameters.

        :return: List of tuples containing parameters and points sorted by
            parameters.
        :rtype: list
        """
        if curve not in self._pnts or curve not in self._params:
            return []
        prms, pnts = self._params[curve], self._pnts[curve]
        data = zip(prms, pnts)
        if not sort:
            return data
        data.sort(key=lambda tup: tup[0])
        return data


class SurfaceIntersector(object):
    """
    Base class for handling surface intersection methods and results.
    """

    def __init__(self, s1, s2, ftol=None):
        self._s1 = s1
        self._s2 = s2
        if ftol is None:
            ftol = Settings.ftol
        self._ftol = ftol
        self._icrvs = []

    @property
    def ncrvs(self):
        return len(self._icrvs)

    @property
    def success(self):
        if self.ncrvs > 0:
            return True
        return False

    @property
    def s1(self):
        return self._s1

    @property
    def s2(self):
        return self._s2

    @property
    def ftol(self):
        return self._ftol

    @property
    def icurves(self):
        return self._icrvs

    def get_icurve(self, indx=0):
        """
        Generate an :class:`.ICurve` for the specified intersection curve.

        :param int indx: Index of intersection curve.

        :return: Intersection curve.
        :rtype: :class:`.ICurve`
        """
        try:
            return self._icrvs[indx]
        except IndexError:
            return None

    def get_icurves(self):
        """
        Get all intersection curves.

        :return: List of intersection curves.
        :rtype: list
        """
        return self._icrvs

    def curve_nearest_point(self, pref, tol=None):
        """
        Get the index of the intersection curve that is nearest to the given
        reference point.

        :param array_like pref: Reference point.
        :param float tol: Distance tolerance.

        :return: Index of curve nearest point.
        :rtype: int
        """
        if isinstance(self, IntersectPlanePlane):
            return 0
        return si_curve_nearest_point(self, pref, tol)


class IntersectSurfacePlane(SurfaceIntersector):
    """
    Surface-plane intersection.

    :param surface: Surface to intersect.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param plane: Plane to intersect surface with.
    :type plane: :class:`.Plane`
    :param float ftol: Surface flatness tolerance.

    :var bool success: Status of intersection results.
    :var int ncrvs: Number of curves in the solution set.
    """

    def __init__(self, surface, plane, ftol=None):
        super(IntersectSurfacePlane, self).__init__(surface, plane, ftol)
        if CheckGeom.is_surface(surface) and CheckGeom.is_plane(plane):
            self._perform(surface, plane)

    def _perform(self, surface, plane):
        """
        Perform the surface-plane intersection.
        """
        results = intersect_surface_plane(surface, plane, self._ftol)
        self._icrvs = _build_si_curves(self, results)


class IntersectSurfaceSurface(SurfaceIntersector):
    """
    Surface-surface intersection.

    :param surface1: Surface 1 to intersect.
    :type surface1: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param surface2: Surface 2 to intersect.
    :type surface2: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float ftol: Surface flatness tolerance.

    :var bool success: Status of intersection results.
    :var int ncrvs: Number of curves in the solution set.
    """

    def __init__(self, surface1, surface2, ftol=None):
        super(IntersectSurfaceSurface, self).__init__(surface1, surface2, ftol)
        if CheckGeom.is_surface(surface1) and CheckGeom.is_surface(surface2):
            self._perform(surface1, surface2)

    def _perform(self, surface1, surface2):
        """
        Perform the surface-surface intersection.
        """
        results = intersect_surface_surface(surface1, surface2, self._ftol)
        self._icrvs = _build_si_curves(self, results)


class IntersectPlanePlane(SurfaceIntersector):
    """
    Plane-Plane intersection.
    """

    def __init__(self, plane1, plane2):
        super(IntersectPlanePlane, self).__init__(plane1, plane2, None)
        if CheckGeom.is_plane(plane1) and CheckGeom.is_plane(plane2):
            self._perform(plane1, plane2)

    def _perform(self, plane1, plane2):
        """
        Perform plane-plane intersection.
        """
        line = intersect_plane_plane(plane1, plane2)
        self._icrvs = [line]

    @property
    def line(self):
        try:
            return self._icrvs[0]
        except IndexError:
            return None


def _build_si_curves(si, results):
    """
    Build intersection curve results.
    """
    ncrvs = results[0]
    # crv_size = results[1]
    crv_ids = results[2]
    pnts = results[3]
    params1 = results[4]
    params2 = results[5]
    icrvs = []
    for i in range(0, ncrvs):
        # Collect points.
        p3d = [Point(pnts[j]) for j in crv_ids[i]]
        p2d_s1 = [Point2D(params1[j]) for j in crv_ids[i]]
        p2d_s2 = [Point2D(params2[j]) for j in crv_ids[i]]
        # Create new ICurve.
        icrv = _build_icurve(si, p3d, p2d_s1, p2d_s2)
        icrvs.append(icrv)
    # Check curves and return.
    return _check_curves(si, icrvs)


def _build_icurve(si, p3d, p2d_s1, p2d_s2):
    """
    Interpolate points and generate NURBS curves for intersection curve.
    """
    # Interpolate 3-D curve.
    cp, uk, p = global_curve_interpolation(p3d, 1, 'chord')
    c3d = NurbsCurve()
    c3d.set_cp(cp)
    c3d.set_knots(uk)
    c3d.set_deg(p)
    # Interpolate 2-D curves.
    cp, uk, p = global_curve_interpolation(p2d_s1, 1, 'chord')
    cp[:, -1] = 0.
    c2d_s1 = NurbsCurve()
    c2d_s1.set_cp(cp)
    c2d_s1.set_knots(uk)
    c2d_s1.set_deg(p)
    cp, uk, p = global_curve_interpolation(p2d_s2, 1, 'chord')
    cp[:, -1] = 0.
    c2d_s2 = NurbsCurve()
    c2d_s2.set_cp(cp)
    c2d_s2.set_knots(uk)
    c2d_s2.set_deg(p)
    # Set knot vectors of 2-D curves to match 3-D.
    c2d_s1.set_knots(c3d.uk)
    c2d_s2.set_knots(c3d.uk)
    # Create new ICurve.
    return ICurve(si.s1, si.s2, c3d, c2d_s1, c2d_s2, si.ftol)


def _check_curves(si, icrvs_in):
    """
    Check intersection curves.
    """
    # Check each intersection curve for overlapping segments. This can
    # sometimes result if intersection a closed surface right at the seam.
    icrvs = []
    for icrv in icrvs_in:
        # Do nothing if a line.
        if CheckGeom.is_line(icrv):
            icrvs.append(icrv)
            continue
        # First check if endpoints are equal.
        p0, p1 = icrv.eval(0.), icrv.eval(1.)
        if p0.is_equal(p1, icrv.ftol):
            icrvs.append(icrv)
            continue
        # TODO Python implementation
        # Project first point to the curve. If the distance is less than a
        # tolerance and the parameter is less than the end, create a new curve.
        # proj = ProjectPointToCurve(p0, icrv, False, True)
        # if not proj.success:
        #     icrvs.append(icrv)
        #     continue
        # for pi, ui, di in proj.results:
        #     if icrv.a < ui < icrv.b and di <= icrv.ftol:
        #         # Extract and build new curve.
        #         icrv = icrv.extract(icrv.a, ui, domain='global')
        #         p3d, p2d_s1, p2d_s2 = icrv.cp, icrv.cp2d_s1, icrv.cp2d_s2
        #         icrv = _build_icurve(si, p3d, p2d_s1, p2d_s2)
        #         break
        icrvs.append(icrv)

    # Sorted based on length.
    icrv_sort = []
    for icrv in icrvs:
        a = icrv.arc_length()
        icrv_sort.append([a, icrv])

    icrv_sort.sort(key=lambda lst: lst[0], reverse=True)

    return [icrv[1] for icrv in icrv_sort]
