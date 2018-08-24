from numpy import array, float64

from pynurbs.geometry.checker import CheckGeom
from pynurbs.geometry.methods.project import (project_point_to_line,
                                              project_point_to_plane)
from pynurbs.geometry.point import Point


class ProjectGeom(object):
    """
    Project geometry.
    """

    @staticmethod
    def _update(point, proj, update=False):
        """
        Return projection or update point location.
        """
        if proj.success and update:
            point.set_xyz(proj.nearest_point.xyz)
        return proj

    @staticmethod
    def point(pnt, geom, ends=True, all_pnts=False, update=False):
        """
        Project a point to a curve or surface.

        :param pnt: Point to project.
        :type pnt: :class:`.Point`
        :param geom: Curve or surface entity to project point to.
        :param bool ends: Option to project point to the nearest end point of a
            curve, or to the nearest boundary curve of a surface if no
            orthogonal projection is found.
        :param bool all_pnts: Option to return only the nearest point during
            the subdivision process (*False*), or all points (*True*).
        :param bool update: Option to update the location of the *point*
            rather than returning a projection object. If *True*, then the
            nearest point will be used to update the point location and a
            boolean will be returned indicating a successful *True* or
            unsuccessful *False* projection.

        :return: Projection object depending on *geom* type.
        """
        pnt = CheckGeom.to_point(pnt)

        # Project point to curve.
        if CheckGeom.is_curve_like(geom):
            proj = ProjectPointToCurve(pnt, geom, ends, all_pnts)
            return ProjectGeom._update(pnt, proj, update)

        # Project point to surface.
        if CheckGeom.is_surface(geom):
            proj = ProjectPointToSurface(pnt, geom, ends, all_pnts)
            return ProjectGeom._update(pnt, proj, update)

        # Project point to plane.
        if CheckGeom.is_plane(geom):
            proj = ProjectPointToPlane(pnt, geom)
            return ProjectGeom._update(pnt, proj, update)
        return False

    @staticmethod
    def invert(point, geom, ends=True):
        """
        Return the parameters of the nearest orthogonal projection to the
        geometry.

        :param point: Point to project.
        :type point: :class:`.Point`
        :param geom: Curve or surface entity to find parameters.
        :param bool ends: Option to project point to the nearest end point of a
            curve, or to the nearest boundary curve of a surface if no
            orthogonal projection is found.

        :return: Parameter(s) on the geometry of the nearest orthogonal
            projection. For a curve a single float *u* is returned,
            for a surface a tuple (u, v) is returned.
        :rtype: float or tuple

        .. note::
            *None* is returned if no points are found or the method fails. A
            tuple (*None*, *None*) is returned for a surface.
        """
        point = CheckGeom.to_point(point)
        if not CheckGeom.is_point(point):
            if CheckGeom.is_curve_like(geom):
                return None
            return None, None

        if CheckGeom.is_curve_like(geom):
            proj = ProjectPointToCurve(point, geom, ends)
            if not proj.success:
                return None
            return proj.nearest_param
        if CheckGeom.is_surface(geom):
            proj = ProjectPointToSurface(point, geom, ends)
            if not proj.success:
                return None, None
            return proj.nearest_param
        if CheckGeom.is_plane(geom):
            proj = ProjectPointToPlane(point, geom)
            if not proj.success:
                return None, None
            return proj.nearest_param
        return None

    @staticmethod
    def invert_array(pnts, geom, ends=True):
        """
        Invert the array of points to the given geometry.

        :param array_like pnts: Points to invert.
        :param geom: Geometry to invert points on.
        :param bool ends: Option to project point to the nearest end point of a
            curve, or to the nearest boundary curve of a surface if no
            orthogonal projection is found.

        :return: Array of parameters of inverted points.
        :rtype: ndarray
        """
        uv_data = []
        for pi in pnts:
            uv = ProjectGeom.invert(pi, geom, ends)
            uv_data.append(uv)
        return array(uv_data, dtype=float64)

    @staticmethod
    def point_on_geom(point, geom, tol=None):
        """
        Check to see if the point is on the geometry.

        :param point:
        :param geom:
        :param tol:

        :return: *True* if ont geometry, *False* if not and parameter(s).
        :rtype: tuple
        """
        point = CheckGeom.to_point(point)
        if CheckGeom.is_curve_like(geom):
            u = ProjectGeom.invert(point, geom, False)
            if u is None:
                return False, None
            p = geom.eval(u, domain='global')
            is_on = p.is_equal(point, tol)
            return is_on, u
        elif CheckGeom.is_surface_like(geom):
            u, v = ProjectGeom.invert(point, geom, False)
            if None in [u, v]:
                return False, None, None
            p = geom.eval(u, v, domain='global')
            is_on = p.is_equal(point, tol)
            return is_on, u, v
        return False, None

    @staticmethod
    def points_to_curve(curve, pnts, copy=True, ends=True):
        """
        Project points to a curve.

        :param curve_like curve: Curve to project points to.
        :param list pnts: List of point_like entities. Will be converted to
            :clas:`.Point` instances if not already.
        :param bool copy: Option to copy the points and project them rather
            than the original points.
        :param bool ends: Option to project point to the nearest end point of a
            curve.

        :return: List of the projected points.
        :rtype: list

        ..note:
        If a projection fails the point location is not updated.
        """
        if not CheckGeom.is_curve_like(curve):
            return []

        # Make sure a list of Point instances is provided.
        pnts = CheckGeom.to_points(pnts)

        # Copy points if desired.
        if copy:
            pnts = [p.copy() for p in pnts]

        # Project each point to the curve and update its location.
        for p in pnts:
            ProjectGeom.point(p, curve, ends, update=True)

        return pnts

    @staticmethod
    def points_to_surface(surface, pnts, copy=True, ends=True):
        """
        Project points to a surface.

        :param surface_like surface: Surface to project points to.
        :param list pnts: List of point_like entities. Will be converted to
            :clas:`.Point` instances if not already.
        :param bool copy: Option to copy the points and project them rather
            than the original points.
        :param bool ends: Option to project point to the nearest boundary
            curve of a surface.

        :return: List of the projected points.
        :rtype: list

        ..note:
        If a projection fails the point location is not updated.
        """
        if not CheckGeom.is_surface_like(surface):
            return []

        # Make sure a list of Point instances is provided.
        pnts = CheckGeom.to_points(pnts)

        # Copy points if desired.
        if copy:
            pnts = [p.copy() for p in pnts]

        # Project each point to the surface and update its location.
        for p in pnts:
            ProjectGeom.point(p, surface, ends, update=True)

        return pnts


class PointProjector(object):
    """
    Base for handling point projections to curves and surfaces.
    """

    def __init__(self):
        self._npts = 0
        self._results = []

    def _set_results(self, nsub, npts, results):
        """
        Set projection results.
        """
        if npts > 0:
            self._subdivisions = nsub
            self._npts = npts
            # Replace the ndarray with a point instance.
            for i in range(npts):
                results[i][0] = Point(results[i][0])
            self._results = results
            # Sort results by distance.
            self._results.sort(key=lambda lst: lst[2])

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
            return None
        return [results[0] for results in self._results]

    @property
    def parameters(self):
        if self._npts <= 0:
            return None
        return [results[1] for results in self._results]

    @property
    def nearest_point(self):
        return self._results[0][0]

    @property
    def nearest_param(self):
        return self._results[0][1]

    @property
    def dmin(self):
        return self._results[0][2]

    @property
    def results(self):
        return self._results

    def point(self, indx=0):
        """
        Return the point result by index.

        :param int indx: Index for point selection.

        :return: Projection point at index.
        :rtype: :class:`.Point`
        """
        if indx > self._npts - 1:
            return self._results[-1][0]
        return self._results[indx][0]

    def parameter(self, indx=0):
        """
        Return the parameter result by index.

        :param int indx: Index for parameter selection.

        :return: Projection parameter at index. For a curve projection a
            single float *u* will be returned. For a surface projection result
            a tuple containing the *u* and *v* parameters will be returned
            (u,v).
        :rtype: float or tuple
        """
        if indx > self._npts - 1:
            return self._results[-1][1]
        return self._results[indx][1]

    def distance(self, indx=0):
        """
        Return the projection distance by index.

        :param int indx: Index for distance selection.

        :return: Projection distance between original point and projection
            result.
        :rtype: float
        """
        if indx > self._npts - 1:
            return self._results[-1][2]
        return self._results[indx][2]


class ProjectPointToPlane(PointProjector):
    """
    Project a point to a plane.

    :param point: Point to project.
    :type point: :class:`.Point`
    :param plane: Plane to project point to.
    :type plane: :class:`.Plane`

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all points in the solution set as :class:`.Point`
        instances.
    :var list parameters: List of all parameters in the solution set.
    :var nearest_point: Nearest point in the solution set.
    :type nearest_point: :class:`.Point`
    :var float nearest_param: Nearest parameter in the solution set.
    :var float min_distance: Distance to nearest point in solution set.
    """

    def __init__(self, point, plane):
        super(ProjectPointToPlane, self).__init__()
        point = CheckGeom.to_point(point)
        if CheckGeom.is_point(point) and CheckGeom.is_plane(plane):
            self._perform(point, plane)

    def _perform(self, point, plane):
        """
        Perform projection.
        """
        nsub, npts, results = project_point_to_plane(point, plane)
        self._set_results(nsub, npts, results)


class ProjectPointToCurve(PointProjector):
    """
    Project a point to a curve.

    :param point: Point to project.
    :type point: :class:`.Point`
    :param curve: Curve to project point to.
    :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param bool ends: Option to project point to nearest end of curve if
        no orthogonal projection is found.
    :param bool all_pnts: Option to return only the nearest point during the
        subdivision process (*False*), or all points (*True*).

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all points in the solution set as :class:`.Point`
        instances.
    :var list parameters: List of all parameters in the solution set.
    :var nearest_point: Nearest point in the solution set.
    :type nearest_point: :class:`.Point`
    :var float nearest_param: Nearest parameter in the solution set.
    :var float min_distance: Distance to nearest point in solution set.
    """

    def __init__(self, point, curve, ends=True, all_pnts=False):
        super(ProjectPointToCurve, self).__init__()
        point = CheckGeom.to_point(point)
        if CheckGeom.is_point(point) and CheckGeom.is_curve_like(curve):
            self._perform(point, curve, ends, all_pnts)

    def _perform(self, point, curve, ends, all_pnts):
        """
        Perform the point to curve projection.
        """
        if CheckGeom.is_icurve(curve):
            raise NotImplementedError('No Python implementation.')
            # nsub, npts, results = project_point_to_icurve(point, curve, ends,
            #                                               all_pnts)
        elif CheckGeom.is_line(curve):
            nsub, npts, results = project_point_to_line(point, curve)
        else:
            raise NotImplementedError('No Python implementation.')
            # nsub, npts, results = project_point_to_curve(point, curve, ends,
            #                                              all_pnts)
        self._set_results(nsub, npts, results)


class ProjectPointToSurface(PointProjector):
    """
    Project a point to a surface.

    :param point: Point to project.
    :type point: :class:`.Point`
    :param surface: Surface to project point to.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param bool ends: Option to project point to surface boundary if no
        orthogonal point is found.
    :param bool all_pnts: Option to return only the nearest point during the
        subdivision process (*False*), or all points (*True*).

    :var bool success: Status of intersection results. This only indicates that
        there is at least one solution found.
    :var int npts: Number of points in the solution set.
    :var list points: List of all points in the solution set as :class:`.Point`
        instances.
    :var list parameters: List of all parameters in the solution set.
    :var nearest_point: Nearest point in the solution set.
    :type nearest_point: :class:`.Point`
    :var tuple nearest_param: Nearest parameter in the solution set. A surface
        projection result will return a tuple containing the *u* and *v*
        parameters (u, v).
    :var float min_distance: Distance to nearest point in solution set.
    """

    def __init__(self, point, surface, ends=True, all_pnts=False):
        super(ProjectPointToSurface, self).__init__()
        point = CheckGeom.to_point(point)
        if CheckGeom.is_point(point) and CheckGeom.is_surface(surface):
            self._perform(point, surface, ends, all_pnts)

    def _perform(self, point, surface, ends, all_pnts):
        """
        Perform the point to surface projection.
        """
        raise NotImplementedError('No Python implementation.')
        # nsub, npts, results = project_point_to_surface(point, surface, ends,
        #                                                all_pnts)
        # self._set_results(nsub, npts, results)
