from numpy import zeros, int32, float64, array

from pynurbs.geometry.checker import CheckGeom
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.methods.tessellate import adaptive_curve_tessellate
from pynurbs.geometry.methods.tessellate import adaptive_surface_tessellate


class TessGeom(object):
    """
    Tessellate geometry.
    """

    @staticmethod
    def curve_adaptive(curve, tol=None):
        """
        Tessellate a curve using adaptive subdivision.

        :param curve: Curve to tessellate.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        :param float tol: Tolerance to use for flatness criteria.

        :return: Adaptive curve tessellation.
        :rtype: :class:`.CurveTessAdaptive`
        """
        if CheckGeom.is_curve_like(curve):
            return CurveTessAdaptive(curve, tol)
        return None

    @staticmethod
    def surface_adaptive(surface, tol=0.01):
        """
        Tessellate a surface using adaptive subdivision.

        :param surface: Surface to tessellate.
        :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
        :param float tol: Tolerance to use for flatness criteria.

        :return: Adaptive surface tessellation.
        :rtype: :class:`.SurfaceTessAdaptive`
        """
        if CheckGeom.is_surface(surface):
            return SurfaceTessAdaptive(surface, tol)
        return None


class CurveTessAdaptive(Geometry):
    """
    Adaptive curve tessellation.

    :param curve: Curve to tessellate.
    :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param float tol: Tolerance for checking curve flatness.

    :var bool success: Status of tessellation results.
    :var int nverts: Number of vertices.
    :var ndarray verts: Array of vertices.
    """

    def __init__(self, curve, tol=None):
        super(CurveTessAdaptive, self).__init__('curve_tessellation')
        self._verts = zeros((0, 3), dtype=float64)
        self._params = zeros(0, dtype=float64)
        if CheckGeom.is_curve_like(curve):
            self._perform(curve, tol)

    def _perform(self, curve, tol):
        """
        Perform curve tessellation.
        """
        if CheckGeom.is_icurve(curve):
            curve = curve.crv3d
        if CheckGeom.is_line(curve):
            verts = array([curve.eval(0., rtype='ndarray'),
                           curve.eval(curve.length, rtype='ndarray')],
                          dtype=float64)
            params = array([0., curve.length], dtype=float64)
        else:
            verts, params = adaptive_curve_tessellate(curve, tol)
        self._verts = verts
        self._params = params

    @property
    def success(self):
        if self.nverts > 0:
            return True
        return False

    @property
    def verts(self):
        return self._verts

    @property
    def params(self):
        return self._params

    @property
    def nverts(self):
        return self._verts.shape[0]


class SurfaceTessAdaptive(Geometry):
    """
    Adaptive surface tessellation.

    :param surface: Surface to tessellate.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float tol: Tolerance for checking surface flatness.

    :var bool success: Status of tessellation results.
    :var int nverts: Number of triangle vertices.
    :var ndarray verts: Array of triangles vertices.
    :var int ntri: Number of triangles.
    :var ndarray triangles: Connectivity of each triangle. Each row is a
        triangle referencing three indices in the *verts* array.
    :var str representation: Representation for graphics.
    """

    def __init__(self, surface, tol=0.01):
        super(SurfaceTessAdaptive, self).__init__('surface_tessellation')
        self._verts = zeros((0, 3), dtype=float64)
        self._triangles = zeros((0, 3), dtype=int32)
        if CheckGeom.is_surface(surface):
            self._perform(surface, tol)

    @property
    def success(self):
        if self.ntri > 0:
            return True
        return False

    @property
    def verts(self):
        return self._verts

    @property
    def nverts(self):
        return self._verts.shape[0]

    @property
    def triangles(self):
        return self._triangles

    @property
    def ntri(self):
        return self._triangles.shape[0]

    def _perform(self, surface, tol=0.01):
        """
        Perform surface tessellation.
        """
        # Try the method three times in case arrays need to be resized.
        verts, triangles = adaptive_surface_tessellate(surface, tol)
        self._verts = verts
        self._triangles = triangles
