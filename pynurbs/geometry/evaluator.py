from numpy import array, float64

from pynurbs.config import Settings
from pynurbs.geometry.checker import CheckGeom
from pynurbs.geometry.creator import CreateGeom
from pynurbs.geometry.methods.calculate import polygon_area2D, triangle_area
from pynurbs.geometry.methods.misc import is_array_type
from pynurbs.geometry.projector import ProjectGeom


class EvalGeom(object):
    """
    Geometry evaluator.
    """

    @staticmethod
    def curve_eval2d(curve, u, sref=None, rtype='Point', domain='local',
                     tol=None):
        """
        Evaluate a curve in the parametric space of the surface.

        :param curve_like curve: Curve to evaluate.
        :param float u: Curve parameter.
        :param surface_like sref: The surface to return the 2-D parameters on.
        :param str rtype: Option to return a NumPy array or Point2D instance
            (rtype = 'Point' or 'ndarray').
        :param str domain: Option to use local (0 <= u <= 1) or global
            (a <= u <= b) domain ('local', 'l', 'global', 'g').
        :param float tol: Tolerance for point refinement (ICurve only).

        :return: Parameters on surface at curve parameter.
        :rtype: tuple
        """
        if not CheckGeom.is_curve_like(curve) or not \
                CheckGeom.is_surface_like(sref):
            return None

        if tol is None:
            tol = Settings.gtol / 100.

        # ICurve already has built-in method.
        if CheckGeom.is_icurve(curve) and curve.has_surf(sref):
            return curve.eval2d(u, rtype, domain, tol, sref)

        # Evaluate curve and invert on surface.
        p3d = curve.eval(u, domain=domain)
        u, v = ProjectGeom.invert(p3d, sref, True)
        if None in [u, v]:
            return None

        # Return desired type.
        if is_array_type(rtype):
            return array([u, v], dtype=float64)
        return CreateGeom.point2d((u, v))

    @staticmethod
    def polygon_area2d(xy):
        """
        Calculate the signed area of a 2D (planar) polygon using the shoelace
        formula.

        :param array_like xy: Ordered array of points defining polygon
            boundary.

        :return: Area of polygon.
        :rtype: float
        """
        return polygon_area2D(xy)

    @staticmethod
    def triangle_area(t):
        """
        Calculate the area of a triangle given by 3 points.

        :param array_like t: Array containing triangle vertices.
        """
        return triangle_area(t)
