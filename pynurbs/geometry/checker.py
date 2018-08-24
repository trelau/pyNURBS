from numpy import ndarray

from pynurbs.geometry.bezier_curve import BezierCurve
from pynurbs.geometry.bezier_surface import BezierSurface
from pynurbs.geometry.geom import Geometry
from pynurbs.geometry.icurve import ICurve
from pynurbs.geometry.line import Line
from pynurbs.geometry.nurbs_curve import NurbsCurve
from pynurbs.geometry.nurbs_surface import NurbsSurface
from pynurbs.geometry.plane import Plane
from pynurbs.geometry.point import Point, Point2D
from pynurbs.geometry.system import System
from pynurbs.geometry.vector import Vector


class CheckGeom(object):
    """
    Geometry checker.
    """

    @staticmethod
    def is_geom(geom):
        """
        Check if the entity is geometry.

        :param geom: Entity.

        :return: *True* if entity is geometry, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Geometry)

    @staticmethod
    def is_point_like(geom):
        """
        Check to see if the entity is point_like, which includes a Point
        instance, a tuple, list, or NumPy array of shape 1 x 3 that contains
        xyz coordinates.

        :param geom: Entity
        :type geom: :class:`.Point`  or array_like

        :return: *True* if the entity is point_like, *False* if not.
        :rtype: bool
        """
        if isinstance(geom, Point):
            return True
        if isinstance(geom, (tuple, list, ndarray)):
            return len(geom) == 3
        return False

    @staticmethod
    def is_point(geom):
        """
        Check to see if the entity is a point.

        :param geom: Geometric entity.
        :return: *True* if the entity is a point, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point)

    @staticmethod
    def is_point2d(geom):
        """
        Check to see if the entity is a Point2D.

        :param geom: Geometric entity.
        :return: *True* if the entity is a Point2D, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Point2D)

    @staticmethod
    def to_point(geom):
        """
        Check to see if the entity is a :class:`.Point` instance and return
        the instance if it is. If the entity is point_like, create a new
        :class:`.Point` instance and return it.

        :param geom: Geometric entity or possible array.

        :return: The Point instance if already a point, or a new Point in
            case it is array_like.
        :rtype: :class:`.Point`
        """
        if isinstance(geom, Point):
            return geom
        elif isinstance(geom, (tuple, list, ndarray)):
            return Point(geom)
        return None

    @staticmethod
    def to_points(geoms):
        """
        Check to see if the entities are a :class:`.Point` instance or
        convert them if they are point_like.

        :param list geoms: List of point_like entities.

        :return: List of :class:`.Point` instances.
        :rtype: list
        """
        return [CheckGeom.to_point(p) for p in geoms]

    @staticmethod
    def is_vector(geom):
        """
        Check to see if the entity is a vector.

        :param geom: Geometric entity.
        :return: *True* if the entity is a vector, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Vector)

    @staticmethod
    def is_plane(geom):
        """
        Check to see if the entity is a plane.

        :param geom: Geometric entity.
        :return: *True* if the entity is a plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Plane)

    @staticmethod
    def is_line(geom):
        """
        Check to see if the entity is a line.

        :param geom: Geometric entity.
        :return: *True* if the entity is a line, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, Line)

    @staticmethod
    def is_curve(geom):
        """
        Check to see if the entity is a curve.

        :param geom: Geometric entity.
        :return: *True* if the entity is a curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, (BezierCurve, NurbsCurve, ICurve))

    @staticmethod
    def is_curve_like(geom):
        """
        Check to see if the entity is a curve or line.

        :param geom: Geometric entity.
        :return: *True* if the entity is a line or curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, (BezierCurve, NurbsCurve, Line, ICurve))

    @staticmethod
    def is_icurve(geom):
        """
        Check to see if the entity is an intersection curve.

        :param geom: Geometric entity.
        :return: *True* if the entity is an intersection curve, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, ICurve)

    @staticmethod
    def is_surface(geom):
        """
        Check to see if the entity is a surface.

        :param geom: Geometric entity.
        :return: *True* if the entity is a surface, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, (BezierSurface, NurbsSurface))

    @staticmethod
    def is_surface_like(geom):
        """
        Check to see if the entity is a surface or a plane.

        :param geom: Geometric entity.
        :return: *True* if the entity is a surface or plane, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, (BezierSurface, NurbsSurface, Plane))

    @staticmethod
    def is_system(geom):
        """
        Check to see if the entity is a system.

        :param geom:

        :return: *True* if system, *False* if not.
        :rtype: bool
        """
        return isinstance(geom, System)
