from numpy import max as np_max
from numpy import min as np_min
from numpy import zeros, float64
from numpy.linalg import norm

from pynurbs.geometry.methods.intersect_bbox import bboxes_intersect


class BoundingBox(object):
    """
    3-dimensional axis-aligned bounding box.

    The BoundingBox class is typically used to test for possible intersection
    of curves and surfaces within subdivision methods.

    :var float xmin: Minimum value along x-axis.
    :var float xmax: Maximum value along x-axis.
    :var float ymin: Minimum value along y-axis.
    :var float ymax: Maximum value along y-axis.
    :var float zmin: Minimum value along z-axis.
    :var float zmax: Maximum value along z-axis.
    :var ndarray bounds: A 3 x 2 NumPy array containing the minimum and maximum
        values along each each. Each row is an axis and the columns are minimum
        and maximum values, respectively.
    """

    def __init__(self):
        self._bounds = zeros((3, 2), dtype=float64)

    @property
    def xmin(self):
        return self._bounds[0, 0]

    @property
    def xmax(self):
        return self._bounds[0, 1]

    @property
    def ymin(self):
        return self._bounds[1, 0]

    @property
    def ymax(self):
        return self._bounds[1, 1]

    @property
    def zmin(self):
        return self._bounds[2, 0]

    @property
    def zmax(self):
        return self._bounds[2, 1]

    @property
    def bounds(self):
        return self._bounds

    def clear(self):
        """
        Clear bounding box.
        """
        self._bounds = zeros((3, 2), dtype=float64)

    def set_bounds(self, bmin, bmax):
        """
        Set the bounds.

        :param ndarray bmin: Minimum bounds.
        :param ndarray bmax: Maximum bounds.
        """
        self._bounds[:, 0] = bmin
        self._bounds[:, 1] = bmax

    def add_curve(self, curve):
        """
        Add a curve to the bounding box.

        :param curve: Curve to add to bounding box.
        :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
        """
        cp = curve.cp
        bmin = np_min(cp, axis=0)
        bmax = np_max(cp, axis=0)
        self.set_bounds(bmin, bmax)

    def add_surface(self, surface):
        """
        Add a surface to the bounding box.

        :param surface: Surface to add to bounding box.
        :type surface: :class:`.BezierCurve` or :class:`.NurbsSurface`
        """
        d1, d2 = surface.n + 1, surface.m + 1
        cp = surface.cp
        cp = cp.reshape(d1 * d2, -1)
        bmin = np_min(cp, axis=0)
        bmax = np_max(cp, axis=0)
        self.set_bounds(bmin, bmax)

    def diagonal_length(self):
        """
        Calculate the diagonal length of the bounding box.


        :return: Length of box diagonal.
        :type: float
        """
        return norm(self._bounds[:, 1] - self._bounds[:, 0])

    def intersects(self, bbox):
        """
        Test if the bounding box intersects with another.

        :param bbox: Another bounding box.
        :type bbox: :class:`.BoundingBox`

        :return: *True* if they intersect, *False* if not.
        :rtype: bool
        """
        return bboxes_intersect(self, bbox)
