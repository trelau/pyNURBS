from numpy import add, array, float64, subtract
from numpy.linalg import norm

from pynurbs.config import Settings
from pynurbs.geometry.geom import Geometry


class Point(Geometry):
    """
    Point.

    3-D point.

    :param array_like xyz: Coordinates of point.

    :var ndarray xyz: Coordinates of point (x, y, z).
    :var float x: x-coordinate of point.
    :var float y: y-coordinate of point.
    :var float z: z-coordinate of point.
    """

    def __init__(self, xyz):
        super(Point, self).__init__('point')
        self._xyz = array(xyz, dtype=float64)
        self._x, self._y, self._z = self._xyz
        self.set_color(1, 1, 0)

    def __str__(self):
        return 'Point = ({0}, {1}, {2})'.format(*self._xyz)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self._xyz, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self._xyz:
            yield elm

    def __len__(self):
        return len(self._xyz)

    def __getitem__(self, item):
        return self._xyz[item]

    def __eq__(self, other):
        return self.is_equal(other, 0.)

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def xyz(self):
        return self._xyz

    @xyz.setter
    def xyz(self, xyz):
        self.set_xyz(xyz)

    @property
    def x(self):
        return self._x

    @x.setter
    def x(self, x):
        self.set_x(x)

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y):
        self.set_y(y)

    @property
    def z(self):
        return self._z

    @z.setter
    def z(self, z):
        self.set_z(z)

    def set_x(self, x):
        """
        Set point x-coordinate.

        :param float x: x-coordinate of point.
        """
        self._x = float64(x)
        self._xyz[0] = self._x

    def set_y(self, y):
        """
        Set point y-coordinate.

        :param float y: y-coordinate of point.
        """
        self._y = float64(y)
        self._xyz[1] = self._y

    def set_z(self, z):
        """
        Set point z-coordinate.

        :param float z: z-coordinate of point.
        """
        self._z = float64(z)
        self._xyz[2] = self._z

    def set_xyz(self, xyz):
        """
        Set point coordinates.

        :param array_like xyz: Point coordinates.
        """
        self._xyz = array(xyz, dtype=float64)
        self._x, self._y, self._z = self._xyz

    def dist2pnt(self, p):
        """
        Calculate the distance between points.

        :param p: Other point.
        :type p: :class:`.Point` or array_like

        :return: Distance to point.
        :rtype: float
        """
        xyz = array(p, dtype=float64)
        return norm(xyz - self._xyz)

    def is_equal(self, p, tol=None):
        """
        Check to see if the points are coincident.

        :param p: Other point.
        :type p: :class:`.Point` or array_like
        :param float tol: Tolerance for checking point coincidence. If *None*
            is provided then the default setting *gtol* will be used.

        :return: *True* if points are coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.gtol
        if self.dist2pnt(p) <= tol:
            return True
        return False


class Point2D(Geometry):
    """
    2-D point.

    :param array_like uv: Coordinates of point.

    :var ndarray uv: Coordinates of point (u, v).
    :var float u: u-coordinate of point.
    :var float v: v-coordinate of point.
    """

    def __init__(self, uv):
        super(Point2D, self).__init__('point')
        self._uv = array(uv, dtype=float64)
        self._u, self._v, = self._uv
        # Use this as a hack in 3-D curves.
        self._uvw = array([self._uv[0], self._uv[1], 0.], dtype=float64)
        self.set_color(1, 1, 0)

    def __str__(self):
        return 'Point2D = ({0}, {1})'.format(*self._uv)

    def __array__(self, dtype=float64, copy=True, order=None, subok=False,
                  ndmin=0):
        return array(self._uvw, dtype=dtype, copy=copy, order=order,
                     subok=subok, ndmin=ndmin)

    def __iter__(self):
        for elm in self._uvw:
            yield elm

    def __len__(self):
        return len(self._uvw)

    def __getitem__(self, item):
        return self._uvw[item]

    def __eq__(self, other):
        return self.is_equal(other, 0.)

    def __add__(self, other):
        return add(self, other)

    def __sub__(self, other):
        return subtract(self, other)

    @property
    def uv(self):
        return self._uv

    @property
    def u(self):
        return self._u

    @property
    def v(self):
        return self._v

    def set_u(self, u):
        """
        Set point u-coordinate.

        :param float u: u-coordinate of point.
        """
        self._u = float64(u)
        self._uv[0] = self._u
        self._uvw[0] = self._u

    def set_v(self, v):
        """
        Set point v-coordinate.

        :param float v: v-coordinate of point.
        """
        self._v = float64(v)
        self._uv[1] = self._v
        self._uvw[1] = self._v

    def set_uv(self, uv):
        """
        Set point coordinates.

        :param array_like uv: Point coordinates.
        """
        self._uv = array(uv, dtype=float64)
        self._u, self._v = self._uv
        self._uvw[0] = self._u
        self._uvw[1] = self._v

    def dist2pnt(self, p):
        """
        Calculate the distance between points.

        :param p: Other point.
        :type p: :class:`.Point2D` or array_like

        :return: Distance to point.
        :rtype: float
        """
        uv = array(p, dtype=float64)
        if uv.shape[0] > 2:
            uv = uv[:2]
        return norm(uv - self._uv)

    def is_equal(self, p, tol=None):
        """
        Check to see if the points are coincident.

        :param p: Other point.
        :type p: :class:`.Point2D` or array_like
        :param float tol: Tolerance for checking point coincidence. If *None*
            is provided then the default setting *ptol* will be used.

        :return: *True* if points are coincident, *False* if not.
        :rtype: bool
        """
        if tol is None:
            tol = Settings.ptol
        if self.dist2pnt(p) <= tol:
            return True
        return False
