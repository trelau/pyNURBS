from pynurbs.geometry.methods.misc import is_array_type
from pynurbs.geometry.point import Point


class System(object):
    """
    Coordinate system (prelim only).
    """

    def __init__(self, origin, vx, vy, vz):
        self._origin = origin
        self._vx = vx
        self._vy = vy
        self._vz = vz

    @property
    def vx(self):
        return self._vx

    @property
    def vy(self):
        return self._vy

    @property
    def vz(self):
        return self._vz

    def eval(self, dx, dy, dz, rtype='Point', *args, **kwargs):
        """
        Evaluate a point relative to the system.

        :param dx:
        :param dy:
        :param dz:
        :param rtype:
        :param args:
        :param kwargs:

        :return:
        """
        pnt = (self._origin.xyz + dx * self._vx.ijk +
               dy * self._vy.ijk + dz * self._vz.ijk)
        if is_array_type(rtype):
            return pnt
        else:
            return Point(pnt)
