from __future__ import print_function

try:
    import mayavi.mlab as mlab

    has_mayavi = True
except ImportError:
    has_mayavi = False
    mlab = None

from pynurbs.graphics.methods.plot_geometry import (plot_geom,
                                                    plot_triangular_mesh)
from pynurbs.geometry.checker import CheckGeom


class Viewer(object):
    """
    View objects using Mayavi visualization tools.
    """
    _geom = set()
    _triplots = []
    _ctol = 0.001
    _stol = 0.005
    _scale = 0.15

    @property
    def ctol(self):
        return self._ctol

    @ctol.setter
    def ctol(self, ctol):
        self._ctol = float(ctol)

    @property
    def stol(self):
        return self._stol

    @stol.setter
    def stol(self, stol):
        self._stol = float(stol)

    @classmethod
    def clear(cls):
        """
        Clear entities from viewer.
        """
        cls._geom.clear()
        cls._triplots = []

    @classmethod
    def set_scale(cls, scale=0.15):
        """
        Set scale for plotting points.

        :param float scale: Scale value.
        """
        cls._scale = float(scale)

    @classmethod
    def show(cls, clear=True):
        """
        Show the scene.
        """
        if not has_mayavi:
            print('Failed to show items. Mayavi is not available')
            return None
        fig = mlab.figure(size=(800, 600), bgcolor=(1., 1., 1.))
        for g in cls._geom:
            plot_geom(g, fig, cls._ctol, cls._stol, cls._scale)
        for verts, tri, rep, color in cls._triplots:
            plot_triangular_mesh(verts, tri, fig, rep, color)
        mlab.show()
        if clear:
            cls.clear()

    @classmethod
    def add_items(cls, *items):
        """
        Add item to the viewer.

        :param items: Item(s) to add.
        """
        for item in items:
            if CheckGeom.is_geom(item):
                cls.add_geom(item)

    @classmethod
    def add_geom(cls, *args):
        """
        Add geometry to the viewer.

        :param args: Geometry entity(s) to add.
        """
        for g in args:
            cls._geom.add(g)

    @classmethod
    def add_triplot(cls, verts, tri, rep='wireframe', color=None):
        cls._triplots.append((verts, tri, rep, color))
