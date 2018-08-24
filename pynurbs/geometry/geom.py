from copy import deepcopy

from numpy.random import rand


class Geometry(object):
    """
    Base class for geometry.

    :param str etype: Type of geometry.

    :var str etype: Type of geometry.
    :var tuple color: RGB color value.
    """

    def __init__(self, etype):
        self._etype = etype
        self._color = tuple(rand(1, 3)[0])
        self._representation = 'surface'
        self._opacity = 1.

    @property
    def etype(self):
        return self._etype

    @property
    def color(self):
        return self._color

    def set_color(self, r, g, b):
        """
        Set color (0. <= r, g, b <= 1.).

        :param float r: Red.
        :param float g: Green.
        :param float b: Blue.
        """
        self._color = (float(r), float(g), float(b))

    def copy(self):
        """
        Return a deepcopy of the geometry instance.

        :return: Copy of geometry instance.
        """
        return deepcopy(self)

    @property
    def representation(self):
        return self._representation

    @property
    def opacity(self):
        return self._opacity

    def set_representation(self, representation='surface'):
        """
        Set the tessellation representation for visualization.

        :param str representation: Representation ('surface', 'wireframe',
            'points', 'mesh', 'fancymesh').
        """
        self._representation = representation

    def set_opacity(self, opacity):
        """
        Set the opacity for graphics.

        :param float opacity: Level of opacity.
        """
        self._opacity = opacity
