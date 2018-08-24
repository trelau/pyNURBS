from pynurbs.geometry.methods.map_surface import surface_distance_map


class SurfaceMap(object):
    """
    Surface mapping class.
    """

    def __init__(self, surface, build=True, n=100, h=None):
        self._s = surface
        self._interp_udist = None
        self._interp_vdist = None
        self._interp_uparam = None
        self._interp_vparam = None
        if build:
            self.build_distance_map(n, h)

    def build_distance_map(self, n=100, h=None):
        """
        Build a distance map of the surface.

        :param int n:
        :param float h:

        :return:
        """
        results = surface_distance_map(self._s, n, h)
        self._interp_udist = results[0]
        self._interp_vdist = results[1]
        self._interp_uparam = results[2]
        self._interp_vparam = results[3]
        return True

    def eval_udist(self, u, v):
        """
        Evaluate u-distance given parameters.

        :param u:
        :param v:

        :return:
        """
        try:
            return float(self._interp_udist(u, v)[0, 0])
        except (AttributeError, IndexError):
            return None

    def eval_vdist(self, u, v):
        """
        Evaluate v-distance given parameters.

        :param u:
        :param v:

        :return:
        """
        try:
            return float(self._interp_vdist(u, v)[0, 0])
        except (AttributeError, IndexError):
            return None

    def eval_uparam(self, s, t):
        """
        Evaluate u-parameter given distances.

        :param s:
        :param t:

        :return:
        """
        return float(self._interp_uparam(s, t))

    def eval_vparam(self, s, t):
        """
        Evaluate t-parameter given distance.

        :param s:
        :param t:

        :return:
        """
        return float(self._interp_vparam(s, t))
