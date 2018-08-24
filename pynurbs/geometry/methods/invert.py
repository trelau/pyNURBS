from __future__ import division

from numpy import dot, array, float64, zeros

from pynurbs.config import Settings
from pynurbs.geometry.methods.project import project_point_to_plane


def invert_point_on_plane(point, plane):
    """
    Invert a point on a plane.
    """
    _, _, proj = project_point_to_plane(point, plane)

    u, v = proj[0][1]
    return u, v


def invert_points_on_plane(pnts, plane):
    """
    Find the parameters of the points on a plane.

    :param pnts: Points to project.
    :type pnts: array_like
    :param plane: Plane to invert points on.
    :type plane: :class:`.Plane`

    :return: Parameters on plane (u, v) as NumPy array.
    :rtype: ndarray
    """
    ptol = Settings.ptol
    vu = plane.vu.ijk
    vv = plane.vv.ijk
    p0 = plane.p0.xyz
    vu2 = dot(vu, vu)
    vv2 = dot(vv, vv)

    pnts = array(pnts, dtype=float64)
    npts = pnts.shape[0]
    params = zeros((npts, 2), dtype=float64)
    for i in range(0, npts):
        pi = pnts[i]
        u = dot(pi - p0, vu) / vu2
        v = dot(pi - p0, vv) / vv2
        if abs(u) <= ptol:
            u = 0.
        if abs(v) <= ptol:
            v = 0.
        params[i, :] = [u, v]

    return params


def invert_point_on_surface(point, surface, ends=True):
    """
    Invert a point on a surface.
    """
    raise NotImplementedError('No Python implementation.')
    # _, npts, proj = project_point_to_surface(point, surface, ends)
    #
    # if npts == 0:
    #     return None, None
    #
    # u, v = proj[0][1]
    # return u, v
