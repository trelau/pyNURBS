import warnings

import mayavi.mlab as mlab
import numpy as np

from pynurbs.geometry.bezier_curve import BezierCurve
from pynurbs.geometry.bezier_surface import BezierSurface
from pynurbs.geometry.icurve import ICurve
from pynurbs.geometry.line import Line
from pynurbs.geometry.methods.tessellate import (adaptive_curve_tessellate,
                                                 adaptive_surface_tessellate)
from pynurbs.geometry.nurbs_curve import NurbsCurve
from pynurbs.geometry.nurbs_surface import NurbsSurface
from pynurbs.geometry.plane import Plane
from pynurbs.geometry.point import Point
from pynurbs.geometry.tessellator import SurfaceTessAdaptive, CurveTessAdaptive
from pynurbs.geometry.vector import Vector


def plot_triangular_mesh(verts, triangles, fig, representation='wireframe',
                         color=None):
    if color is None:
        color = tuple(np.random.rand(1, 3)[0])
    # Generate arrays.
    x = verts[:, 0]
    y = verts[:, 1]
    try:
        z = verts[:, 2]
    except IndexError:
        z = np.zeros(x.shape, dtype=np.float64)
    mlab.triangular_mesh(x, y, z, triangles, figure=fig,
                         representation=representation,
                         color=color)


def plot_geom(geom, fig, ctol=0.001, stol=0.01, scale=0.15):
    """
    Method to plot and show a list of geometry.

    :param geom: Geometry entity to plot.
    :param figure fig: Existing figure object to plot to.
    :param float ctol: Tolerance for adaptive curve tessellation.
    :param float stol: Tolerance for adaptive surface tessellation.
    :param float scale: Option to set scale of point instances.
    """
    if isinstance(geom, Point):
        _plot_point(geom, fig, scale)
    if isinstance(geom, Vector):
        _plot_vector(geom, fig)
    if isinstance(geom, Plane):
        _plot_plane(geom, fig)
    if isinstance(geom, (BezierCurve, NurbsCurve)):
        _plot_curve(geom, fig, ctol)
    if isinstance(geom, ICurve):
        c = geom.crv3d
        c.set_color(*geom.color)
        _plot_curve(c, fig, ctol)
    if isinstance(geom, BezierSurface):
        _plot_bezier_surface(geom, fig, stol)
    if isinstance(geom, NurbsSurface):
        _plot_nurbs_surface(geom, fig, stol)
    if isinstance(geom, SurfaceTessAdaptive):
        _plot_surface_tess(geom, fig)
    if isinstance(geom, CurveTessAdaptive):
        _plot_curve_tess(geom, fig)
    if isinstance(geom, Line):
        _plot_line(geom, fig)


def _plot_point(p, fig, scale=0.15):
    x, y, z = p.xyz
    mlab.points3d(x, y, z, color=p.color, figure=fig, scale_factor=scale)


def _plot_vector(vec, fig):
    x, y, z = vec.origin.xyz
    u, v, w = vec.vxyz
    mlab.quiver3d(x, y, z, u, v, w, color=vec.color, figure=fig, mode='arrow')


def _plot_line(line, fig):
    x, y, z = line.p0.xyz
    u, v, w = line.v.vxyz
    mlab.quiver3d(x, y, z, u, v, w, color=line.color, figure=fig, mode='arrow')


def _plot_plane(plane, fig):
    p0 = plane.p0
    _plot_point(p0, fig)
    x, y, z = plane.p0.xyz
    u, v, w = plane.vn.vxyz
    mlab.quiver3d(x, y, z, u, v, w, color=(0, 0, 1), figure=fig, mode='arrow')

    xyz = np.zeros((5, 3))
    xyz[0] = plane.eval(0.5, 0.5, rtype='ndarray')
    xyz[1] = plane.eval(-0.5, 0.5, rtype='ndarray')
    xyz[2] = plane.eval(-0.5, -0.5, rtype='ndarray')
    xyz[3] = plane.eval(0.5, -0.5, rtype='ndarray')
    xyz[4] = plane.eval(0.5, 0.5, rtype='ndarray')
    mlab.plot3d(xyz[:, 0], xyz[:, 1], xyz[:, 2], figure=fig, color=(0, 0, 0))


def _plot_curve(c, fig, tol):
    xyz, _ = adaptive_curve_tessellate(c, tol)
    mlab.plot3d(xyz[:, 0], xyz[:, 1], xyz[:, 2], figure=fig, color=c.color)


def _plot_bezier_surface(s, fig, tol):
    verts, tris = adaptive_surface_tessellate(s, tol)
    x = verts[:, 0]
    y = verts[:, 1]
    z = verts[:, 2]
    if verts.shape[0] >= 3:
        mlab.triangular_mesh(x, y, z, tris, figure=fig, color=s.color,
                             opacity=s.opacity)
    else:
        warnings.warn('Unable to tessellate surface in graphics module.',
                      RuntimeWarning)


def _plot_nurbs_surface(s, fig, tol):
    verts, tris = adaptive_surface_tessellate(s, tol)
    x = verts[:, 0]
    y = verts[:, 1]
    z = verts[:, 2]
    if verts.shape[0] >= 3:
        mlab.triangular_mesh(x, y, z, tris, figure=fig, color=s.color,
                             opacity=s.opacity)
    else:
        warnings.warn('Unable to tessellate surface in graphics module.',
                      RuntimeWarning)


def _plot_curve_tess(tess, fig):
    verts = tess.verts
    mlab.plot3d(verts[:, 0], verts[:, 1], verts[:, 2], figure=fig,
                color=tess.color)


def _plot_surface_tess(tess, fig):
    verts, triangles = tess.verts, tess.triangles
    x = verts[:, 0]
    y = verts[:, 1]
    z = verts[:, 2]
    mlab.triangular_mesh(x, y, z, triangles, figure=fig,
                         representation=tess.representation,
                         color=tess.color, opacity=tess.opacity)
