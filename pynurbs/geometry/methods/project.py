# TODO Projection methods in Python
from __future__ import division

from numpy import dot
from numpy.linalg import norm


def project_point_to_line(point, line):
    """
    Project a point to a line.

    :param point:
    :param line:

    :return:
    """
    vp = point.xyz - line.p0.xyz
    u = dot(vp, line.v.ijk)
    pp = line.eval(u, rtype='ndarray')
    d = point.dist2pnt(pp)
    return 0, 1, [[pp, u, d]]


def project_point_to_plane(point, plane):
    """
    Project a point to a plane.

    :param point: Point to project.
    :type point: :class:`.Point`
    :param plane: Plane to project point to.
    :type plane: :class:`.Plane`

    :return: Parameters on plane and projected point (u, v, pnt).
    :rtype: tuple
    """
    vu = plane.vu.ijk
    vv = plane.vv.ijk
    p0 = plane.p0.xyz
    pxyz = point.xyz
    u = dot(pxyz - p0, vu) / dot(vu, vu)
    v = dot(pxyz - p0, vv) / dot(vv, vv)
    pnt = p0 + u * vu + v * vv
    d = norm(pnt - pxyz)
    results = [[pnt, (u, v), d]]
    return 0, 1, results

# def project_point_to_curve(point, curve, ends=True, all_pnts=False):
#     """
#     Project a point to a curve using recursive subdivision combined with
#     root finding methods.
#
#     :param Point point: Point to project.
#     :param curve: Curve to project point to.
#     :type curve: BezierCurve or NurbsCurve
#     :param bool ends: If no potential projections are found, this option
#         will project the point to the nearest end point of the curve.
#     :param bool all_pnts: Option to return only the nearest point during the
#         subdivision process (*False*), or all points (*True*).
#
#     :return: Number of projections and results (npts, results).
#     :rtype: tuple
#     """
#     # Global variables.
#     pxyz = point.xyz
#
#     # Fortran method ----------------------------------------------------------
#     n, p, uk, cpw = curve.data
#     lib_proj.project_point_to_curve(pxyz, n, p, uk, cpw, ends, all_pnts)
#     # Access results
#     npts = int(lib_out.npts_)
#     points = array(lib_out.points_[:npts], dtype=float64, copy=True)
#     params = array(lib_out.params1d_[:npts], dtype=float64, copy=True)
#     dist = array(lib_out.dproj_[:npts], dtype=float64, copy=True)
#     results = []
#     for i in range(0, npts):
#         results.append([points[i], params[i], dist[i]])
#     return 0, npts, results


# def project_point_to_icurve(point, icurve, ends=True, all_pnts=False):
#     """
#     Project a point to an intersection curve.
#     """
#     # Global variables.
#     pxyz = point.xyz
#     gtol = Settings.gtol
#     ptol = Settings.ptol
#
#     # Project point to 3-D curve.
#     crv3d = icurve.crv3d
#     nsub, npts, r3d = project_point_to_curve(point, crv3d)
#     if npts == 0:
#         return nsub, npts, r3d
#
#     # Define objective function for point refinment.
#     def _obj(x):
#         pi = ceval(x[0], rtype='ndarray', domain='global')
#         cu = cderiv(x[0], rtype='ndarray', domain='global')
#         return array([dot(cu, pi - pxyz)], dtype=float64)
#
#     # Refine points.
#     ceval = icurve.eval
#     cderiv = icurve.deriv
#     check_param = icurve.check_param
#     results = []
#     for u0 in [r[1] for r in r3d]:
#         # SciPy root finding method.
#         sol = root(_obj, u0, method='hybr')
#         # if not sol.success:
#         # msg = '. '.join(['Unsuccessful solution in curve projection',
#         #                  sol.message])
#         # warn(msg, RuntimeWarning)
#         # continue
#         # Get solution.
#         u = check_param(sol.x[0])
#         p = ceval(u, rtype='ndarray', domain='global')
#         # Add to results only if the point is unique.
#         unique = True
#         for ptest, utest, _ in results:
#             if norm(ptest - p) <= gtol and abs(utest - u) <= ptol:
#                 unique = False
#                 break
#         if not unique:
#             continue
#         d = norm(p - pxyz)
#         results.append([p, u, d])
#
#     # Return results.
#     npts = len(results)
#     if npts > 0:
#         return nsub, npts, results
#
#     # Select nearest endpoint if no candidate curves are found.
#     if ends:
#         u0, u1 = crv3d.a, crv3d.b
#         p0 = ceval(u0, rtype='ndarray', domain='global')
#         p1 = ceval(u1, rtype='ndarray', domain='global')
#         d0 = norm(pxyz - p0)
#         d1 = norm(pxyz - p1)
#         if d0 <= d1:
#             return 0, 1, [[p0, u0, d0]]
#         else:
#             return 0, 1, [[p1, u1, d1]]


# def project_point_to_surface(point, surface, ends=True, all_pnts=False):
#     """
#     Project a point to a surface using recursive subdivision combined with
#     root finding methods.
#
#     :param Point point: Point to project.
#     :param surface: Surface to project point to.
#     :type surface: BezierSurface or NurbsSurface
#     :param bool ends: Option to project point to nearest boundary curve
#         if no orthogonal projection is found within the initial surface.
#     :param bool all_pnts: Option to return only the nearest point during the
#         subdivision process (*False*), or all points (*True*).
#
#     :return: Number of subdivisions, projections and results
#         (nsub, npts, results).
#     :rtype: tuple
#     """
#
#     # Global variables.
#     pxyz = point.xyz
#
#     # Fortran method ----------------------------------------------------------
#     n, p, uk, m, q, vk, cpw = surface.data
#     lib_proj.project_point_to_surface(pxyz, n, p, uk, m, q, vk, cpw,
#                                       ends, all_pnts)
#     # Access results.
#     npts = int(lib_out.npts_)
#     points = array(lib_out.points_[:npts], dtype=float64, copy=True)
#     params = array(lib_out.params2d_[:npts], dtype=float64, copy=True)
#     dist = array(lib_out.dproj_[:npts], dtype=float64, copy=True)
#     results = []
#     for i in range(0, npts):
#         results.append([points[i], tuple(params[i]), dist[i]])
#     return 0, npts, results
