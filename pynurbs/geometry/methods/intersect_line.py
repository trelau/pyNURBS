from __future__ import division

from numpy import array, cross, dot, float64
from numpy.linalg import norm
from scipy.optimize import minimize

from pynurbs.config import Settings
from pynurbs.geometry.methods.geom_utils import (is_surface_flat,
                                                 invert_point_in_triangle)
from pynurbs.geometry.methods.intersect_bbox import bbox_ray_intersect


def intersect_line_line(line1, line2, itol=None):
    """
    Intersect two lines.
    """
    pi1, pi2 = line1.p0.xyz, line1.p0.xyz + line1.v.ijk
    pi3, pi4 = line2.p0.xyz, line2.p0.xyz + line2.v.ijk
    d4321 = dot(pi4 - pi3, pi2 - pi1)
    d1321 = dot(pi1 - pi3, pi2 - pi1)
    d4343 = dot(pi4 - pi3, pi4 - pi3)
    d2121 = dot(pi2 - pi1, pi2 - pi1)
    denom = d2121 * d4343 - d4321 * d4321
    if abs(denom) <= 1.0e-12:
        return 0, 0, []
    d1343 = dot(pi1 - pi3, pi4 - pi3)
    numer = d1343 * d4321 - d1321 * d4343
    # Parameters on each line.
    u12 = numer / denom
    u34 = (d1343 + u12 * d4321) / d4343
    # Point on each line.
    p1 = line1.eval(u12, rtype='ndarray')
    p2 = line2.eval(u34, rtype='ndarray')
    d = norm(p1 - p2)
    if d > itol:
        return 0, 0, []
    pi = 0.5 * (p1 + p2)
    return 0, 1, [[(u12, u34), pi]]


def intersect_line_plane(line, plane, itol=None, t0=None, t1=None):
    """
    Intersect a line and a plane.
    """
    pnorm = plane.vn.ijk
    vline = line.v.ijk
    denom = dot(pnorm, vline)
    # Check for parellel line.
    if abs(denom) <= 1.0e-12:
        return 0, 0, []
    if itol is None:
        itol = Settings.gtol
    p0_plane = plane.p0.xyz
    p0_line = line.p0.xyz
    ui = dot(pnorm, p0_plane - p0_line)
    pi = line.eval(ui, rtype='ndarray')
    if t0 is not None:
        plow = line.eval(t0, rtype='ndarray')
        if ui < t0 and norm(pi - plow) > itol:
            return 0, 0, []
    if t1 is not None:
        phigh = line.eval(t1, rtype='ndarray')
        if ui > t1 and norm(pi - phigh) > itol:
            return 0, 0, []
    return 1, 1, [[ui, pi]]


def intersect_line_surface(line, surface, itol=None, t0=None, t1=None):
    """
    Find the intersection points between a line and a surface.

    :param line:
    :param surface:
    :param itol:
    :param t0:
    :param t1:

    :return:
    """
    # Global parameters.
    ftol = Settings.ftol
    if itol is None:
        itol = Settings.gtol

    # Get line origin and direction vector.
    p0 = line.p0.xyz
    vline = line.v.vxyz
    vline /= norm(vline)

    # Step 1: Checking box-ray intersection for possible intersection.
    bbox = surface.get_bbox()
    if not bbox_ray_intersect(bbox, p0, vline, t0, t1, itol):
        return 0, 0, []

    # Step 2: Decompose each into Bezier segments.
    bezier_patches = []
    _, _, surfs = surface.decompose()
    for row in surfs:
        for s in row:
            bezier_patches.append(s)

    # Step 3: Define methods for recursive subdivision.
    def _subdivide(si):
        """
        Recursive subdivision to find possible intersections.
        """
        nsub[0] += 1
        cpi = si.cp
        is_flat = _is_flat(si.n, si.m, cpi)
        if is_flat:
            _intersect(si, cpi)
        else:
            si1234 = si.split(0.5, 0.5)
            for sii in si1234:
                if _candidate_intersect(sii):
                    _subdivide(sii)

    def _is_flat(ni, mi, cpi):
        """
        Check surface flatness.
        """
        return is_surface_flat(ni, mi, cpi, ftol)

    def _candidate_intersect(si):
        """
        Check for potential intersection.
        """
        return bbox_ray_intersect(si.get_bbox(), p0, vline, t0, t1, itol)

    def _intersect(si, cpi):
        """
        Intersect line and plane.
        """
        ps0 = cpi[0, 0]
        vu, vv = cpi[-1, 0] - ps0, cpi[0, -1] - ps0
        if abs(norm(vu)) <= 1.0e-12 or abs(norm(vv)) <= 1.0e-12:
            return None
        pnorm = cross(vu, vv)
        pnorm /= norm(pnorm)
        denom = dot(pnorm, vline)
        if abs(denom) <= 1.0e-12:
            # If parallel use a starting point at midpoints.
            ti = line.local_to_global_param(0.5)
            ui = si.local_to_global_param('u', 0.5)
            vi = si.local_to_global_param('v', 0.5)
            tuv0_list.append((ti, ui, vi))
        else:
            ti = dot(pnorm, ps0 - p0) / denom
            pi = p0 + ti * vline
            tri3d = array([cpi[0, 0], cpi[-1, 0], cpi[0, -1]], dtype=float64)
            tri2d = array([[si.au, si.av],
                           [si.bu, si.av],
                           [si.au, si.bv]], dtype=float64)
            ui, vi = invert_point_in_triangle(pi, tri3d, tri2d, False)
            tuv0_list.append((ti, ui, vi))

    # Step 4: Perform recursive subdivision.
    tuv0_list = []
    nsub = [0]
    for sbez in bezier_patches:
        if _candidate_intersect(sbez):
            _subdivide(sbez)

    # Step 6: Refine points.
    def _obj(x):
        factor = 1.
        # Penalize objective function if variable is outside domain.
        if x[1] < au or x[1] > bu:
            factor = 1000.
        if x[2] < av or x[2] > bv:
            factor = 1000.
        pci = ceval(x[0], rtype='ndarray', domain='global')
        psi = seval(x[1], x[2], rtype='ndarray', domain='global')
        return norm(pci - psi) * factor

    ceval = line.eval
    seval = surface.eval
    results = []
    au, bu, av, bv = surface.au, surface.bu, surface.av, surface.bv
    tol = itol / 100.
    for t0, u0, v0 in tuv0_list:
        tuv0 = array([t0, u0, v0], dtype=float64)
        sol = minimize(_obj, tuv0, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        t, u, v = sol.x
        u, v = surface.check_params(u, v)
        pc = ceval(t, rtype='ndarray', domain='global')
        ps = seval(u, v, rtype='ndarray', domain='global')
        if norm(pc - ps) > itol:
            continue
        unique = True
        for _, p12 in results:
            if norm(p12 - pc) <= itol:
                unique = False
                break
        if not unique:
            continue
        results.append([(t, u, v), pc])

    npts = len(results)
    return nsub[0], npts, results
