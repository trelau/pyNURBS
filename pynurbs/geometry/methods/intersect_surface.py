from __future__ import division, print_function

from math import ceil

from numpy import array, cross, dot, float64, int32, zeros, mean
from numpy.linalg import norm
from scipy.optimize import minimize

from pynurbs.config import Settings
from pynurbs.geometry.methods.calculate import triangle_area
from pynurbs.geometry.methods.geom_utils import (barycentric_params,
                                                 is_surface_flat,
                                                 angle_between_vecs)
from pynurbs.geometry.methods.intersect_bbox import (bboxes_intersect,
                                                     bbox_intersects_plane)
from pynurbs.geometry.methods.intersect_triangle import (
    intersect_triangle_plane, intersect_triangles)
from pynurbs.geometry.methods.invert import invert_points_on_plane
from pynurbs.geometry.methods.tessellate import tessellate_cell


class _Cell(object):
    """
    Surface cell for Python surface intersection methods.
    """

    def __init__(self):
        self.sid = 0
        self.position = 0
        self.is_cand = False
        self.has_child = False
        self.u0 = 0.
        self.u1 = 0.
        self.v0 = 0.
        self.v1 = 0.
        self.parent = None
        self.ne = None
        self.se = None
        self.sw = None
        self.nw = None
        self.n = None
        self.e = None
        self.s = None
        self.w = None


def intersect_surface_plane(surface, plane, ftol):
    """
    Find the intersection curve(s) between a surface and a plane.

    :param surface: Surface to intersect.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param Plane plane: Intersection plane.
    :param float ftol: Surface flatness tolerance.

    :return: Surface intersection results.
    :rtype: tuple
    """
    # Global parameters.
    gtol = Settings.gtol
    ptol = Settings.ptol
    p0 = plane.p0.xyz
    pnorm = plane.vn.ijk
    vx = plane.vu.vxyz
    vy = plane.vv.vxyz

    # Parameters for subdivision.
    cell_list = []
    ncells = [0]
    nsub = [0]
    ncand = [0]

    # Step 1: Define methods for recursive subdivision.
    def _subdivide(si, ci):
        """
        Recursive subdivision for potential intersection segments.
        """
        nsub[0] += 1

        # Store surface parameters.
        ci.u0 = si.au
        ci.u1 = si.bu
        ci.v0 = si.av
        ci.v1 = si.bv

        # Store cell.
        cell_list.append(ci)

        # Check flatness.
        cpi = si.cp
        is_flat = is_surface_flat(si.n, si.m, cpi, ftol)
        if is_flat and nsub[0] > 1:
            # Save candidate surfaces.
            ci.is_cand = True
            ncand[0] += 1
        else:
            # Split surface into four patches.
            ui = 0.5 * (ci.u0 + ci.u1)
            vi = 0.5 * (ci.v0 + ci.v1)
            if si.p == 1 and si.uk.size > 4:
                mid = int(ceil((si.uk.size - 1) / 2))
                ui = si.uk[mid]
            if si.q == 1 and si.vk.size > 4:
                mid = int(ceil((si.vk.size - 1) / 2))
                vi = si.vk[mid]
            s1, s2, s3, s4 = si.split(ui, vi, domain='global')

            # Check potential intersection
            t1, t2, t3, t4 = map(_candidate_intersect, [s1, s2, s3, s4])

            # New cells.
            c_ne, c_se, c_sw, c_nw = ci.ne, ci.se, ci.sw, ci.nw

            # Assign children (if any) in case cells are duplicated.
            if c_ne is None:
                c_ne = _Cell()
            if c_se is None:
                c_se = _Cell()
            if c_sw is None:
                c_sw = _Cell()
            if c_nw is None:
                c_nw = _Cell()

            if t1:
                # Surface 1: Cell 1 (SW)
                if c_sw.sid == 0:
                    c_sw.sid = ncells[0] + 1
                    ncells[0] += 1
                ci.sw = c_sw
                c_sw.parent = ci
                c_sw.position = 1
                ci.has_child = True

            if t2:
                # Surface 1: Cell 2 (NW)
                if c_nw.sid == 0:
                    c_nw.sid = ncells[0] + 1
                    ncells[0] += 1
                ci.nw = c_nw
                c_nw.parent = ci
                c_nw.position = 2
                ci.has_child = True

            if t3:
                # Surface 1: Cell 3 (SE)
                if c_se.sid == 0:
                    c_se.sid = ncells[0] + 1
                    ncells[0] += 1
                ci.se = c_se
                c_se.parent = ci
                c_se.position = 3
                ci.has_child = True

            if t4:
                # Surface 1: Cell 4 (NE)
                if c_ne.sid == 0:
                    c_ne.sid = ncells[0] + 1
                    ncells[0] += 1
                ci.ne = c_ne
                c_ne.parent = ci
                c_ne.position = 4
                ci.has_child = True

            # Adjacent cells.
            if t1:
                # Surface 1: Cell 1 (SW)
                c_sw.n = c_nw
                c_sw.e = c_se
                c_sw.s = ci.s
                c_sw.w = ci.w

            if t2:
                # Surface 1: Cell 2 (NW)
                c_nw.n = ci.n
                c_nw.e = c_ne
                c_nw.s = c_sw
                c_nw.w = ci.w

            if t3:
                # Surface 1: Cell 3 (SW)
                c_se.n = c_ne
                c_se.e = ci.e
                c_se.s = ci.s
                c_se.w = c_sw

            if t4:
                # Surface 1: Cell 4 (NE)
                c_ne.n = ci.n
                c_ne.e = ci.e
                c_ne.s = c_se
                c_ne.w = c_nw

            # Subdivide.
            if t1:
                _subdivide(s1, c_sw)
            if t2:
                _subdivide(s2, c_nw)
            if t3:
                _subdivide(s3, c_se)
            if t4:
                _subdivide(s4, c_ne)

    def _candidate_intersect(si):
        """
        Check for possible intersection using bounding box.
        """
        bbox = si.get_bbox()
        return bbox_intersects_plane(bbox, plane, gtol)

    # Step 2: Use recursive subdivision to find candidate surfaces.
    c0 = _Cell()
    if _candidate_intersect(surface):
        _subdivide(surface, c0)

    # Check for no intersection.
    if ncand[0] == 0:
        return 0, [], [], [], []

    children, parent, acells, position, candidates, cell_params = \
        _build_tess_data(cell_list)

    # Methods for intersection.
    def _intersect(icell):
        """
        Intersect cell using triangle-plane intersection.
        """
        # Get parameters of tessellated cell.
        ntri, triangles = tessellate_cell(icell, children, acells, position,
                                          parent, cell_params)
        # Intersect each triangle with the plane.
        for ii in range(0, ntri):
            uv0, uv1, uv2 = triangles[ii]
            ti[0] = seval(*uv0, rtype='ndarray', domain='global')
            ti[1] = seval(*uv1, rtype='ndarray', domain='global')
            ti[2] = seval(*uv2, rtype='ndarray', domain='global')
            # Check for degenerate triangle.
            area = triangle_area(ti)
            if area <= 1.0e-12:
                continue
            ni, pi = intersect_triangle_plane(ti, p0, pnorm, gtol)
            if ni == 2:
                edges.append([ptotal[0], ptotal[0] + 1])
                ptotal[0] += 2
                etotal[0] += 1
                points.append(pi[0])
                points.append(pi[1])
                # Get parameters.
                # Surface
                ub, vb = barycentric_params(pi[0], ti)
                ur, vr = (1. - ub - vb) * uv0 + ub * uv1 + vb * uv2
                params1.append([ur, vr])
                points2d_s1.append([ur, vr, 0.])
                ub, vb = barycentric_params(pi[1], ti)
                ur, vr = (1. - ub - vb) * uv0 + ub * uv1 + vb * uv2
                params1.append([ur, vr])
                points2d_s1.append([ur, vr, 0.])

                # Viewer.add_items(surface.eval(ur, vr, domain='global'))

                # Plane
                uv_plane = invert_points_on_plane(pi, plane)
                params2.append([uv_plane[0, 0], uv_plane[0, 1]])
                params2.append([uv_plane[1, 0], uv_plane[1, 1]])
                points2d_s2.append([uv_plane[0, 0], uv_plane[0, 1], 0.])
                points2d_s2.append([uv_plane[1, 0], uv_plane[1, 1], 0.])

    # Tessellate potential surfaces and intersect.
    ptotal = [0]
    etotal = [0]
    points = []
    points2d_s1 = []
    points2d_s2 = []
    edges = []
    params1 = []
    params2 = []
    seval = surface.eval
    ti = zeros((3, 3), dtype=float64)
    for i in range(0, ncells[0] + 1):
        if candidates[i] == 1:
            _intersect(i)

    # Check for no intersections.
    if ptotal[0] <= 1:
        return 0, [], [], [], []

    # Find the topology of the intersection curves.
    arr_points = array(points, dtype=float64)
    arr_params1 = array(params1, dtype=float64)
    arr_params2 = array(params2, dtype=float64)
    arr_points2d_s1 = array(points2d_s1, dtype=float64)
    arr_points2d_s2 = array(points2d_s2, dtype=float64)
    ncrvs, all_crv_ids = _trace_curves(ptotal[0], etotal[0], arr_points2d_s1,
                                       arr_points2d_s2, edges, ptol)
    # Refine intersection points.
    tol = gtol / 100.
    for i in range(ncrvs):
        for indx in all_crv_ids[i]:
            u, v = arr_params1[indx]
            u, v, up, vp, p = refine_spi_point(surface, plane, u, v, tol)
            arr_params1[indx, :] = [u, v]
            arr_params2[indx, :] = [up, vp]
            arr_points[indx, :] = p

    # Build a list that specifies if a point should not be filtered.
    # no_filter = []
    # get_mult = surface.get_mult
    # p = surface.p
    # q = surface.q
    # for u, v in arr_params1:
    #     if p >= get_mult(u, 'u') or q >= get_mult(v, 'v'):
    #         no_filter.append(True)
    #     else:
    #         no_filter.append(False)

    # Filter out points.
    # no_filter = [False] * arr_points.shape[0]
    crv_size, crv_ids = _filter_points(ncrvs, all_crv_ids, arr_points, gtol)
    # crv_ids = all_crv_ids
    # crv_size = [len(row) for row in crv_ids]

    # Return results.
    return ncrvs, crv_size, crv_ids, arr_points, arr_params1, arr_params2


def intersect_surface_surface(surface1, surface2, ftol):
    """
    Find the intersection curve(s) between two surfaces.

    :param surface1: Surface 1 to intersect.
    :type surface1: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param surface2: Surface 2 to intersect.
    :type surface1: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float ftol: Surface flatness tolerance.

    :return: Surface intersection results.
    :rtype: tuple
    """
    # Global parameters
    gtol = Settings.gtol
    ptol = Settings.ptol

    cell_list1 = []
    cell_list2 = []
    ssi_list = []
    nsub = [0]
    ncand = [0]
    ncells1 = [0]
    ncells2 = [0]

    def _subdivide(si1, si2, ci1, ci2):
        """
        Recursive subdivision for potential intersection segments.
        """
        nsub[0] += 1

        # Store surface parameters.
        # Surface 1
        ci1.u0 = si1.au
        ci1.u1 = si1.bu
        ci1.v0 = si1.av
        ci1.v1 = si1.bv

        # Surface 2
        ci2.u0 = si2.au
        ci2.u1 = si2.bu
        ci2.v0 = si2.av
        ci2.v1 = si2.bv

        # Store cells
        cell_list1.append(ci1)
        cell_list2.append(ci2)

        # Check flatness.
        cpi1 = si1.cp
        cpi2 = si2.cp
        is_flat1 = _is_flat(si1.n, si1.m, cpi1)
        is_flat2 = _is_flat(si2.n, si2.m, cpi2)
        if is_flat1 and is_flat2 and nsub[0] > 1:
            # Save candidate surfaces.
            ci1.is_cand = True
            ci2.is_cand = True
            ssi_list.append([ci1, ci2])
            ncand[0] += 1
        else:
            if not is_flat1 and is_flat2:
                # Subdivide surface 1.
                ui = 0.5 * (ci1.u0 + ci1.u1)
                vi = 0.5 * (ci1.v0 + ci1.v1)
                # If linear, use middle interior knot.
                if si1.p == 1 and si1.uk.size > 4:
                    mid = int(ceil((si1.uk.size - 1) / 2))
                    ui = si1.uk[mid]
                if si1.q == 1 and si1.vk.size > 4:
                    mid = int(ceil((si1.vk.size - 1) / 2))
                    vi = si1.vk[mid]
                s11, s12, s13, s14 = si1.split(ui, vi, domain='global')

                # Check potential intersection
                t11 = _candidate_intersect(s11, si2)
                t12 = _candidate_intersect(s12, si2)
                t13 = _candidate_intersect(s13, si2)
                t14 = _candidate_intersect(s14, si2)

                # New cells
                c1_ne, c1_se, c1_sw, c1_nw = ci1.ne, ci1.se, ci1.sw, ci1.nw

                # Assign children (if any) in case cells are duplicated.
                if c1_ne is None:
                    c1_ne = _Cell()
                if c1_se is None:
                    c1_se = _Cell()
                if c1_sw is None:
                    c1_sw = _Cell()
                if c1_nw is None:
                    c1_nw = _Cell()

                # Assign children
                if t11:
                    # Surface 1: Cell 1 (SW)
                    if c1_sw.sid == 0:
                        c1_sw.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.sw = c1_sw
                    c1_sw.parent = ci1
                    c1_sw.position = 1
                    ci1.has_child = True

                if t12:
                    # Surface 1: Cell 2 (NW)
                    if c1_nw.sid == 0:
                        c1_nw.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.nw = c1_nw
                    c1_nw.parent = ci1
                    c1_nw.position = 2
                    ci1.has_child = True

                if t13:
                    # Surface 1: Cell 3 (SE)
                    if c1_se.sid == 0:
                        c1_se.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.se = c1_se
                    c1_se.parent = ci1
                    c1_se.position = 3
                    ci1.has_child = True

                if t14:
                    # Surface 1: Cell 4 (NE)
                    if c1_ne.sid == 0:
                        c1_ne.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.ne = c1_ne
                    c1_ne.parent = ci1
                    c1_ne.position = 4
                    ci1.has_child = True

                # Adjacent cells
                if t11:
                    # Surface 1: Cell 1 (SW)
                    c1_sw.n = c1_nw
                    c1_sw.e = c1_se
                    c1_sw.s = ci1.s
                    c1_sw.w = ci1.w

                if t12:
                    # Surface 1: Cell 2 (NW)
                    c1_nw.n = ci1.n
                    c1_nw.e = c1_ne
                    c1_nw.s = c1_sw
                    c1_nw.w = ci1.w

                if t13:
                    # Surface 1: Cell 3 (SW)
                    c1_se.n = c1_ne
                    c1_se.e = ci1.e
                    c1_se.s = ci1.s
                    c1_se.w = c1_sw

                if t14:
                    # Surface 1: Cell 4 (NE)
                    c1_ne.n = ci1.n
                    c1_ne.e = ci1.e
                    c1_ne.s = c1_se
                    c1_ne.w = c1_nw

                # Subdivide.
                if t11:
                    _subdivide(s11, si2, c1_sw, ci2)
                if t12:
                    _subdivide(s12, si2, c1_nw, ci2)
                if t13:
                    _subdivide(s13, si2, c1_se, ci2)
                if t14:
                    _subdivide(s14, si2, c1_ne, ci2)

            elif is_flat1 and not is_flat2:
                # Subdivide surface 2.
                ui = 0.5 * (ci2.u0 + ci2.u1)
                vi = 0.5 * (ci2.v0 + ci2.v1)
                # If linear, use middle interior knot.
                if si2.p == 1 and si2.uk.size > 4:
                    mid = int(ceil((si2.uk.size - 1) / 2))
                    ui = si2.uk[mid]
                if si2.q == 1 and si2.vk.size > 4:
                    mid = int(ceil((si2.vk.size - 1) / 2))
                    vi = si2.vk[mid]
                s21, s22, s23, s24 = si2.split(ui, vi, domain='global')

                # Check potential intersection
                t21 = _candidate_intersect(si1, s21)
                t22 = _candidate_intersect(si1, s22)
                t23 = _candidate_intersect(si1, s23)
                t24 = _candidate_intersect(si1, s24)

                # New cells
                c2_ne, c2_se, c2_sw, c2_nw = ci2.ne, ci2.se, ci2.sw, ci2.nw

                # Assign children (if any) in case cells are duplicated.
                if c2_ne is None:
                    c2_ne = _Cell()
                if c2_se is None:
                    c2_se = _Cell()
                if c2_sw is None:
                    c2_sw = _Cell()
                if c2_nw is None:
                    c2_nw = _Cell()

                # Assign children
                if t21:
                    # Surface 2: Cell 1 (SW)
                    if c2_sw.sid == 0:
                        c2_sw.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.sw = c2_sw
                    c2_sw.parent = ci2
                    c2_sw.position = 1
                    ci2.has_child = True

                if t22:
                    # Surface 2: Cell 2 (NW)
                    if c2_nw.sid == 0:
                        c2_nw.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.nw = c2_nw
                    c2_nw.parent = ci2
                    c2_nw.position = 2
                    ci2.has_child = True

                if t23:
                    # Surface 2: Cell 3 (SE)
                    if c2_se.sid == 0:
                        c2_se.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.se = c2_se
                    c2_se.parent = ci2
                    c2_se.position = 3
                    ci2.has_child = True

                if t24:
                    # Surface 2: Cell 4 (NE)
                    if c2_ne.sid == 0:
                        c2_ne.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.ne = c2_ne
                    c2_ne.parent = ci2
                    c2_ne.position = 4
                    ci2.has_child = True

                # Adjacent cells
                if t21:
                    # Surface 2: Cell 1 (SW)
                    c2_sw.n = c2_nw
                    c2_sw.e = c2_se
                    c2_sw.s = ci2.s
                    c2_sw.w = ci2.w

                if t22:
                    # Surface 2: Cell 2 (NW)
                    c2_nw.n = ci2.n
                    c2_nw.e = c2_ne
                    c2_nw.s = c2_sw
                    c2_nw.w = ci2.w

                if t23:
                    # Surface 2: Cell 3 (SW)
                    c2_se.n = c2_ne
                    c2_se.e = ci2.e
                    c2_se.s = ci2.s
                    c2_se.w = c2_sw

                if t24:
                    # Surface 2: Cell 4 (NE)
                    c2_ne.n = ci2.n
                    c2_ne.e = ci2.e
                    c2_ne.s = c2_se
                    c2_ne.w = c2_nw

                # Subdivide
                if t21:
                    _subdivide(si1, s21, ci1, c2_sw)
                if t22:
                    _subdivide(si1, s22, ci1, c2_nw)
                if t23:
                    _subdivide(si1, s23, ci1, c2_se)
                if t24:
                    _subdivide(si1, s24, ci1, c2_ne)

            else:
                # Subdivide each surface into four patches.
                # Surface 1.
                ui = 0.5 * (ci1.u0 + ci1.u1)
                vi = 0.5 * (ci1.v0 + ci1.v1)
                # If linear, use middle interior knot.
                if si1.p == 1 and si1.uk.size > 4:
                    mid = int(ceil((si1.uk.size - 1) / 2))
                    ui = si1.uk[mid]
                if si1.q == 1 and si1.vk.size > 4:
                    mid = int(ceil((si1.vk.size - 1) / 2))
                    vi = si1.vk[mid]
                s11, s12, s13, s14 = si1.split(ui, vi, domain='global')

                # Surface 2.
                ui = 0.5 * (ci2.u0 + ci2.u1)
                vi = 0.5 * (ci2.v0 + ci2.v1)
                # If linear, use middle interior knot.
                if si2.p == 1 and si2.uk.size > 4:
                    mid = int(ceil((si2.uk.size - 1) / 2))
                    ui = si2.uk[mid]
                if si2.q == 1 and si2.vk.size > 4:
                    mid = int(ceil((si2.vk.size - 1) / 2))
                    vi = si2.vk[mid]
                s21, s22, s23, s24 = si2.split(ui, vi, domain='global')

                # Check potential intersection
                t11_t21 = _candidate_intersect(s11, s21)
                t11_t22 = _candidate_intersect(s11, s22)
                t11_t23 = _candidate_intersect(s11, s23)
                t11_t24 = _candidate_intersect(s11, s24)
                t12_t21 = _candidate_intersect(s12, s21)
                t12_t22 = _candidate_intersect(s12, s22)
                t12_t23 = _candidate_intersect(s12, s23)
                t12_t24 = _candidate_intersect(s12, s24)
                t13_t21 = _candidate_intersect(s13, s21)
                t13_t22 = _candidate_intersect(s13, s22)
                t13_t23 = _candidate_intersect(s13, s23)
                t13_t24 = _candidate_intersect(s13, s24)
                t14_t21 = _candidate_intersect(s14, s21)
                t14_t22 = _candidate_intersect(s14, s22)
                t14_t23 = _candidate_intersect(s14, s23)
                t14_t24 = _candidate_intersect(s14, s24)

                # New cells
                c1_ne, c1_se, c1_sw, c1_nw = ci1.ne, ci1.se, ci1.sw, ci1.nw
                c2_ne, c2_se, c2_sw, c2_nw = ci2.ne, ci2.se, ci2.sw, ci2.nw

                # Assign children (if any) in case cells are duplicated.
                if c1_ne is None:
                    c1_ne = _Cell()
                if c1_se is None:
                    c1_se = _Cell()
                if c1_sw is None:
                    c1_sw = _Cell()
                if c1_nw is None:
                    c1_nw = _Cell()

                if c2_ne is None:
                    c2_ne = _Cell()
                if c2_se is None:
                    c2_se = _Cell()
                if c2_sw is None:
                    c2_sw = _Cell()
                if c2_nw is None:
                    c2_nw = _Cell()

                # Assign children
                if t11_t21 or t11_t22 or t11_t23 or t11_t24:
                    # Surface 1: Cell 1 (SW)
                    if c1_sw.sid == 0:
                        c1_sw.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.sw = c1_sw
                    c1_sw.parent = ci1
                    c1_sw.position = 1
                    ci1.has_child = True

                if t12_t21 or t12_t22 or t12_t23 or t12_t24:
                    # Surface 1: Cell 2 (NW)
                    if c1_nw.sid == 0:
                        c1_nw.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.nw = c1_nw
                    c1_nw.parent = ci1
                    c1_nw.position = 2
                    ci1.has_child = True

                if t13_t21 or t13_t22 or t13_t23 or t13_t24:
                    # Surface 1: Cell 3 (SE)
                    if c1_se.sid == 0:
                        c1_se.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.se = c1_se
                    c1_se.parent = ci1
                    c1_se.position = 3
                    ci1.has_child = True

                if t14_t21 or t14_t22 or t14_t23 or t14_t24:
                    # Surface 1: Cell 4 (NE)
                    if c1_ne.sid == 0:
                        c1_ne.sid = ncells1[0] + 1
                        ncells1[0] += 1
                    ci1.ne = c1_ne
                    c1_ne.parent = ci1
                    c1_ne.position = 4
                    ci1.has_child = True

                if t11_t21 or t12_t21 or t13_t21 or t14_t21:
                    # Surface 2: Cell 1 (SW)
                    if c2_sw.sid == 0:
                        c2_sw.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.sw = c2_sw
                    c2_sw.parent = ci2
                    c2_sw.position = 1
                    ci2.has_child = True

                if t11_t22 or t12_t22 or t13_t22 or t14_t22:
                    # Surface 2: Cell 2
                    if c2_nw.sid == 0:
                        c2_nw.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.nw = c2_nw
                    c2_nw.parent = ci2
                    c2_nw.position = 2
                    ci2.has_child = True

                if t11_t23 or t12_t23 or t13_t23 or t14_t23:
                    # Surface 2: Cell 3 (SE)
                    if c2_se.sid == 0:
                        c2_se.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.se = c2_se
                    c2_se.parent = ci2
                    c2_se.position = 3
                    ci2.has_child = True

                if t11_t24 or t12_t24 or t13_t24 or t14_t24:
                    # Surface 2: Cell 4
                    if c2_ne.sid == 0:
                        c2_ne.sid = ncells2[0] + 1
                        ncells2[0] += 1
                    ci2.ne = c2_ne
                    c2_ne.parent = ci2
                    c2_ne.position = 4
                    ci2.has_child = True

                # Adjacent cells
                if t11_t21 or t11_t22 or t11_t23 or t11_t24:
                    # Surface 1: Cell 1 (SW)
                    c1_sw.n = c1_nw
                    c1_sw.e = c1_se
                    c1_sw.s = ci1.s
                    c1_sw.w = ci1.w

                if t12_t21 or t12_t22 or t12_t23 or t12_t24:
                    # Surface 1: Cell 2 (NW)
                    c1_nw.n = ci1.n
                    c1_nw.e = c1_ne
                    c1_nw.s = c1_sw
                    c1_nw.w = ci1.w

                if t13_t21 or t13_t22 or t13_t23 or t13_t24:
                    # Surface 1: Cell 3 (SW)
                    c1_se.n = c1_ne
                    c1_se.e = ci1.e
                    c1_se.s = ci1.s
                    c1_se.w = c1_sw

                if t14_t21 or t14_t22 or t14_t23 or t14_t24:
                    # Surface 1: Cell 4 (NE)
                    c1_ne.n = ci1.n
                    c1_ne.e = ci1.e
                    c1_ne.s = c1_se
                    c1_ne.w = c1_nw

                if t11_t21 or t12_t21 or t13_t21 or t14_t21:
                    # Surface 2: Cell 1 (SW)
                    c2_sw.n = c2_nw
                    c2_sw.e = c2_se
                    c2_sw.s = ci2.s
                    c2_sw.w = ci2.w

                if t11_t22 or t12_t22 or t13_t22 or t14_t22:
                    # Surface 2: Cell 2 (NW)
                    c2_nw.n = ci2.n
                    c2_nw.e = c2_ne
                    c2_nw.s = c2_sw
                    c2_nw.w = ci2.w

                if t11_t23 or t12_t23 or t13_t23 or t14_t23:
                    # Surface 2: Cell 3 (SE)
                    c2_se.n = c2_ne
                    c2_se.e = ci2.e
                    c2_se.s = ci2.s
                    c2_se.w = c2_sw

                if t11_t24 or t12_t24 or t13_t24 or t14_t24:
                    # Surface 2: Cell 4 (NW)
                    c2_ne.n = ci2.n
                    c2_ne.e = ci2.e
                    c2_ne.s = c2_se
                    c2_ne.w = c2_nw

                # Subdivide.
                if t11_t21:
                    _subdivide(s11, s21, c1_sw, c2_sw)
                if t11_t22:
                    _subdivide(s11, s22, c1_sw, c2_nw)
                if t11_t23:
                    _subdivide(s11, s23, c1_sw, c2_se)
                if t11_t24:
                    _subdivide(s11, s24, c1_sw, c2_ne)

                if t12_t21:
                    _subdivide(s12, s21, c1_nw, c2_sw)
                if t12_t22:
                    _subdivide(s12, s22, c1_nw, c2_nw)
                if t12_t23:
                    _subdivide(s12, s23, c1_nw, c2_se)
                if t12_t24:
                    _subdivide(s12, s24, c1_nw, c2_ne)

                if t13_t21:
                    _subdivide(s13, s21, c1_se, c2_sw)
                if t13_t22:
                    _subdivide(s13, s22, c1_se, c2_nw)
                if t13_t23:
                    _subdivide(s13, s23, c1_se, c2_se)
                if t13_t24:
                    _subdivide(s13, s24, c1_se, c2_ne)

                if t14_t21:
                    _subdivide(s14, s21, c1_ne, c2_sw)
                if t14_t22:
                    _subdivide(s14, s22, c1_ne, c2_nw)
                if t14_t23:
                    _subdivide(s14, s23, c1_ne, c2_se)
                if t14_t24:
                    _subdivide(s14, s24, c1_ne, c2_ne)

    def _candidate_intersect(si1, si2):
        """
        Check for possible intersection using bounding box.
        """
        bb1 = si1.get_bbox()
        bb2 = si2.get_bbox()
        return bboxes_intersect(bb1, bb2, gtol)

    def _is_flat(ni, mi, cpi):
        """
        Check surface flatness.
        """
        return is_surface_flat(ni, mi, cpi, ftol)

    # Use recursive subdivision to find candidate surfaces.
    c1 = _Cell()
    c2 = _Cell()
    if _candidate_intersect(surface1, surface2):
        _subdivide(surface1, surface2, c1, c2)

    # Check for no intersection.
    if ncand[0] == 0:
        return 0, [], [], [], [], []

    # Build arrays to tessellate data.
    children1, parent1, acells1, position1, candidates1, cell_params1 = \
        _build_tess_data(cell_list1)

    children2, parent2, acells2, position2, candidates2, cell_params2 = \
        _build_tess_data(cell_list2)

    # Build candidate arrays for intersecting cells.
    candidates = zeros((nsub[0], nsub[0]), dtype=int32)
    for c1, c2 in ssi_list:
        candidates[c1.sid, c2.sid] = 1

    # It's possible that a parent surface and its children may be in the
    # potential intersection list. This can sometimes cause issues in the
    # tessellation algorithm leading to non-congruent edges. For this reason,
    # if the parent of a surface is flat and also in the potential intersection
    # list, replace the child surface with its parent.
    for i in range(0, ncells1[0] + 1):
        for j in range(0, ncells2[0] + 1):
            if candidates[i, j] == 1:
                # Erase intersection.
                candidates[i, j] = 0
                cell1, cell2 = i, j

                # Surface 1.
                while candidates1[parent1[cell1]] == 1:
                    cell1 = parent1[cell1]
                    children1[cell1, :] = 0
                    if cell1 == 0:
                        break

                # Surface 2.
                while candidates2[parent2[cell2]] == 1:
                    cell2 = parent2[cell2]
                    children2[cell2, :] = 0
                    if cell2 == 0:
                        break

                # Reset intersection.
                candidates[cell1, cell2] = 1

    # from .tessellate import tessellate_data
    # from ...graphics import Viewer
    #
    # vert, tri = tessellate_data(surface1, ncells1[0], candidates1, children1,
    #                             acells1, position1, parent1, cell_params1,
    #                             gtol)
    # Viewer.add_triplot(vert, tri)
    #
    # vert, tri = tessellate_data(surface2, ncells2[0], candidates2, children2,
    #                             acells2, position2, parent2, cell_params2,
    #                             gtol)
    # Viewer.add_triplot(vert, tri)

    def _intersect(icell1, icell2):
        """
        Intersect cells using triangle-triangle intersection.
        """
        # Get parameters of tessellated cells.
        ntri1, triangles1 = tessellate_cell(icell1, children1, acells1,
                                            position1, parent1, cell_params1)
        ntri2, triangles2 = tessellate_cell(icell2, children2, acells2,
                                            position2, parent2, cell_params2)
        # Intersect triangles.
        for ii in range(0, ntri1):
            uv10, uv11, uv12 = triangles1[ii]
            ti1[0, :] = s1eval(*uv10, rtype='ndarray', domain='global')
            ti1[1, :] = s1eval(*uv11, rtype='ndarray', domain='global')
            ti1[2, :] = s1eval(*uv12, rtype='ndarray', domain='global')
            # Check for degenerate triangle.
            a1 = triangle_area(ti1)
            if a1 <= 1.0e-12:
                continue
            for jj in range(0, ntri2):
                uv20, uv21, uv22 = triangles2[jj]
                ti2[0, :] = s2eval(*uv20, rtype='ndarray', domain='global')
                ti2[1, :] = s2eval(*uv21, rtype='ndarray', domain='global')
                ti2[2, :] = s2eval(*uv22, rtype='ndarray', domain='global')
                # Check for degenerate triangle.
                a2 = triangle_area(ti2)
                if a2 <= 1.0e-12:
                    continue
                ni, pi = intersect_triangles(ti1, ti2, gtol)
                if ni == 2:
                    edges.append([ptotal[0], ptotal[0] + 1])
                    ptotal[0] += 2
                    etotal[0] += 1
                    points.append(pi[0])
                    points.append(pi[1])
                    # Get parameters.
                    # Surface 1.
                    ub, vb = barycentric_params(pi[0], ti1)
                    ur, vr = (1. - ub - vb) * uv10 + ub * uv11 + vb * uv12
                    params1.append([ur, vr])
                    points2d_s1.append([ur, vr, 0.])
                    ub, vb = barycentric_params(pi[1], ti1)
                    ur, vr = (1. - ub - vb) * uv10 + ub * uv11 + vb * uv12
                    params1.append([ur, vr])
                    points2d_s1.append([ur, vr, 0.])

                    # Surface 2.
                    ub, vb = barycentric_params(pi[0], ti2)
                    ur, vr = (1. - ub - vb) * uv20 + ub * uv21 + vb * uv22
                    params2.append([ur, vr])
                    points2d_s2.append([ur, vr, 0.])
                    ub, vb = barycentric_params(pi[1], ti2)
                    ur, vr = (1. - ub - vb) * uv20 + ub * uv21 + vb * uv22
                    params2.append([ur, vr])
                    points2d_s2.append([ur, vr, 0.])

    # Tessellate potential surfaces and intersect.
    s1eval = surface1.eval
    s2eval = surface2.eval
    etotal = [0]
    ptotal = [0]
    edges = []
    points = []
    points2d_s1 = []
    points2d_s2 = []
    params1 = []
    params2 = []
    ti1 = zeros((3, 3), dtype=float64)
    ti2 = zeros((3, 3), dtype=float64)
    for i in range(0, ncells1[0] + 1):
        for j in range(0, ncells2[0] + 1):
            if candidates[i, j] == 1:
                _intersect(i, j)

    # Check for no intersection.
    if ptotal[0] <= 1:
        return 0, [], [], [], [], []

    # Find the topology of the intersection curves.
    arr_points = array(points, dtype=float64)
    arr_params1 = array(params1, dtype=float64)
    arr_params2 = array(params2, dtype=float64)
    arr_points2d_s1 = array(points2d_s1, dtype=float64)
    arr_points2d_s2 = array(points2d_s2, dtype=float64)
    ncrvs, all_crv_ids = _trace_curves(ptotal[0], etotal[0], arr_points2d_s1,
                                       arr_points2d_s2, edges, ptol)

    # Refine intersection points.
    tol = gtol / 100.
    for i in range(ncrvs):
        for indx in all_crv_ids[i]:
            u1, v1 = arr_params1[indx]
            u2, v2 = arr_params2[indx]
            u1, v1, u2, v2, p = refine_ssi_point(surface1, surface2, u1,
                                                 v1, u2, v2, tol)
            arr_params1[indx, :] = [u1, v1]
            arr_params2[indx, :] = [u2, v2]
            arr_points[indx, :] = p

    # Filter out points.
    # no_filter = [False] * arr_points.shape[0]
    crv_size, crv_ids = _filter_points(ncrvs, all_crv_ids, arr_points, gtol)
    # crv_ids = all_crv_ids
    # crv_size = [len(row) for row in crv_ids]

    # Return
    return ncrvs, crv_size, crv_ids, arr_points, arr_params1, arr_params2


def _build_tess_data(cells):
    """
    Build arrays for tessellation.
    """
    n = len(cells)
    children = zeros((n, 5), dtype=int32)
    parent = zeros(n, dtype=int32)
    acells = zeros((n, 5), dtype=int32)
    position = zeros(n, dtype=int32)
    candidates = zeros(n, dtype=int32)
    cell_params = zeros((n, 4, 2), dtype=float64)

    for ci in cells:
        sid = ci.sid

        # Candidate
        if ci.is_cand:
            candidates[sid] = 1

        # Parameters
        cell_params[sid, 0, :] = [ci.u0, ci.v0]
        cell_params[sid, 1, :] = [ci.u1, ci.v0]
        cell_params[sid, 2, :] = [ci.u1, ci.v1]
        cell_params[sid, 3, :] = [ci.u0, ci.v1]

        # Parent
        if ci.parent is not None:
            parent[sid] = ci.parent.sid

        # Position
        if ci.position > 0:
            position[sid] = ci.position

        # Children
        if ci.sw is not None:
            children[sid, 0] = 1
            children[sid, 1] = ci.sw.sid
        if ci.nw is not None:
            children[sid, 0] = 1
            children[sid, 2] = ci.nw.sid
        if ci.se is not None:
            children[sid, 0] = 1
            children[sid, 3] = ci.se.sid
        if ci.ne is not None:
            children[sid, 0] = 1
            children[sid, 4] = ci.ne.sid

        # Adjacent cells
        if ci.n is not None:
            acells[sid, 0] = 4
            acells[sid, 4] = ci.n.sid
        if ci.e is not None:
            acells[sid, 0] = 4
            acells[sid, 3] = ci.e.sid
        if ci.s is not None:
            acells[sid, 0] = 4
            acells[sid, 1] = ci.s.sid
        if ci.w is not None:
            acells[sid, 0] = 4
            acells[sid, 2] = ci.w.sid

    return children, parent, acells, position, candidates, cell_params


def _trace_curves(ptotal, etotal, points2d_s1, points2d_s2, edges, tol):
    """
    Trace an unsorted collection of edges.
    """
    # Equivalence points2d_s1 tracking point use id's.
    points_use = zeros(ptotal, dtype=int32)
    vert_to_point = zeros(ptotal, dtype=int32)
    verts1 = zeros((ptotal, 3), dtype=float64)
    verts1[0] = points2d_s1[0]
    verts2 = zeros((ptotal, 3), dtype=float64)
    verts2[0] = points2d_s2[0]
    nverts = 1
    for i in range(1, ptotal):
        p1 = points2d_s1[i]
        p2 = points2d_s2[i]
        unique = True
        for j in range(0, nverts):
            v1 = verts1[j]
            v2 = verts2[j]
            if norm(v1 - p1) <= tol and norm(v2 - p2) <= tol:
                unique = False
                points_use[i] = j
                vert_to_point[j] = i
                break
        if unique:
            verts1[nverts] = points2d_s1[i]
            verts2[nverts] = points2d_s2[i]
            points_use[i] = nverts
            vert_to_point[nverts] = i
            nverts += 1

    # Build new edges
    new_edges = []
    point_count = zeros(nverts, dtype=int32)
    adj_edges = zeros((nverts, etotal), dtype=int32)
    visited = zeros(nverts, dtype=int32)
    point_to_point = zeros((nverts, nverts), dtype=int32)
    eid = 0
    for e in edges:
        pid1 = points_use[e[0]]
        pid2 = points_use[e[1]]
        # Remove duplicate edges if present.
        if (point_to_point[pid1, pid2] == 1 or
                point_to_point[pid2, pid1] == 1 or pid1 == pid2):
            continue
        point_to_point[pid1, pid2] = 1
        point_to_point[pid2, pid1] = 1
        new_edges.append([pid1, pid2])
        adj_edges[pid1, point_count[pid1]] = eid
        adj_edges[pid2, point_count[pid2]] = eid
        point_count[pid1] += 1
        point_count[pid2] += 1
        eid += 1

    # Process the curves until all points are visited.
    ncrvs = 0
    crv_ids = []
    process_curves = True
    while process_curves:
        # Try to find a point with only one adjacent edge. If all have more
        # than one adjacent edge it implies a closed curve. In that case start
        # anywhere.
        point_found = False
        is_closed = False
        pid1 = 0
        for i in range(nverts):
            # Try to find a single starting point.
            if point_count[i] == 1 and visited[i] == 0:
                pid1 = i
                point_found = True
                break
        # Select the first unvisited point if no single point was found.
        if not point_found:
            for i in range(nverts):
                if point_count[i] > 0 and visited[i] == 0:
                    pid1 = i
                    point_found = True
                    is_closed = True
                    break

        # Trace the topology of the curve using a DFS search.
        crv_id = []
        if point_found:
            # Non-recursive DFS.
            stack = [pid1]
            nstack = 1
            while nstack > 0:
                pid = stack.pop()
                nstack -= 1
                visited[pid] = 1
                crv_id.append(vert_to_point[pid])
                edge_list = adj_edges[pid, 0:point_count[pid]]
                for ei in edge_list:
                    for vi in new_edges[ei]:
                        if visited[vi] == 0:
                            stack.append(vi)
                            visited[vi] = 1
                            nstack += 1
            # Append first point if a closed curve was traced.
            if is_closed and visited[0] == 1:
                crv_id.append(crv_id[0])
            crv_ids.append(crv_id)
            ncrvs += 1
        else:
            process_curves = False

    return ncrvs, crv_ids


# def _filter_points(ncrvs, all_crv_ids, arr_points, gap_tol, dist_tol,
#                    no_filter):
#     """
#     Filter points based on flatness criteria and distance.
#     """
#     crv_ids = []
#     crv_size = []
#     for i in range(ncrvs):
#         ids = all_crv_ids[i]
#         nids = len(ids)
#         if nids == 2:
#             crv_ids.append([ids[0], ids[1]])
#             continue
#         i0 = 0
#         i1 = 2
#         fcrv = [ids[i0]]
#         while i1 < nids:
#             is_flat = False
#             dline = norm(arr_points[ids[i1]] - arr_points[ids[i0]])
#             for k in range(i0 + 1, i1):
#                 v0 = arr_points[ids[k]] - arr_points[ids[i0]]
#                 v1 = arr_points[ids[i1]] - arr_points[ids[k]]
#
#                 # Check for reversed points using the dot product and angle
#                 # between the vectors.
#                 dp = dot(v0, v1)
#                 if dp < 0.:
#                     angle = angle_between_vecs(v0, v1)
#                     if angle > 170.:
#                         is_flat = True
#                         break
#
#                 # Check for not filtered flag.
#                 if no_filter[ids[k]]:
#                     break
#
#                 # Check minimum distance.
#                 d = norm(v0)
#                 if d <= dist_tol:
#                     is_flat = True
#                     break
#
#                 # Check maximum gap.
#                 gap = norm(cross(v0, v1)) / dline
#                 if gap < gap_tol:
#                     is_flat = True
#                     break
#             if is_flat:
#                 i1 += 1
#             else:
#                 i0 = i1 - 1
#                 i1 = i0 + 2
#                 fcrv.append(ids[i0])
#         # Append last point or replace it if previous point is coincident.
#         if norm(arr_points[fcrv[-1]] - arr_points[ids[-1]]) > dist_tol:
#             fcrv.append(ids[-1])
#         else:
#             fcrv[-1] = ids[-1]
#         # Append curve length and id's.
#         crv_size.append(len(fcrv))
#         crv_ids.append(fcrv)
#     return crv_size, crv_ids


def _filter_points(ncrvs, all_crv_ids, arr_points, gtol):
    """
    Filter points based on flatness criteria and distance.
    """
    crv_ids = []
    crv_size = []
    for i in range(ncrvs):
        ids = all_crv_ids[i]
        nids = len(ids)
        if nids == 2:
            crv_ids.append([ids[0], ids[1]])
            continue
        i0 = 0
        i1 = 1
        i2 = 2
        fcrv = [ids[0]]
        while i2 < nids:
            is_flat = False
            v0 = arr_points[ids[i1]] - arr_points[ids[i0]]
            v1 = arr_points[ids[i2]] - arr_points[ids[i1]]

            # Check minimum distance.
            d = norm(v0)
            if d <= gtol:
                is_flat = True
            # Check for reversed points using the dot product and angle
            # between the vectors.
            if not is_flat:
                dp = dot(v0, v1)
                if dp < 0.:
                    angle = angle_between_vecs(v0, v1)
                    if angle > 170.:
                        is_flat = True
            # Adjust indices and/or add point to curve.
            if is_flat:
                i1 += 1
                i2 = i1 + 1
            else:
                fcrv.append(ids[i1])
                i0 = i1
                i1 = i0 + 1
                i2 = i1 + 1

        # Append last point or replace it if previous point is coincident.
        if norm(arr_points[fcrv[-1]] - arr_points[ids[-1]]) > gtol:
            fcrv.append(ids[-1])
        else:
            fcrv[-1] = ids[-1]
        # Append curve length and id's.
        crv_size.append(len(fcrv))
        crv_ids.append(fcrv)
    return crv_size, crv_ids


def refine_spi_point(surface, plane, u, v, tol):
    """
    Refine surface-plane intersection point.

    :param surface: Intersected surface.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param plane: Intersection plane.
    :type plane: :class:`.Plane`
    :param float u: Initial parameter.
    :param float v: Initial parameter.
    :param float tol: Refinement tolerance.

    :return: Refined parameters and point as NumPy array (u, v, pnt).
    :rtype: tuple
    """
    # Global parameters.
    nq = plane.vn.ijk
    origin = plane.p0.xyz
    vx = plane.vu.vxyz
    vy = plane.vv.vxyz

    seval = surface.eval
    sderiv = surface.deriv
    umin, umax = surface.au, surface.bu
    vmin, vmax = surface.av, surface.bv
    get_mult = surface.get_mult
    p = surface.p
    q = surface.q
    k = 0
    # Initial values
    p0 = seval(u, v, rtype='ndarray', domain='global')
    vp = p0 - origin
    dd = dot(vp, nq)
    q0 = p0 - dd * nq
    u0, v0 = u, v
    while k < 100:
        # Point on surface.
        p0 = seval(u, v, rtype='ndarray', domain='global')
        # Project point to plane.
        vp = p0 - origin
        dd = dot(vp, nq)
        q0 = p0 - dd * nq
        if norm(p0 - q0) <= tol:
            break
        # Surface unit normal.
        su = sderiv(u, v, 1, 0, rtype='ndarray', domain='global')
        sv = sderiv(u, v, 0, 1, rtype='ndarray', domain='global')
        denom = norm(cross(su, sv))
        if denom <= 1.0e-12:
            break
        np = cross(su, sv) / denom
        # Intersection of all three planes.
        dp = dot(np, p0)
        dq = dot(nq, q0)
        denom = norm(cross(np, nq))
        if denom <= 1.0e-12:
            break
        nn = cross(np, nq) / denom
        pq0 = mean([p0, q0], axis=0)
        dn = dot(nn, pq0)
        xi = (dp * cross(nq, nn) + dq * cross(nn, np) +
              dn * cross(np, nq)) / (dot(cross(np, nq), nn))
        # New increments.
        dp0 = xi - p0
        ru = cross(su, np)
        rv = cross(sv, np)
        # Check to see if current parameter is on a isoparameter of
        # the surface. If it is and its multiplicity is equal to the
        # degree, constrain the refinement process along the
        # isoparameter direction.
        dpq = dot(nq, p0 - q0)
        if p <= get_mult(u, 'u'):
            # Adjust v only.
            du = 0.
            if dot(nq, sv) * dpq >= 0:
                dv = -abs(dot(ru, dp0) / dot(ru, sv))
            else:
                dv = abs(dot(ru, dp0) / dot(ru, sv))
        elif q <= get_mult(v, 'v'):
            dv = 0.
            # Adjust u only.
            if dot(nq, su) * dpq >= 0.:
                du = -abs(dot(rv, dp0) / dot(rv, su))
            else:
                du = abs(dot(rv, dp0) / dot(rv, su))
        else:
            du = dot(rv, dp0) / dot(rv, su)
            dv = dot(ru, dp0) / dot(ru, sv)
        u += du
        v += dv
        # Check parameters.
        if u < umin:
            u = umin
        elif u > umax:
            u = umax
        if v < vmin:
            v = vmin
        elif v > vmax:
            v = vmax
        k += 1
    if k >= 100 or norm(p0 - q0) > tol:
        # Attempt Nelder-Mead.
        if p <= get_mult(u0, 'u'):
            # Adjust v only.
            u, v = _refine_spi_nm(surface, plane, u0, v0, 'v', tol)
        elif q <= get_mult(v0, 'v'):
            # Adjust u only.
            u, v = _refine_spi_nm(surface, plane, u0, v0, 'u', tol)
        else:
            # Adjust both.
            u, v = _refine_spi_nm(surface, plane, u0, v0, 'uv', tol)
        # Check parameters.
        if u < umin:
            u = umin
        elif u > umax:
            u = umax
        if v < vmin:
            v = vmin
        elif v > vmax:
            v = vmax
        p0 = surface.eval(u, v, domain='global', rtype='ndarray')
        up, vp = invert_points_on_plane([p0], plane)[0]
        q0 = plane.eval(up, vp, domain='global', rtype='ndarray')
        d = norm(p0 - q0)
        if d > tol and Settings.warnings:
            print('WARNING: Distance in SPI refinement exceeds tolerance.',
                  'Distance=', d)

    # Invert the refined point on the plane to get parameters.
    up, vp = invert_points_on_plane([p0], plane)[0]
    pi = mean([p0, q0], axis=0)
    return u, v, up, vp, pi


def _refine_spi_nm(surface, plane, u, v, d, tol):
    """
    Refine using Nelder-Mead optimization.
    """

    def _obj(x):
        if d == 'u':
            p0 = seval(x[0], v, domain='global', rtype='ndarray')
        elif d == 'v':
            p0 = seval(u, x[0], domain='global', rtype='ndarray')
        else:
            p0 = seval(x[0], x[1], domain='global', rtype='ndarray')
        vp = p0 - origin
        dd = dot(vp, nq)
        q0 = p0 - dd * nq
        return norm(q0 - p0)

    seval = surface.eval
    origin = plane.p0.xyz
    nq = plane.vn.ijk

    if d == 'u':
        x0 = u
    elif d == 'v':
        x0 = v
    else:
        x0 = [u, v]
    sol = minimize(_obj, x0, method='Nelder-Mead', tol=tol,
                   options={'ftol': tol})
    if d == 'u':
        return sol.x[0], v
    elif d == 'v':
        return u, sol.x[0]
    else:
        return sol.x


def refine_ssi_point(s1, s2, u1, v1, u2, v2, tol):
    """
    Refine surface-surface intersection point.

    :param s1: Intersected surface.
    :type s1: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param s2: Ohter intersected surface.
    :type s2: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float u1: Starting parameter for surface 1.
    :param float v1: Starting parameter for surface 1.
    :param float u2: Starting parameter for surface 2.
    :param float v2: Starting parameter for surface 2.
    :param float tol: Refinement tolerance.

    :return: Refined parameters and point as NumPy array (u1, v1, u2, v2, pnt).
    :rtype: tuple
    """
    # Methods for loop.
    s1eval = s1.eval
    s2eval = s2.eval
    s1deriv = s1.deriv
    s2deriv = s2.deriv
    umin1, umax1, vmin1, vmax1 = s1.au, s1.bu, s1.av, s1.bv
    umin2, umax2, vmin2, vmax2 = s2.au, s2.bu, s2.av, s2.bv
    # Initial values.
    k = 0
    p0 = s1eval(u1, v1, rtype='ndarray', domain='global')
    q0 = s2eval(u2, v2, rtype='ndarray', domain='global')
    d0 = norm(p0 - q0)
    u10, v10, u20, v20 = u1, v1, u2, v2
    while k < 100:
        if d0 <= tol:
            break
        # Surface unit normals.
        su1 = s1deriv(u1, v1, 1, 0, rtype='ndarray', domain='global')
        sv1 = s1deriv(u1, v1, 0, 1, rtype='ndarray', domain='global')
        denom = norm(cross(su1, sv1))
        if denom <= 1.0e-12:
            break
        np = cross(su1, sv1) / denom
        su2 = s2deriv(u2, v2, 1, 0, rtype='ndarray', domain='global')
        sv2 = s2deriv(u2, v2, 0, 1, rtype='ndarray', domain='global')
        denom = norm(cross(su2, sv2))
        if denom <= 1.0e-12:
            break
        nq = cross(su2, sv2) / denom
        # Intersection of all three planes.
        dp = dot(np, p0)
        dq = dot(nq, q0)
        denom = norm(cross(np, nq))
        if denom <= 1.0e-12:
            break
        nn = cross(np, nq) / denom
        pq0 = mean([p0, q0], axis=0)
        dn = dot(nn, pq0)
        xi = (dp * cross(nq, nn) + dq * cross(nn, np) +
              dn * cross(np, nq)) / (dot(cross(np, nq), nn))
        # New increments for surface 1.
        dp0 = xi - p0
        ru1 = cross(su1, np)
        rv1 = cross(sv1, np)
        u1 += dot(rv1, dp0) / dot(rv1, su1)
        v1 += dot(ru1, dp0) / dot(ru1, sv1)
        # Check parameters.
        if u1 < umin1:
            u1 = umin1
        elif u1 > umax1:
            u1 = umax1
        if v1 < vmin1:
            v1 = vmin1
        elif v1 > vmax1:
            v1 = vmax1
        # New increments for surface 2.
        dq0 = xi - q0
        ru2 = cross(su2, nq)
        rv2 = cross(sv2, nq)
        u2 += dot(rv2, dq0) / dot(rv2, su2)
        v2 += dot(ru2, dq0) / dot(ru2, sv2)
        # Check parameters.
        if u2 < umin2:
            u2 = umin2
        elif u2 > umax2:
            u2 = umax2
        if v2 < vmin2:
            v2 = vmin2
        elif v2 > vmax2:
            v2 = vmax2
        # New location.
        p0 = s1eval(u1, v1, rtype='ndarray', domain='global')
        q0 = s2eval(u2, v2, rtype='ndarray', domain='global')
        d0 = norm(p0 - q0)
        k += 1
    if k >= 100 or d0 > tol:
        # Attempt Nelder-Mead.
        u1, v1, u2, v2 = _refine_ssi_nm(s1, s2, u10, v10, u20, v20, tol)
        # Check parameters.
        if u1 < umin1:
            u1 = umin1
        elif u1 > umax1:
            u1 = umax1
        if v1 < vmin1:
            v1 = vmin1
        elif v1 > vmax1:
            v1 = vmax1
        if u2 < umin2:
            u2 = umin2
        elif u2 > umax2:
            u2 = umax2
        if v2 < vmin2:
            v2 = vmin2
        elif v2 > vmax2:
            v2 = vmax2
        p0 = s1eval(u1, v1, rtype='ndarray', domain='global')
        q0 = s2eval(u2, v2, rtype='ndarray', domain='global')
        d0 = norm(p0 - q0)
        if d0 > tol and Settings.warnings:
            print('WARNING: Distance in SSI refinement exceeds tolerance.',
                  'Distance=', d0)

    pi = mean([p0, q0], axis=0)
    return u1, v1, u2, v2, pi


def _refine_ssi_nm(surface1, surface2, u1, v1, u2, v2, tol):
    """
    Refine using Nelder-Mead optimization.
    """

    def _obj(x):
        # factor = 1.
        # if x[0] < surface1.au or x[0] > surface1.bu:
        #     factor = 1000.
        # elif x[1] < surface1.av or x[1] > surface1.bv:
        #     factor = 1000.
        # elif x[2] < surface2.au or x[2] > surface2.bu:
        #     factor = 1000.
        # elif x[3] < surface2.av or x[3] > surface2.bv:
        #     factor = 1000.

        p0 = s1eval(x[0], x[1], domain='global', rtype='ndarray')
        q0 = s2eval(x[2], x[3], domain='global', rtype='ndarray')
        return norm(q0 - p0)

    s1eval = surface1.eval
    s2eval = surface2.eval
    x0 = array([u1, v1, u2, v2], dtype=float64)
    sol = minimize(_obj, x0, method='Nelder-Mead', tol=tol,
                   options={'ftol': tol})
    return sol.x
