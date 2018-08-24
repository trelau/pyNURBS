from __future__ import division

from numpy import zeros, int32, float64, array, ceil

from pynurbs.config import Settings
from pynurbs.geometry.methods.geom_utils import is_curve_flat

# Array for looking up adjacent cell based on cell position and edge.
_adj_cell_pos1 = zeros((5, 5), dtype=int32)
_adj_cell_pos1[1, :] = [0, 2, 3, 1, 1]
_adj_cell_pos1[2, :] = [0, 2, 4, 2, 1]
_adj_cell_pos1[3, :] = [0, 4, 3, 1, 3]
_adj_cell_pos1[4, :] = [0, 4, 4, 2, 3]

# Array for looking up two adjacent cell positions.
_adj_cell_pos2 = zeros((5, 2), dtype=int32)
_adj_cell_pos2[1, :] = [2, 4]
_adj_cell_pos2[2, :] = [4, 3]
_adj_cell_pos2[3, :] = [1, 2]
_adj_cell_pos2[4, :] = [3, 1]

# Array for looking up adjacent edge parameters.
_adj_edge_param = zeros((5, 2), dtype=int32)
_adj_edge_param[1, :] = [3, 2]
_adj_edge_param[2, :] = [2, 1]
_adj_edge_param[3, :] = [0, 3]
_adj_edge_param[4, :] = [1, 0]


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


def adaptive_curve_tessellate(curve, tol=None):
    """
    Adaptive curve tessellation.

    :param curve: Curve to tessellate.
    :param float tol: Tolerance for checking curve flatness.

    :return: Status of tessellation method, number of subdivisions,
        vertices, and edges (success, nsub, verts, edges).
    :rtype: tuple
    """
    if tol is None:
        tol = Settings.ftol

    # Use control points if linear.
    # if curve.p == 1:
    #     return curve.cp, curve.uk[1:-1]

    # Parameters for subdivision.
    nvert = [0]
    verts = []
    params = []

    # Step 1: Define methods for recursive subdivision.
    def _subdivide(ci):
        """
        Recursive subdivision.
        """
        # Check flatness.
        cpi = ci.cp
        if not is_curve_flat(ci.n, cpi, tol):
            # Split curve.
            if ci.p == 1:
                mid = int(ceil((ci.uk.size - 1) / 2))
                ui = ci.uk[mid]
            else:
                ui = 0.5 * (ci.a + ci.b)
            c1, c2 = ci.split(ui, domain='global')
            # Subdivide.
            map(_subdivide, [c1, c2])
        else:
            params.append(ci.a)
            verts.append(curve.eval(ci.a, rtype='ndarray', domain='global'))
            nvert[0] += 2

    # Step 2: Use recursive subdivision to generate flat curves.
    _subdivide(curve)

    # Add last point.
    params.append(curve.b)
    verts.append(curve.eval(curve.b, rtype='ndarray', domain='global'))
    nvert[0] += 2

    # Convert to array.
    verts = array(verts, dtype=float64)
    params = array(params, dtype=float64)

    # Return results.
    return verts, params


def adaptive_surface_tessellate(surface, tol=None):
    """
    Adaptive surface tessellation.

    :param surface: Surface to tessellate.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float tol: Tolerance for checking surface flatness.

    :return: Status of tessellation method, number of subdivisions, vertices,
        and triangles
        (success, nsub, verts, triangles).
    """
    # Global parameters.
    gtol = Settings.gtol
    ptol = Settings.ptol
    if tol is None:
        ftol = Settings.ftol
    else:
        ftol = tol

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
        return True

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


def tessellate_data(surface, ncells, candidates, children, acells, position,
                    parent, cell_params, tol):
    """
    Tessellate cell data. Primarily for debugging.
    """
    # Parameters for tessellation.
    ptotal = 0
    tri_total = 0
    points = []
    conn = []
    seval = surface.eval

    for i in range(0, ncells + 1):
        if candidates[i] != 1:
            continue
        # Get parameters of tessellated cell.
        ntri, triangles = tessellate_cell(i, children, acells, position,
                                          parent, cell_params)
        # Determine triangle points and store.
        for ii in range(0, ntri):
            uv0, uv1, uv2 = triangles[ii]
            t0 = seval(*uv0, rtype='ndarray', domain='global')
            t1 = seval(*uv1, rtype='ndarray', domain='global')
            t2 = seval(*uv2, rtype='ndarray', domain='global')
            map(points.append, [t0, t1, t2])
            conn.append([ptotal, ptotal + 1, ptotal + 2])
            ptotal += 3
            tri_total += 1

    # Step 4: Equivalence vertices and build vertex and connectivity array.
    points_arr = array(points, dtype=float64)
    points_use = zeros(ptotal, dtype=int32)
    verts = zeros((ptotal, 3), dtype=float64)
    verts[0] = points_arr[0]
    nverts = 1
    gtol2 = tol * tol
    for i in range(1, ptotal):
        p = points_arr[i]
        unique = True
        for j in range(0, nverts):
            v = verts[j]
            d2 = ((v[0] - p[0]) ** 2. +
                  (v[1] - p[1]) ** 2. +
                  (v[2] - p[2]) ** 2.)
            if d2 <= gtol2:
                unique = False
                points_use[i] = j
                break
        if unique:
            verts[nverts] = points_arr[i]
            points_use[i] = nverts
            nverts += 1

    verts = verts[:nverts]
    triangles = zeros((tri_total, 3), dtype=int32)
    for i in range(0, tri_total):
        row = conn[i]
        triangles[i, :] = [points_use[row[0]],
                           points_use[row[1]],
                           points_use[row[2]]]

    # Return vertices and triangle connectivity.
    return verts, triangles


def tessellate_cell(csn, children, acells, position, parent, cell_params):
    """
    Tessellate a cell.

    :param int csn: Cell number.
    :param ndarray children: Array specifying children of each cell.
    :param ndarray acells: Array specifying the adjacent cells of each cell.
    :param ndarray position: Array specifying the position of each cell.
    :param ndarray parent: Array specifying the parent of each cell.
    :param ndarray cell_params: Array specifying the corner parameters of each
        cell.

    :return: The number of triangles and the parameter values for each of the
        three triangle vertices (ntri, triangles). The *triangles* list will
        be *ntri* x 3, where each row is a triangle and contains a list of
        three tuples containing two parameter values on the original surface.
        These parameter values represent a vertex of the triangle. The order of
        the parameter values should result in a normal vector oriented the same
        as the original surface.
    :rtype: tuple

    Reference: Anderson, J., Khamayseh, A., and Jean, B. "Adaptive Resolution
    Refinement," Technical Report, Los Alamos National Laboratory.
    """
    # Determine number of interior points and parameters on each edge.
    nprms = []
    edge_prms = []
    edge_order = [1, 3, 4, 2]
    simple = False
    for i in edge_order:
        adj_cells = _find_neighbors(csn, i, children, acells, position, parent)
        prms = _find_edge_params(i, adj_cells, cell_params)
        nprms.append(len(prms))
        edge_prms.append(prms)
        if len(prms) > 1:
            simple = True

    # Use simple triangulation if any edge has more than one interior point.
    if simple:
        # Make a single list of parameters in counter-clockwise order.
        # Edge 1
        all_params = [cell_params[csn, 0]]
        for uv in edge_prms[0]:
            all_params.append(uv)
        all_params.append(cell_params[csn, 1])
        # Edge 3
        for uv in edge_prms[1]:
            all_params.append(uv)
        all_params.append(cell_params[csn, 2])
        # Edge 4
        for uv in edge_prms[2]:
            all_params.append(uv)
        all_params.append(cell_params[csn, 3])
        # Edge 2
        for uv in edge_prms[3]:
            all_params.append(uv)
        all_params.append(cell_params[csn, 0])
        # Middle parameter.
        uv0 = cell_params[csn, 0]
        uv1 = cell_params[csn, 2]
        uvc = 0.5 * (uv0 + uv1)
        # Generate triangles
        triangles = []
        for i in range(len(all_params) - 1):
            uv0 = all_params[i]
            uv1 = all_params[i + 1]
            triangles.append([uv0, uv1, uvc])
        return len(triangles), triangles

    # Use predefined triangles.
    # Determine triangulation case by the number of interior points on each
    # edge.
    triangles = []
    triapp = triangles.append
    cprms = cell_params[csn, :]
    eprms = [row[0] for row in edge_prms if len(row) > 0]
    case = nprms[0] + nprms[1] * 10 + nprms[2] * 100 + nprms[3] * 1000

    # Case 0
    if case == 0:
        triapp([cprms[0], cprms[1], cprms[2]])
        triapp([cprms[2], cprms[3], cprms[0]])
        return len(triangles), triangles

    # Case 1
    if case == 1:
        triapp([cprms[0], eprms[0], cprms[3]])
        triapp([eprms[0], cprms[1], cprms[2]])
        triapp([cprms[2], cprms[3], eprms[0]])
        return len(triangles), triangles

    # Case 2
    if case == 10:
        triapp([cprms[0], cprms[1], eprms[0]])
        triapp([eprms[0], cprms[2], cprms[3]])
        triapp([cprms[3], cprms[0], eprms[0]])
        return len(triangles), triangles

    # Case 3
    if case == 11:
        triapp([cprms[0], eprms[0], cprms[3]])
        triapp([eprms[0], cprms[1], eprms[1]])
        triapp([eprms[1], cprms[2], cprms[3]])
        triapp([cprms[3], eprms[0], eprms[1]])
        return len(triangles), triangles

    # Case 4
    if case == 100:
        triapp([cprms[0], cprms[1], eprms[0]])
        triapp([cprms[1], cprms[2], eprms[0]])
        triapp([eprms[0], cprms[3], cprms[0]])
        return len(triangles), triangles

    # Case 5
    if case == 101:
        triapp([cprms[0], eprms[0], cprms[3]])
        triapp([eprms[0], eprms[1], cprms[3]])
        triapp([eprms[0], cprms[1], eprms[1]])
        triapp([cprms[1], cprms[2], eprms[1]])
        return len(triangles), triangles

    # Case 6
    if case == 110:
        triapp([cprms[0], cprms[1], eprms[0]])
        triapp([eprms[0], cprms[2], eprms[1]])
        triapp([eprms[1], cprms[0], eprms[0]])
        triapp([eprms[1], cprms[3], cprms[0]])
        return len(triangles), triangles

    # Case 7
    if case == 111:
        triapp([cprms[0], eprms[0], eprms[2]])
        triapp([eprms[0], cprms[1], eprms[1]])
        triapp([eprms[1], cprms[2], eprms[2]])
        triapp([eprms[2], eprms[0], eprms[1]])
        triapp([eprms[2], cprms[3], cprms[0]])
        return len(triangles), triangles

    # Case 8
    if case == 1000:
        triapp([cprms[0], cprms[1], eprms[0]])
        triapp([cprms[1], cprms[2], eprms[0]])
        triapp([cprms[2], cprms[3], eprms[0]])
        return len(triangles), triangles

    # Case 9
    if case == 1001:
        triapp([cprms[0], eprms[0], eprms[1]])
        triapp([eprms[0], cprms[1], cprms[2]])
        triapp([cprms[2], eprms[1], eprms[0]])
        triapp([cprms[2], cprms[3], eprms[1]])
        return len(triangles), triangles

    # Case 10
    if case == 1010:
        triapp([cprms[0], cprms[1], eprms[0]])
        triapp([eprms[0], cprms[2], eprms[1]])
        triapp([cprms[2], cprms[3], eprms[1]])
        triapp([eprms[1], cprms[0], eprms[0]])
        return len(triangles), triangles

    # Case 11
    if case == 1011:
        triapp([cprms[0], eprms[0], eprms[2]])
        triapp([eprms[0], eprms[1], eprms[2]])
        triapp([eprms[0], cprms[1], eprms[1]])
        triapp([eprms[1], cprms[2], eprms[2]])
        triapp([cprms[2], cprms[3], eprms[2]])
        return len(triangles), triangles

    # Case 12
    if case == 1100:
        triapp([cprms[0], cprms[1], eprms[1]])
        triapp([cprms[1], cprms[2], eprms[0]])
        triapp([eprms[0], eprms[1], cprms[1]])
        triapp([eprms[0], cprms[3], eprms[1]])
        return len(triangles), triangles

    # Case 13
    if case == 1101:
        triapp([cprms[0], eprms[0], eprms[2]])
        triapp([eprms[0], cprms[1], cprms[2]])
        triapp([cprms[2], eprms[1], eprms[0]])
        triapp([eprms[1], eprms[2], eprms[0]])
        triapp([eprms[1], cprms[3], eprms[2]])
        return len(triangles), triangles

    # Case 14
    if case == 1110:
        triapp([cprms[0], cprms[1], eprms[0]])
        triapp([eprms[0], eprms[2], cprms[0]])
        triapp([eprms[0], cprms[2], eprms[1]])
        triapp([eprms[1], eprms[2], eprms[0]])
        triapp([eprms[1], cprms[3], eprms[2]])
        return len(triangles), triangles

    # Case 15
    if case == 1111:
        triapp([cprms[0], eprms[0], eprms[3]])
        triapp([eprms[0], eprms[1], eprms[3]])
        triapp([eprms[0], cprms[1], eprms[1]])
        triapp([eprms[1], cprms[2], eprms[2]])
        triapp([eprms[2], eprms[3], eprms[1]])
        triapp([eprms[2], cprms[3], eprms[3]])
        return len(triangles), triangles
    # Return empty list.
    return 0, []


def _find_neighbors(csn, esn, children, acells, position, parent):
    """
    Find adjacent cells of a given cell.
    """
    # Find adjacent cell.
    acell = acells[csn, esn]

    # If adjacent cell has no children or is at a subdivision level equal to
    # or higher than the current cell, there are no interior points.
    if children[acell, 0] == 0 or acell == 0:
        return []

    # Traverse up the tree from the current cell until you find a common
    # parent. Track the parent cell numbers along the way.
    parents = [csn]
    pcell = csn
    adj_pcell = parent[acell]
    while parent[pcell] != adj_pcell:
        parents.append(parent[pcell])
        pcell = parent[pcell]

    # Traverse back down the tree from the highest parent cell to the
    # current cell, tracking adjacent cells along the way. If a cell is reached
    # before you get back to current cell, then the current cell has no
    # interior points.
    parents.reverse()
    for icell in parents[1:]:
        pos = position[icell]
        acell = children[acell, _adj_cell_pos1[pos, esn]]
        if acell == 0:
            return []

    # For each subdivision level beyond current cell, gather the two
    # adjacent cells next to specified edge. These will have interior points
    # on the current edge.
    adj_cells = []
    stack = [acell]
    nstack = 1
    pos1 = _adj_cell_pos2[esn, 0]
    pos2 = _adj_cell_pos2[esn, 1]
    while nstack > 0:
        acell = stack.pop()
        nstack -= 1
        if acell == 0 or children[acell, 0] == 0:
            adj_cells.append(acell)
        else:
            stack.append(children[acell, pos2])
            stack.append(children[acell, pos1])
            nstack += 2
    return adj_cells


def _find_edge_params(esn, adj_cells, cell_params):
    """
    Find interior points along an edge in a cell.
    """
    nadj = len(adj_cells)
    if nadj == 1:
        return []
    params = []
    for i in range(0, nadj):
        icell = adj_cells[i]
        if icell == 0:
            continue
        if i == 0:
            # Get single corner.
            indx = _adj_edge_param[esn, 1]
            uv = cell_params[icell, indx]
            params.append(uv)
        elif 0 < i < nadj - 1:
            if adj_cells[i - 1] == 0:
                # Get both corners if previous cell was empty.
                indx = _adj_edge_param[esn, 0]
                uv = cell_params[icell, indx]
                params.append(uv)
                indx = _adj_edge_param[esn, 1]
                uv = cell_params[icell, indx]
                params.append(uv)
            else:
                # Get singe corner.
                indx = _adj_edge_param[esn, 1]
                uv = cell_params[icell, indx]
                params.append(uv)
        elif i == nadj - 1:
            # Get parameter if previous cell was empty.
            if adj_cells[i - 1] == 0:
                # Get singe corner.
                indx = _adj_edge_param[esn, 0]
                uv = cell_params[icell, indx]
                params.append(uv)
    return params
