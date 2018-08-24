from __future__ import division

from numpy import array, cross, dot, float64
from numpy.linalg import norm
from scipy.optimize import minimize

from pynurbs.config import Settings
from pynurbs.geometry.methods.geom_utils import (invert_point_in_triangle,
                                                 is_curve_flat,
                                                 is_surface_flat)
from pynurbs.geometry.methods.intersect_bbox import (bbox_intersects_plane,
                                                     bbox_ray_intersect,
                                                     bboxes_intersect)
from pynurbs.geometry.methods.project import project_point_to_line


def intersect_curve_line(curve, line, itol=None, t0=None, t1=None):
    """
    Find the intersection points of a line and a curve.

    :param line:
    :param curve:
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
    vline = line.v.ijk

    # Step 1: Check to see if endpoints of lines are within itol of line.
    pc0 = curve.eval(0.)
    pc1 = curve.eval(1.)

    u0_list = []
    _, _, proj = project_point_to_line(pc0, line)
    if proj[0][2] <= itol:
        u0_list.append((curve.a, proj[0][1]))

    _, _, proj = project_point_to_line(pc1, line)
    if proj[0][2] <= itol:
        u0_list.append((curve.b, proj[0][1]))

    # Step 2: Decompose curve into Bezier segments.
    _, bezier_curves = curve.decompose()

    # Step 3: Define methods for recursive subdivision.
    def _subdivide(ci):
        """
        Recursive subdivision to find possible intersections.
        """
        nsub[0] += 1
        cpi = ci.cp
        if _is_flat(ci.n, cpi):
            _intersect(ci, cpi)
        else:
            c1, c2 = ci.split(0.5)
            if _candidate_intersect(c1):
                _subdivide(c1)
            if _candidate_intersect(c2):
                _subdivide(c2)

    def _is_flat(ni, cpi):
        """
        Check curve flatness.
        """
        return is_curve_flat(ni, cpi, ftol)

    def _candidate_intersect(ci):
        """
        Check for potential intersection.
        """
        return bbox_ray_intersect(ci.get_bbox(), p0, vline, t0, t1, itol)

    def _intersect(ci, cpi):
        """
        Intersect lines.
        """
        pi1, pi2 = cpi[0], cpi[-1]
        d4321 = dot(vline, pi2 - pi1)
        d1321 = dot(pi1 - p0, pi2 - pi1)
        d4343 = dot(vline, vline)
        d2121 = dot(pi2 - pi1, pi2 - pi1)
        denom = d2121 * d4343 - d4321 * d4321
        if abs(denom) <= 1.0e-12:
            # If parallel add a point halfway between as starting point.
            u12 = ci.local_to_global_param(0.5)
            u0_list.append((u12, line.v.mag / 2.))
        else:
            d1343 = dot(pi1 - p0, vline)
            numer = d1343 * d4321 - d1321 * d4343
            u12 = numer / denom
            u34 = (d1343 + u12 * d4321) / d4343
            u12 = ci.local_to_global_param(u12)
            u0_list.append((u12, u34))

    # Step 4: Perform recursive subdivision.
    nsub = [0]
    for cbez in bezier_curves:
        if _candidate_intersect(cbez):
            _subdivide(cbez)

    # Step 5: Refine initial points.
    def _obj(x):
        factor = 1.
        # Penalize objective function if variable is outside domain.
        if x[0] < a or x[0] > b:
            factor = 1000.
        pic = ceval(x[0], rtype='ndarray', domain='global')
        pil = p0 + x[1] * vline
        return norm(pic - pil) * factor

    ceval = curve.eval
    check_param = curve.check_param
    results = []
    a, b = curve.a, curve.b
    tol = itol / 100.
    for u01, u02 in u0_list:
        u012 = array([u01, u02], dtype=float64)
        sol = minimize(_obj, u012, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        u1, u2 = sol.x
        u1 = check_param(u1)
        p1 = ceval(u1, rtype='ndarray', domain='global')
        p2 = p0 + u2 * vline
        if norm(p1 - p2) > itol:
            continue
        unique = True
        for _, p12 in results:
            if norm(p12 - p1) <= itol:
                unique = False
                break
        if not unique:
            continue
        results.append([(u1, u2), p1])

    npts = len(results)
    return nsub[0], npts, results


def intersect_icurve_line(icrv, line, itol=None, t0=None, t1=None):
    """
    Find the intersection points between an intersection curve and a line.

    :param icrv:
    :param line:
    :param itol:
    :param t0:
    :param t1:

    :return:
    """
    # Global parameters.
    if itol is None:
        itol = Settings.gtol

    # Get line origin and direction vector.
    p0 = line.p0.xyz
    vline = line.v.ijk

    # Intersect 3-D for initial parameters.
    crv3d = icrv.crv3d
    itol_ = max(icrv.ftol, itol)
    nsub, npts, r3d = intersect_curve_line(crv3d, line, itol_, t0, t1)
    if npts == 0:
        return nsub, npts, r3d

    # Define objective function for point refinement.
    def _obj(x):
        # factor = 1.
        # Penalize objective function if variable is outside domain.
        # if x[0] < a or x[0] > b:
        #     factor = 1000.
        pic = ceval(x[0], rtype='ndarray', domain='global')
        pil = p0 + x[1] * vline
        return norm(pic - pil)

    ceval = icrv.eval
    check_param = icrv.check_param
    results = []
    # a, b = icrv.a, icrv.b
    tol = itol / 100.
    for u01, u02 in [r[0] for r in r3d]:
        u012 = array([u01, u02], dtype=float64)
        sol = minimize(_obj, u012, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        u1, u2 = sol.x
        u1 = check_param(u1)
        p1 = ceval(u1, rtype='ndarray', domain='global')
        p2 = p0 + u2 * vline
        if norm(p1 - p2) > itol:
            continue
        unique = True
        for _, p12 in results:
            if norm(p12 - p1) <= itol:
                unique = False
                break
        if not unique:
            continue
        # Put line then curve parameter.
        results.append([(u1, u2), 0.5 * (p1 + p2)])

    npts = len(results)
    return nsub, npts, results


def intersect_curve_curve(curve1, curve2, itol=None):
    """
    Find the intersection points of two curves.

    :param curve1: Curve 1 to intersect.
    :type curve1: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param curve2: Curve 2 to intersect.
    :type curve2: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param float itol: Intersection tolerance.

    :return: Number of subdivisions, intersection points, and results as a
        list of lists nsub, npts, [[(u1, u2), pi], ...].
    :rtype: tuple
    """
    # Global parameters
    ftol = Settings.ftol
    if itol is None:
        itol = Settings.gtol

    # Step 1: Check bounding boxes for possible intersection.
    bb1 = curve1.get_bbox()
    bb2 = curve2.get_bbox()
    if not bboxes_intersect(bb1, bb2, itol):
        return 0, 0, []

    # Step 2: Decompose each curve into Bezier segments.
    _, bezier_curves1 = curve1.decompose()
    _, bezier_curves2 = curve2.decompose()

    # Step 3: Build initial intersection list.
    cci_list = []
    for c1 in bezier_curves1:
        bb1 = c1.get_bbox()
        for c2 in bezier_curves2:
            bb2 = c2.get_bbox()
            if bboxes_intersect(bb1, bb2, itol):
                cci_list.append((c1, c2))

    # Step 4: Define methods for recursive subdivision.
    def _subdivide(ci1, ci2):
        """
        Recursive subdivision to find possible intersections.
        """
        subdivisions[0] += 1
        cpi1 = ci1.cp
        cpi2 = ci2.cp
        is_flat1 = _is_flat(ci1.n, cpi1)
        is_flat2 = _is_flat(ci2.n, cpi2)
        if is_flat1 and is_flat2:
            _intersect(ci1, ci2, cpi1, cpi2)
        else:
            cs12 = ci1.split(0.5)
            cs34 = ci2.split(0.5)
            for csi1 in cs12:
                for csi2 in cs34:
                    if _candidate_intersect(csi1, csi2):
                        _subdivide(csi1, csi2)

    def _is_flat(ni, cpi):
        """
        Check curve flatness.
        """
        return is_curve_flat(ni, cpi, ftol)

    def _intersect(ci1, ci2, cpi1, cpi2):
        """
        Intersect lines.
        """
        pi1, pi2 = cpi1[0], cpi1[-1]
        pi3, pi4 = cpi2[0], cpi2[-1]
        d4321 = dot(pi4 - pi3, pi2 - pi1)
        d1321 = dot(pi1 - pi3, pi2 - pi1)
        d4343 = dot(pi4 - pi3, pi4 - pi3)
        d2121 = dot(pi2 - pi1, pi2 - pi1)
        denom = d2121 * d4343 - d4321 * d4321
        if abs(denom) <= 1.0e-12:
            # If parallel add a point halfway between as starting point.
            ui1 = ci1.local_to_global_param(0.5)
            ui2 = ci2.local_to_global_param(0.5)
            u0_list.append((ui1, ui2))
        else:
            d1343 = dot(pi1 - pi3, pi4 - pi3)
            numer = d1343 * d4321 - d1321 * d4343
            u12 = numer / denom
            u34 = (d1343 + u12 * d4321) / d4343
            ui1 = ci1.local_to_global_param(u12)
            ui2 = ci2.local_to_global_param(u34)
            u0_list.append((ui1, ui2))

    def _candidate_intersect(ci1, ci2):
        """
        Check if the two curves intersect.
        """
        bbi1 = ci1.get_bbox()
        bbi2 = ci2.get_bbox()
        return bboxes_intersect(bbi1, bbi2, itol)

    # Step 5: Perform recursive subdivision to find initial points.
    u0_list = []
    subdivisions = [0]
    for c1, c2 in cci_list:
        _subdivide(c1, c2)

    # Refine initial points.
    def _obj(x):
        factor = 1.
        # Penalize objective function if variable is outside domain.
        if x[0] < a1 or x[0] > b1:
            factor = 1000.
        if x[1] < a2 or x[1] > b2:
            factor = 1000.
        pi1 = c1eval(x[0], rtype='ndarray', domain='global')
        pi2 = c2eval(x[1], rtype='ndarray', domain='global')
        return norm(pi2 - pi1) * factor

    def _add_pnt(u1i, u2i, pi):
        """
        Add point to intersection results.
        """
        unique = True
        for _, pii in results:
            if norm(pii - pi) <= itol:
                unique = False
                break
        if not unique:
            return False
        results.append([(u1i, u2i), pi])
        return True

    # Manually check curve endpoint intersection to improve robustness.
    results = []
    p10 = curve1.eval(0., rtype='ndarray')
    p11 = curve1.eval(1., rtype='ndarray')
    p20 = curve2.eval(0., rtype='ndarray')
    p21 = curve2.eval(1., rtype='ndarray')
    if norm(p10 - p20) <= itol:
        _add_pnt(curve1.a, curve2.a, 0.5 * (p10 + p20))
    elif norm(p10 - p21) <= itol:
        _add_pnt(curve1.a, curve2.b, 0.5 * (p10 + p21))
    if norm(p11 - p20) <= itol:
        _add_pnt(curve1.b, curve2.a, 0.5 * (p11 + p20))
    elif norm(p11 - p21) <= itol:
        _add_pnt(curve1.b, curve2.b, 0.5 * (p11 + p21))

    c1eval = curve1.eval
    c2eval = curve2.eval
    check_param1 = curve1.check_param
    check_param2 = curve2.check_param
    a1, b1 = curve1.a, curve1.b
    a2, b2 = curve2.a, curve2.b
    tol = itol / 100.
    for u01, u02 in reversed(u0_list):
        u012 = array([u01, u02], dtype=float64)
        sol = minimize(_obj, u012, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        u1, u2 = sol.x
        u1 = check_param1(u1)
        u2 = check_param2(u2)
        p1 = c1eval(u1, rtype='ndarray', domain='global')
        p2 = c2eval(u2, rtype='ndarray', domain='global')
        if norm(p1 - p2) > itol:
            continue
        _add_pnt(u1, u2, 0.5 * (p1 + p2))

    npts = len(results)
    nsub = subdivisions[0]
    return nsub, npts, results


def intersect_curve_icurve(curve, icrv, itol=None):
    """
    Find the intersection points between a 3-D curve and intersection curve.

    :param curve:
    :param icrv:
    :param itol:

    :return:
    """
    # Global parameters
    if itol is None:
        itol = Settings.gtol

    # Intersect 3-D curves.
    crv3d = icrv.crv3d
    itol_ = max(icrv.ftol, itol)
    nsub, npts, r3d = intersect_curve_curve(curve, crv3d, itol_)
    if npts == 0:
        return nsub, npts, r3d

    # Refine initial points.
    def _obj(x):
        # factor = 1.
        # Penalize objective function if variable is outside domain.
        # if x[0] < a1 or x[0] > b1:
        #     factor = 1000.
        # if x[1] < a2 or x[1] > b2:
        #     factor = 1000.
        pi1 = c1eval(x[0], rtype='ndarray', domain='global')
        pi2 = c2eval(x[1], rtype='ndarray', domain='global')
        return norm(pi2 - pi1)

    c1eval = curve.eval
    c2eval = icrv.eval
    check_param1 = curve.check_param
    check_param2 = icrv.check_param
    results = []
    # a1, b1 = curve.a, curve.b
    # a2, b2 = icrv.a, icrv.b
    tol = itol / 100.
    for u01, u02 in [r[0] for r in r3d]:
        u012 = array([u01, u02], dtype=float64)
        sol = minimize(_obj, u012, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        u1, u2 = sol.x
        u1 = check_param1(u1)
        u2 = check_param2(u2)
        p1 = c1eval(u1, rtype='ndarray', domain='global')
        p2 = c2eval(u2, rtype='ndarray', domain='global')
        if norm(p1 - p2) > itol:
            continue
        pi = array(0.5 * (p1 + p2), dtype=float64)
        unique = True
        for _, pii in results:
            if norm(pii - pi) <= itol:
                unique = False
                break
        if unique:
            results.append([(u1, u2), pi])

    npts = len(results)
    return nsub, npts, results


def intersect_icurve_icurve(icrv1, icrv2, itol=None):
    """
    Find the intersection points of two intersection curves.

    :param icrv1:
    :param icrv2:
    :param itol:

    :return:
    """
    # Global parameters
    if itol is None:
        itol = Settings.gtol

    # Intersect 3-D curves.
    crv3d1, crv3d2 = icrv1.crv3d, icrv2.crv3d
    itol_ = max(icrv1.ftol + icrv2.ftol, itol)
    nsub, npts, r3d = intersect_curve_curve(crv3d1, crv3d2, itol_)
    if npts == 0:
        return nsub, npts, r3d

    # Objective function for refinment
    def _obj(x):
        # factor = 1.
        # Penalize objective function if variable is outside domain.
        # if x[0] < a1 or x[0] > b1:
        #     print 'x1', x, a1, b1
        #     factor = 1000.
        # if x[1] < a2 or x[1] > b2:
        #     print 'x2', x, a1, b2
        #     factor = 1000.
        pi1 = c1eval(x[0], rtype='ndarray', domain='global')
        pi2 = c2eval(x[1], rtype='ndarray', domain='global')
        return norm(pi2 - pi1)

    c1eval = icrv1.eval
    c2eval = icrv2.eval
    check_param1 = icrv1.check_param
    check_param2 = icrv2.check_param
    results = []
    # a1, b1 = icrv1.a, icrv1.b
    # a2, b2 = icrv2.a, icrv2.b
    tol = itol / 100.
    for u01, u02 in [r[0] for r in r3d]:
        u012 = array([u01, u02], dtype=float64)
        sol = minimize(_obj, u012, method='Nelder-Mead', tol=tol)
        u1, u2 = sol.x
        u1 = check_param1(u1)
        u2 = check_param2(u2)
        p1 = c1eval(u1, rtype='ndarray', domain='global')
        p2 = c2eval(u2, rtype='ndarray', domain='global')
        if norm(p1 - p2) > itol:
            continue
        pi = array(0.5 * (p1 + p2), dtype=float64)
        unique = True
        for _, pii in results:
            if norm(pii - pi) <= itol:
                unique = False
                break
        if unique:
            results.append([(u1, u2), pi])

    npts = len(results)
    return nsub, npts, results


def intersect_curve_plane(curve, plane, itol=None):
    """
    Find the intersection points of a curve and a plane.

    :param curve: Curve to intersect.
    :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param Plane plane: Intersection plane.
    :param float itol: Intersection tolerance.

    :return: Number of subdivisions, intersections, and list containing
        intersection results nsub, npts, [[ui, pi], ...].
    :rtype: list
    """
    # Global parameters.
    p0 = plane.p0.xyz
    pnorm = plane.vn.ijk
    ftol = Settings.ftol
    if itol is None:
        itol = Settings.gtol

    # Step 1: Check bounding box and plane for possible intersection.
    if not bbox_intersects_plane(curve.get_bbox(), plane, itol):
        return 0, 0, []

    # Step 2: Decompose curve into Bezier segments.
    _, bezier_curves = curve.decompose()

    # Step 3: Define methods to subdivide, intersect, and check potential
    # curve-plane intersections.
    def _subdivide(ci):
        """
        Recursive subdivision for potential intersection points.
        """
        subdivisions[0] += 1
        cpi = ci.cp
        if _is_flat(ci.n, cpi):
            _intersect(ci, cpi)
        else:
            c1, c2 = ci.split(0.5)
            if _candidate_intersect(c1):
                _subdivide(c1)
            if _candidate_intersect(c2):
                _subdivide(c2)

    def _is_flat(ni, cpi):
        """
        Check the flatness of the curve.
        """
        return is_curve_flat(ni, cpi, ftol)

    def _intersect(ci, cpi):
        """
        Approximate intersection point using line-plane intersection.
        """
        pi0, pi1 = cpi[0], cpi[-1]
        denom = dot(pnorm, pi1 - pi0)
        if abs(denom) <= 1.0e-12:
            # If parallel use a starting point at midpoints.
            ui = ci.local_to_global_param(0.5)
            u0_list.append(ui)
        else:
            ui = dot(pnorm, p0 - pi0) / denom
            ui = ci.local_to_global_param(ui)
            u0_list.append(ui)

    def _candidate_intersect(ci):
        """
        Check for possible intersection using bounding box and plane.
        """
        bbox = ci.get_bbox()
        return bbox_intersects_plane(bbox, plane, itol)

    # Step 4: Discard initial Bezier curves that do not contain potential
    # intersection points.
    candidate_curves = []
    while len(bezier_curves) > 0:
        c = bezier_curves.pop(0)
        if _candidate_intersect(c):
            candidate_curves.append(c)

    # Step 5: Use recursive subdivision to find potential intersection points
    # and starting parameters.
    u0_list = []
    subdivisions = [0]
    for c in candidate_curves:
        _subdivide(c)

    # Step 6: Use solver to refine potential intersection points and
    # parameters.
    def _obj(x):
        factor = 1.
        # Penalize objective function if variable is outside domain.
        if x[0] < a or x[0] > b:
            factor = 1000.
        pc = ceval(x[0], rtype='ndarray', domain='global')
        pp = pc - pnorm * dot(pnorm, pc - p0)
        return norm(pc - pp) * factor

    results = []
    ceval = curve.eval
    check_param = curve.check_param
    a, b = curve.a, curve.b
    tol = itol / 100.
    for u0 in u0_list:
        sol = minimize(_obj, u0, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        u = sol.x[0]
        u = check_param(u)
        p = ceval(u, rtype='ndarray', domain='global')
        # Project point to plane and use midpoint.
        p_pln = p - pnorm * dot(pnorm, p - p0)
        p = 0.5 * (p_pln + p)
        unique = True
        for _, pi in results:
            if norm(p - pi) <= itol:
                unique = False
                break
        if not unique:
            continue
        results.append([u, p])

    npts = len(results)
    nsub = subdivisions[0]
    return nsub, npts, results


def intersect_icurve_plane(icrv, plane, itol=None):
    """
    Find the intersection points of an intersection curve and a plane.

    :param icrv:
    :param plane:
    :param itol:
    :return:
    """
    # Global parameters.
    p0 = plane.p0.xyz
    pnorm = plane.vn.ijk
    if itol is None:
        itol = Settings.gtol

    # Intersect 3-D curve.
    itol_ = max(icrv.ftol, itol)
    nsub, npts, r3d = intersect_curve_plane(icrv.crv3d, plane, itol_)
    if npts == 0:
        return nsub, npts, r3d

    # Objective function for refinement.
    def _obj(x):
        # factor = 1.
        # Penalize objective function if variable is outside domain.
        # if x[0] < a or x[0] > b:
        #     factor = 1000.
        pc = ceval(x[0], rtype='ndarray', domain='global')
        pp = pc - pnorm * dot(pnorm, pc - p0)
        return norm(pc - pp)

    ceval = icrv.eval
    check_param = icrv.check_param
    results = []
    # a, b = icrv.a, icrv.b
    tol = itol / 100.
    for u0 in [r[0] for r in r3d]:
        sol = minimize(_obj, u0, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        u = sol.x[0]
        u = check_param(u)
        p = ceval(u, rtype='ndarray', domain='global')
        # Project point to plane and use midpoint.
        p_pln = p - pnorm * dot(pnorm, p - p0)
        p = 0.5 * (p_pln + p)
        unique = True
        for _, pi in results:
            if norm(p - pi) <= itol:
                unique = False
                break
        if not unique:
            continue
        results.append([u, p])

    npts = len(results)
    return nsub, npts, results


def intersect_curve_surface(curve, surface, itol=None):
    """
    Find the intersection points of a curve and a surface.

    :param curve: Curve to intersect.
    :type curve: :class:`.BezierCurve` or :class:`.NurbsCurve`
    :param surface: Surface to intersect.
    :type surface: :class:`.BezierSurface` or :class:`.NurbsSurface`
    :param float itol: Intersection tolerance.

    :return: Number of subdivisions, intersections, and list containing
        intersection results nsub, npts, [[(cui, su, sv), pi], ...].
    :rtype: list
    """
    # Global parameters
    ftol = Settings.ftol
    if itol is None:
        itol = Settings.gtol

    # Step 1: Check bounding boxes for possible intersection.
    bb1 = curve.get_bbox()
    bb2 = surface.get_bbox()
    if not bboxes_intersect(bb1, bb2, itol):
        return 0, 0, []

    # Step 2: Decompose each into Bezier segments.
    _, bezier_curves = curve.decompose()
    bezier_patches = []
    _, _, surfs = surface.decompose()
    for row in surfs:
        for s in row:
            bezier_patches.append(s)

    # Step 3: Build initial intersection list.
    csi_list = []
    for c in bezier_curves:
        bb1 = c.get_bbox()
        for s in bezier_patches:
            bb2 = s.get_bbox()
            if bboxes_intersect(bb1, bb2, itol):
                csi_list.append((c, s))

    # Step 4: Define methods for recursive subdivision.
    def _subdivide(ci, si):
        """
        Recursive subdivision to find possible intersections.
        """
        nsub[0] += 1
        cpi = ci.cp
        spi = si.cp
        is_flat1 = _is_curve_flat(ci.n, cpi)
        is_flat2 = _is_surface_flat(si.n, si.m, spi)
        if is_flat1 and is_flat2:
            _intersect(ci, cpi, si, spi)
        else:
            ci12 = ci.split(0.5)
            si1234 = si.split(0.5, 0.5)
            for cii in ci12:
                for sii in si1234:
                    if _candidate_intersect(cii, sii):
                        _subdivide(cii, sii)

    def _is_curve_flat(ni, cpi):
        """
        Check curve flatness.
        """
        return is_curve_flat(ni, cpi, ftol)

    def _is_surface_flat(ni, mi, cpi):
        """
        Check surface flatness.
        """
        return is_surface_flat(ni, mi, cpi, ftol)

    def _intersect(ci, cpi, si, spi):
        """
        Intersect line and plane.
        """
        p0 = spi[0, 0]
        vu, vv = spi[-1, 0] - p0, spi[0, -1] - p0
        p1, p2 = cpi[0], cpi[-1]
        pnorm = cross(vu, vv)
        pnorm /= norm(pnorm)
        denom = dot(pnorm, p2 - p1)
        if abs(denom) <= 1.0e-12:
            # If parallel use a starting point at midpoints.
            ti = ci.local_to_global_param(0.5)
            ui = si.local_to_global_param('u', 0.5)
            vi = si.local_to_global_param('v', 0.5)
            tuv0_list.append((ti, ui, vi))
        else:
            ti = dot(pnorm, p0 - p1) / denom
            ti = ci.local_to_global_param(ti)
            pi = p1 + ti * (p2 - p1)
            tri3d = array([spi[0, 0], spi[-1, 0], spi[0, -1]], dtype=float64)
            tri2d = array([[si.au, si.av],
                           [si.bu, si.av],
                           [si.au, si.bv]], dtype=float64)
            ui, vi = invert_point_in_triangle(pi, tri3d, tri2d, False)
            tuv0_list.append((ti, ui, vi))

    def _candidate_intersect(ci, si):
        """
        Check if the curve and surface intersect.
        """
        bbi1 = ci.get_bbox()
        bbi2 = si.get_bbox()
        return bboxes_intersect(bbi1, bbi2, itol)

    # Step 5: Perform recursive subdivision to find initial points.
    tuv0_list = []
    nsub = [0]
    for c, s in csi_list:
        _subdivide(c, s)

    # Step 6: Refine points.
    def _obj(x):
        factor = 1.
        # Penalize objective function if variable is outside domain.
        if x[0] < a or x[0] > b:
            factor = 1000.
        if x[1] < au or x[1] > bu:
            factor = 1000.
        if x[2] < av or x[2] > bv:
            factor = 1000.
        pci = ceval(x[0], rtype='ndarray', domain='global')
        psi = seval(x[1], x[2], rtype='ndarray', domain='global')
        return norm(pci - psi) * factor

    ceval = curve.eval
    seval = surface.eval
    results = []
    a, b = curve.a, curve.b
    au, bu, av, bv = surface.au, surface.bu, surface.av, surface.bv
    tol = itol / 100.
    for t0, u0, v0 in tuv0_list:
        tuv0 = array([t0, u0, v0], dtype=float64)
        sol = minimize(_obj, tuv0, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        t, u, v = sol.x
        t = curve.check_param(t)
        u, v = surface.check_params(u, v)
        pc = ceval(t, rtype='ndarray', domain='global')
        ps = seval(u, v, rtype='ndarray', domain='global')
        if norm(pc - ps) > itol:
            continue
        pii = array(0.5 * (pc + ps), dtype=float64)
        unique = True
        for _, p12 in results:
            if norm(p12 - pii) <= itol:
                unique = False
                break
        if not unique:
            continue
        results.append([(t, u, v), pii])

    npts = len(results)
    return nsub[0], npts, results


def intersect_icurve_surface(icrv, surface, itol=None):
    """
    Find the intersection points of an intersection curve and a surface.

    :param icrv:
    :param surface:
    :param itol:
    :return:
    """
    if itol is None:
        itol = Settings.gtol

    # Intersect 3-D curve.
    itol_ = max(icrv.ftol, itol)
    nsub, npts, r3d = intersect_curve_surface(icrv.crv3d, surface, itol_)
    if npts == 0:
        return nsub, npts, r3d

    # Objective function for refinement.
    def _obj(x):
        factor = 1.
        # Penalize objective function if variable is outside domain.
        # if x[0] < a or x[0] > b:
        #     factor = 1000.
        # if x[1] < au or x[1] > bu:
        #     factor = 1000.
        # if x[2] < av or x[2] > bv:
        #     factor = 1000.
        pci = ceval(x[0], rtype='ndarray', domain='global')
        psi = seval(x[1], x[2], rtype='ndarray', domain='global')
        return norm(psi - pci) * factor

    ceval = icrv.eval
    seval = surface.eval
    crv_check = icrv.check_param
    srf_check = surface.check_params
    results = []
    # a, b = icrv.a, icrv.b
    # au, bu, av, bv = surface.au, surface.bu, surface.av, surface.bv
    tol = itol / 100.
    for t0, u0, v0 in [r[0] for r in r3d]:
        tuv0 = array([t0, u0, v0], dtype=float64)
        sol = minimize(_obj, tuv0, method='Nelder-Mead', tol=tol,
                       options={'ftol': tol})
        t, u, v = sol.x
        t = crv_check(t)
        u, v = srf_check(u, v)
        pc = ceval(t, rtype='ndarray', domain='global')
        ps = seval(u, v, rtype='ndarray', domain='global')
        if norm(pc - ps) > itol:
            continue
        pi = array(0.5 * (pc + ps), dtype=float64)
        unique = True
        for _, p12 in results:
            if norm(p12 - pi) <= itol:
                unique = False
                break
        if not unique:
            continue
        results.append([(t, u, v), pi])

    npts = len(results)
    return nsub, npts, results
