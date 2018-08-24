from __future__ import division

from numpy import (ndarray, array, float64, zeros, int32, amax, mean, dot,
                   identity)
from scipy.linalg import lu_factor, lu_solve

from pynurbs.config import Settings
from pynurbs.geometry.methods.compare_floats import CompareFloats as CmpFlt
from pynurbs.geometry.methods.evaluate import (find_span, basis_funs,
                                               find_span_mult)
from pynurbs.geometry.methods.parameterize import (uniform, chord_length,
                                                   centripetal)


def global_curve_interpolation(qp, p=3, method='chord'):
    """
    Global curve interpolation through n + 1 points.

    :param qp: List or array of points.
    :etype qp: Points or array_like
    :param int p: Degree. Will be adjusted if not enough interpolation points
        are provided (p = npoints - 1).
    :param str method: Option to specify method for selecting parametersexit
        ("uniform", "chord", or "centripetal").

    :return: NURBS curve data interpolating points (cp, uk, p).
    :rtype: tuple

    *Reference:* "The NURBS Book" Section 9.2.1.
    """
    # Get qp as array if not already.
    if not isinstance(qp, ndarray):
        qp = array(qp, dtype=float64)

    # Determine n from shape of qp.
    n = qp.shape[0] - 1

    # Check that enough points are provided for desired degree.
    if n < p:
        p = n

    # Compute parameters between [0, 1].
    m = n + p + 1
    if method.lower() in ['u', 'uniform']:
        u = uniform(qp, 0., 1.)
    elif method.lower() in ['ch', 'chord']:
        u = chord_length(qp, 0., 1.)
    else:
        u = centripetal(qp, 0., 1.)

    # Compute knot vector by averaging.
    uk = zeros(m + 1, dtype=float64)
    uk[m - p:] = 1.0
    for j in range(1, n - p + 1):
        temp = 0.
        for i in range(j, j + p):
            temp += u[i]
        uk[j + p] = 1.0 / p * temp

    # Set up system of linear equations.
    if p > 1:
        a = zeros((n + 1, n + 1), dtype=float64)
        for i in range(0, n + 1):
            span = find_span(n, p, u[i], uk)
            a[i, span - p: span + 1] = basis_funs(span, u[i], p, uk)
    else:
        a = identity(n + 1, dtype=float64)

    # Solve for [a][cp] = [qp] using LU decomposition.
    lu, piv = lu_factor(a, overwrite_a=True, check_finite=False)
    cp = lu_solve((lu, piv), qp, trans=0, overwrite_b=True, check_finite=True)

    return cp, uk, p


# def global_surface_interpolation():
#     """
#     Global surface interpolation through (n + 1) x (m + 1) points.
#     """
#     pass

def interpolate_curves(curves, q=3, method='chord', inplace=False,
                       auto_reverse=True):
    """
    Create a surface by interpolating the list of section curves.

    :param list curves: List of NURBS curves to skin. Order of the curves in
        the list will determine the longitudinal (*v*) direction.
    :param int q: Degree in v-direction. Will be adjusted if not enough
        skinning curves are provided (q = ncurves - 1).
    :param str method: Option to specify method for selecting parameters
        ("uniform", "chord", or "centripetal").
    :param bool inplace: Option to modify the curves in the list in-place. If
        this option is *False*, then new curves will be created. If this option
        is *True*, then the curves in the list will be modified in-place.
    :param bool auto_reverse: Option to check direction of curves and
        reverse if necessary.

    :return: NURBS surface data (cpw, uk, vk, p, q).
    :rtype: tuple

    *Reference:* "The NURBS Book" Section 10.3.
    """
    # Step 0: Adjust desired degree in case not enough curves are provided.
    if len(curves) - 1 < q:
        q = len(curves) - 1

    # Step 1: Make sure each curve in a NURBS curve.
    nurbs_curves = []
    for c in curves:
        if inplace:
            nurbs_curves.append(c)
        else:
            nurbs_curves.append(c.copy())
    if not nurbs_curves:
        return None

    # Check direction of curves by comparing the derivatives at the midpoint.
    if auto_reverse:
        cu_0 = nurbs_curves[0].deriv(0.5, 1, rtype='ndarray')
        for c in nurbs_curves[1:]:
            cu_i = c.deriv(0.5, 1, rtype='ndarray')
            if dot(cu_0, cu_i) < 0.:
                c.reverse()
            cu_0 = cu_i

    # Step 2: Make sure each curve has a common degree.
    pmax = 0
    for c in nurbs_curves:
        if c.p > pmax:
            pmax = c.p
    for c in nurbs_curves:
        if c.p < pmax:
            c.elevate_degree(pmax, inplace=True)

    # Step 3: Make sure each curve has a similar knot vector.
    # Find each unique interior knot value in every curve.
    a, b = [], []
    for c in nurbs_curves:
        a.append(c.a)
        b.append(c.b)
    amin = min(a)
    bmax = max(b)
    for c in nurbs_curves:
        c.set_domain(amin, bmax)
    all_knots = []
    for c in nurbs_curves:
        nu, _, knots = c.get_mult_and_knots()
        for uk in knots[:nu]:
            all_knots.append(uk)
    all_knots.sort()
    unique_knots = []
    m = len(all_knots) - 1
    i = 0
    while i <= m:
        unique_knots.append(all_knots[i])
        # while i <= m and all_knots[i] == unique_knots[-1]:
        while i <= m and CmpFlt.eq(all_knots[i], unique_knots[-1],
                                   Settings.ptol):
            i += 1
    # For each curve, find the multiplicity and then the max multiplicity of
    # each unique knot value.
    all_mult = []
    for c in nurbs_curves:
        row = []
        for ui in unique_knots:
            _, mult = find_span_mult(c.n, c.p, ui, c.uk)
            row.append(mult)
        all_mult.append(row)
    mult_arr = array(all_mult, dtype=int32)
    max_mult = amax(mult_arr, axis=0)
    # For each curve, insert the knot the required number of times, if any.
    for mult, ui in zip(max_mult, unique_knots):
        for c in nurbs_curves:
            _, c_mult = find_span_mult(c.n, c.p, ui, c.uk)
            if c_mult < mult:
                c.insert_knot(ui, mult - c_mult, inplace=True, domain='global')

    # Step 4: Compute v-direction parameters between [0, 1].
    # Find parameters between each curve by averaging each segment.
    temp = array([c.cpw for c in nurbs_curves], dtype=float64)
    pnts_matrix = temp.transpose((1, 0, 2))
    c0 = nurbs_curves[0]
    n = c0.n
    m = len(nurbs_curves) - 1
    v_matrix = zeros((n + 1, m + 1), dtype=float64)
    for i in range(0, n + 1):
        if method.lower() in ['u', 'uniform']:
            v = uniform(pnts_matrix[i, :], 0., 1.)
        elif method.lower() in ['ch', 'chord']:
            v = chord_length(pnts_matrix[i, :], 0., 1.)
        else:
            v = centripetal(pnts_matrix[i, :], 0., 1.)
        v_matrix[i] = v
    # Average each column.
    v = mean(v_matrix, axis=0, dtype=float64)
    v[0] = 0.0
    v[-1] = 1.0

    # Step 5: Compute knot vector by averaging.
    uk = c0.uk
    s = m + q + 1
    vk = zeros(s + 1, dtype=float64)
    vk[s - q:] = 1.0
    for j in range(1, m - q + 1):
        temp = 0.
        for i in range(j, j + q):
            temp += v[i]
        vk[j + q] = 1.0 / q * temp

    # Step 6: Perform n + 1 interpolations for v-direction.
    cpw = zeros((n + 1, m + 1, 4), dtype=float64)
    for i in range(0, n + 1):
        qp = pnts_matrix[i]
        # Set up system of linear equations.
        a = zeros((m + 1, m + 1), dtype=float64)
        for j in range(0, m + 1):
            span = find_span(m, q, v[j], vk)
            a[j, span - q: span + 1] = basis_funs(span, v[j], q, vk)
        # Solve for [a][cp] = [qp] using LU decomposition.
        lu, piv = lu_factor(a, overwrite_a=True, check_finite=False)
        cpw[i, :] = lu_solve((lu, piv), qp, trans=0, overwrite_b=True,
                             check_finite=True)

    return cpw, uk, vk, pmax, q
