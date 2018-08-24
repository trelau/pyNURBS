from __future__ import division

from numpy import array, zeros, floor, float64, int32

from pynurbs.config import Settings
from pynurbs.geometry.methods.compare_floats import CompareFloats as CmpFlt


def bin_coeff(n, k):
    """
    Compute binomial coefficients.

    :param int n: Non-negative integer (*n* > 0).
    :param int k: Non-negative integer (*k* < *n*).

    :returns: Binomial coefficient.
    :rtype: float
    """
    if k < 0. or k > n:
        return 0.
    if k == 0 or k == n:
        return 1.
    k = min(k, n - k)
    c = 1.
    for i in range(k):
        c *= (n - i) / (i + 1)
    return c


def horner1(a, n, u):
    """
    Compute point on power basis curve.

    :param ndarray a: *n* x 3 array of points (x, y, z) where *n* is the
        number of points
    :param n: Degree of power basis curve.
    :param u: Parametric point to evaluate (0 <= u <= 1).

    :return: Point on power basis curve.
    :rtype: ndarray

    *Reference:* Algorithm A1.1 from "The NURBS Book".
    """
    u = CmpFlt.check_bounds(u, 0., 1.)
    pnt = a[n]
    for i in range(n - 1, -1, -1):
        pnt = pnt * u + a[i]
    return pnt


def all_bernstein(n, u):
    """
    Compute all *n* th-degree Bernstein polynomials.

    :param int n: Degree.
    :param float u: Parametric point (0 <= u <= 1).

    :return: Bernstein polynomials which are nonzero at *u*.
    :rtype: ndarray

    *Reference:* Algorithm A1.3 from "The NURBS Book".
    """
    u = CmpFlt.check_bounds(u, 0., 1.)
    bernstein = [0.0] * (n + 1)
    bernstein[0] = 1.0
    u1 = 1.0 - u
    for j in range(1, n + 1):
        saved = 0.0
        for k in range(0, j):
            temp = bernstein[k]
            bernstein[k] = saved + u1 * temp
            saved = u * temp
        bernstein[j] = saved
    return array(bernstein, dtype=float64)


def point_on_bezier_curve(cpw, n, u):
    """
    Compute point on Bezier curve.

    :param ndarray cpw: Control points.
    :param int n: Degree.
    :param u: Parametric point (0 <= u <= 1).

    :return: Point on Bezier curve.
    :rtype: ndarray

    *Reference:* Algorithm A1.4 from "The NURBS Book".
    """
    bernstein = all_bernstein(n, u)
    pnt = zeros(4, dtype=float64)
    for k in range(0, n + 1):
        pnt += bernstein[k] * cpw[k]
    return pnt


def deCasteljau1(cpw, n, u):
    """
    Compute point on a Bezier curve using deCasteljau algorithm.

    :param ndarray cpw: Control points.
    :param int n: Number of control points - 1.
    :param float u: Parametric point (0 <= u <= 1).

    :return: Point on Bezier curve.
    :rtype: ndarray

    *Reference:* Algorithm A1.5 from "The NURBS Book".
    """
    pnt = array(cpw, copy=True)
    u = CmpFlt.check_bounds(u, 0., 1.)
    for k in range(1, n + 1):
        for i in range(0, n - k + 1):
            pnt[i] = (1.0 - u) * pnt[i] + u * pnt[i + 1]
    return pnt[0]


def horner2(a, n, m, u, v):
    """
    Compute a point on a power basis surface.

    :param ndarray a: *n* x *m* array of points (x, y, z) where *n* is the
        number of points in the u-direction and *m* is the number of points
        in the v-direction.
    :param int n: Degree of power basis surface in u-direction.
    :param int m: Degree of power basis surface in v-direction.
    :param float u: Parametric point to evaluate (0 <= u <= 1).
    :param float v: Parametric point to evaluate (0 <= v <= 1).

    :return: Point on power basis surface.
    :rtype: ndarray

    *Reference:* Algorithm A1.6 from "The NURBS Book".
    """
    b = zeros((n + 1, 3), dtype=float64)
    for i in range(0, n + 1):
        b[i] = horner1(a[i], m, v)
    return horner1(b, n, u)


def deCasteljau2(cpw, n, m, u, v):
    """
    Compute a point on a Bezier surface.

    :param ndarray cpw: *n* x *m* array of points (x, y, z) where *n* is the
        number of points in the u-direction and *m* is the number of points
        in the v-direction.
    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param float u: Parametric point to evaluate (0 <= u <= 1).
    :param float v: Parametric point to evaluate (0 <= v <= 1).

    :return: Point on Bezier surface.
    :rtype: ndarray

    *Reference:* Algorithm A1.7 from "The NURBS Book".
    """
    if n <= m:
        temp = zeros((m + 1, 4), dtype=float64)
        for j in range(0, m + 1):
            # Swapped location of j-index compared to text since it appears
            # to be incorrect. The text shows the j-th row, but the number of
            # rows goes up to n, not m. Same applies to below.
            temp[j] = deCasteljau1(cpw[:, j], n, u)
        return deCasteljau1(temp, m, v)
    else:
        temp = zeros((n + 1, 4), dtype=float64)
        for i in range(0, n + 1):
            temp[i] = deCasteljau1(cpw[i, :], m, v)
    return deCasteljau1(temp, n, u)


def find_span(n, p, u, uk):
    """
    Determine the knot span index.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param float u: Parameter.
    :param ndarray uk: Knot vector.

    :return: Knot span.
    :rtype: int

    *Reference:* Algorithm A2.1 from "The NURBS Book".
    """
    ptol = Settings.ptol

    # Special case
    # if u >= uk[n + 1]:
    if CmpFlt.ge(u, uk[n + 1], ptol):
        return n
    # if u <= uk[p]:
    if CmpFlt.le(u, uk[p], ptol):
        return p
    # Do binary search
    low = p
    high = n + 1
    mid = int(floor((low + high) / 2.))
    # while u < uk[mid] or u >= uk[mid + 1]:
    while CmpFlt.lt(u, uk[mid], ptol) or CmpFlt.ge(u, uk[mid + 1], ptol):
        # if u < uk[mid]:
        if CmpFlt.lt(u, uk[mid], ptol):
            high = mid
        else:
            low = mid
        mid = int(floor((low + high) / 2.))
    return mid


def basis_funs(i, u, p, uk):
    """
    Compute the non-vanishing basis functions.

    :param int i: Knot span index.
    :param float u: Parameter.
    :param int p: Degree.
    :param ndarray uk: Knot vector.

    :return: Non-vanishing basis functions.
    :rtype: ndarray

    *Reference:* Algorithm A2.2 from "The NURBS Book"
    """
    u = CmpFlt.check_bounds(u, uk[0], uk[-1])
    bf = [0.0] * (p + 1)
    bf[0] = 1.0
    left = [0.0] * (p + 1)
    right = [0.0] * (p + 1)
    for j in range(1, p + 1):
        left[j] = u - uk[i + 1 - j]
        right[j] = uk[i + j] - u
        saved = 0.0
        for r in range(0, j):
            temp = bf[r] / (right[r + 1] + left[j - r])
            bf[r] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        bf[j] = saved
    return array(bf, dtype=float64)


def ders_basis_funs(i, u, p, n, uk):
    """
    Compute nonzero basis functions and their derivatives. First section is
    A2.2 modified to store functions and knot differences.

    :param int i: Knot span index.
    :param float u: Parameter.
    :param int p: Degree.
    :param int n: Number of derivatives to compute (*n* <= *p*).
    :param ndarray uk: Knot vector.

    :return: Two dimensional array ders[k][j] where k is the kth derivative.
    :rtype: ndarray

    *Reference:* Algorithm A2.3 from "The NURBS Book"
    """
    u = CmpFlt.check_bounds(u, uk[0], uk[-1])
    ndu = [[0. for _ in range(p + 1)] for _ in range(p + 1)]
    ndu[0][0] = 1.0
    left = [0.0] * (p + 1)
    right = [0.0] * (p + 1)
    for j in range(1, p + 1):
        left[j] = u - uk[i + 1 - j]
        right[j] = uk[i + j] - u
        saved = 0.0
        for r in range(0, j):
            ndu[j][r] = right[r + 1] + left[j - r]
            temp = ndu[r][j - 1] / ndu[j][r]
            ndu[r][j] = saved + right[r + 1] * temp
            saved = left[j - r] * temp
        ndu[j][j] = saved
    ders = [[0. for _ in range(p + 1)] for _ in range(n + 1)]
    a = [[0. for _ in range(p + 1)] for _ in range(2)]
    for j in range(0, p + 1):
        ders[0][j] = ndu[j][p]
    for r in range(0, p + 1):
        s1 = 0
        s2 = 1
        a[0][0] = 1.
        for k in range(1, n + 1):
            d = 0.
            rk = r - k
            pk = p - k
            if r >= k:
                a[s2][0] = a[s1][0] / ndu[pk + 1][rk]
                d = a[s2][0] * ndu[rk][pk]
            if rk >= -1:
                j1 = 1
            else:
                j1 = -rk
            if r - 1 <= pk:
                j2 = k - 1
            else:
                j2 = p - r
            for j in range(j1, j2 + 1):
                a[s2][j] = (a[s1][j] - a[s1][j - 1]) / ndu[pk + 1][rk + j]
                d += a[s2][j] * ndu[rk + j][pk]
            if r <= pk:
                a[s2][k] = -a[s1][k - 1] / ndu[pk + 1][r]
                d += a[s2][k] * ndu[r][pk]
            ders[k][r] = d
            j = s1
            s1 = s2
            s2 = j
    r = p
    for k in range(1, n + 1):
        for j in range(0, p + 1):
            ders[k][j] *= r
        r *= (p - k)
    return array(ders, dtype=float64)


def curve_point(n, p, uk, cpw, u):
    """
    Compute curve point.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: knot sequence.
    :param ndarray cpw: Control points.
    :param float u: Parametric point.

    :return: Point on curve.
    :rtype: ndarray

    *Reference:* Algorithm A3.1 from "The NURBS Book"
    """
    span = find_span(n, p, u, uk)
    bf = basis_funs(span, u, p, uk)
    pnt = zeros(4, dtype=float64)
    for i in range(0, p + 1):
        pnt += bf[i] * cpw[span - p + i]
    return pnt[:-1] / pnt[-1]


def curve_derivs_alg1(n, p, uk, cpw, u, d):
    """
    Compute curve derivatives.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Array of control points.
    :param float u: Parametric point.
    :param int d: The maximum derivative to compute.

    :return: The curve derivatives in array ck[], where ck[k] is the *k*-th
      derivative and 0 <= k <= d.
    :rtype: ndarray

    *Reference:* Algorithm A3.2 from "The NURBS Book"
    """
    du = min(d, p)
    ck = zeros((d + 1, 4), dtype=float64)
    span = find_span(n, p, u, uk)
    nders = ders_basis_funs(span, u, p, d, uk)
    for k in range(0, du + 1):
        for j in range(0, p + 1):
            ck[k] += nders[k, j] * cpw[span - p + j]
    return ck


def curve_deriv_cpts(n, p, uk, cpw, d, r1=None, r2=None):
    """
    Compute control points of curve derivatives.

    :param  int n: Number of control points - 1.
    :param int p: Degree.
    :param uk: Knot vector.
    :param cpw: Control points.
    :param d: The maximum derivative to compute.
    :param r1: Lower bound (r1 = 0 if None provided).
    :param r2: Upper bound (r2 = n if None provided).

    :return: Control points of curve derivatives.
    :rtype: ndarray

    *Reference:* Algorithm A3.3 from "The NURBS Book"
    """
    if r1 is None:
        r1 = 0
    if r2 is None:
        r2 = n
    r = r2 - r1
    pk = zeros((d + 1, r + 1, 4), dtype=float64)
    for i in range(r + 1):
        pk[0, r] = cpw[r1 + i]
    for k in range(1, d + 1):
        temp = p - k + 1
        for i in range(0, r - k + 1):
            pk[k, i] = temp * (pk[k - 1, i + 1] - pk[k - 1, i]) / (
                    uk[r1 + i + p + 1] - uk[r1 + i + k])
    return pk


def surface_point(n, p, uk, m, q, vk, cpw, u, v):
    """
    Compute surface point.

    :param int n: Number of control points - 1 for u-direction.
    :param int p: Degree for u-direction.
    :param ndarray uk: Knot vector for u-direction.
    :param int m: Number of control points - 1 for v-direction.
    :param int q: Degree for v-direction.
    :param ndarray vk: Knot vector for v-direction.
    :param ndarray cpw: *n* x *m* array of control points.
    :param float u: Parametric point for u-direction.
    :param float v: Parametric point for v-direction.

    :return: Point on surface.
    :rtype: ndarray

    *Reference:* Algorithm A3.5 from "The NURBS Book"
    """
    uspan = find_span(n, p, u, uk)
    nu = basis_funs(uspan, u, p, uk)
    vspan = find_span(m, q, v, vk)
    nv = basis_funs(vspan, v, q, vk)
    uind = uspan - p
    pnt = zeros(4, dtype=float64)
    for l in range(0, q + 1):
        temp = zeros(4, dtype=float64)
        vind = vspan - q + l
        for k in range(0, p + 1):
            temp += nu[k] * cpw[uind + k, vind]
        pnt += nv[l] * temp
    return pnt[:-1] / pnt[-1]


def surface_derivs_alg1(n, p, uk, m, q, vk, cpw, u, v, d):
    """
    Compute surface derivatives.

    :param int n: Number of control points - 1 for u-direction.
    :param int p: Degree for u-direction.
    :param ndarray uk: Knot vector for u-direction.
    :param int m: Number of control points - 1 for v-direction.
    :param int q: Degree for v-direction.
    :param ndarray vk: Knot vector for v-direction.
    :param ndarray cpw: *n* x *m* array of control points.
    :param float u: Parametric point for u-direction.
    :param float v: Parametric point for v-direction.
    :param int d: Maximum derivative to calculate.

    :return: Array skl[][], where skl[k][l] is the derivative of S(u,v)
        with respect to *u* k times, and *v* l times and 0 <= k <= d.
    :rtype: ndarray

    *Reference:* Algorithm A3.6 from "The NURBS Book"
    """
    du = min(d, p)
    dv = min(d, q)
    skl = zeros((d + 1, d + 1, 4), dtype=float64)
    uspan = find_span(n, p, u, uk)
    nu = ders_basis_funs(uspan, u, p, du, uk)
    vspan = find_span(m, q, v, vk)
    nv = ders_basis_funs(vspan, v, q, dv, vk)
    temp = zeros((q + 1, 4), dtype=float64)
    for k in range(0, du + 1):
        for s in range(0, q + 1):
            temp[s] = [0., 0., 0., 0.]
            for r in range(0, p + 1):
                temp[s] += nu[k, r] * cpw[uspan - p + r, vspan - q + s]
        dd = min(d - k, dv)
        for l in range(0, dd + 1):
            for s in range(0, q + 1):
                skl[k, l] += nv[l, s] * temp[s]
    return skl


def rat_curve_derivs(n, p, uk, cpw, u, d):
    """
    Compute rational curve derivatives.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Array of control points.
    :param float u: Parametric point.
    :param int d: The maximum derivative to compute.

    :return: The curve derivatives in array ck[], where ck[k] is the *k*-th
      derivative and 0 <= k <= d.
    :rtype: ndarray

    *Reference:* Algorithm A4.2 from "The NURBS Book"
    """
    cders = curve_derivs_alg1(n, p, uk, cpw, u, d)
    aders, wders = cders[:, :-1], cders[:, -1]
    ck = zeros((d + 1, 3), dtype=float64)
    for k in range(0, d + 1):
        v = aders[k]
        for i in range(1, k + 1):
            binc = bin_coeff(k, i)
            v -= binc * wders[i] * ck[k - i]
        ck[k] = v / wders[0]
    return ck


def rat_surface_derivs(n, p, uk, m, q, vk, cpw, u, v, d):
    """
    Compute surface derivatives.

    :param int n: Number of control points - 1 for u-direction.
    :param int p: Degree for u-direction.
    :param ndarray uk: Knot vector for u-direction.
    :param int m: Number of control points - 1 for v-direction.
    :param int q: Degree for v-direction.
    :param ndarray vk: Knot vector for v-direction.
    :param ndarray cpw: *n* x *m* array of control points.
    :param float u: Parametric point for u-direction.
    :param float v: Parametric point for v-direction.
    :param int d: Maximum derivative to calculate.

    :return: Array skl[][], where skl[k][l] is the derivative of S(u,v)
        with respect to *u* k times, and *v* l times.
    :rtype: ndarray

    *Reference:* Algorithm A4.4 from "The NURBS Book"
    """
    sders = surface_derivs_alg1(n, p, uk, m, q, vk, cpw,
                                u, v, d)
    aders, wders = sders[:, :, :-1], sders[:, :, -1]
    skl = zeros((d + 1, d + 1, 3), dtype=float64)
    for k in range(0, d + 1):
        for l in range(0, d - k + 1):
            v1 = aders[k, l]
            for j in range(1, l + 1):
                v1 -= bin_coeff(l, j) * wders[0, j] * skl[k, l - j]
            for i in range(1, k + 1):
                v1 -= bin_coeff(k, i) * wders[i, 0] * skl[k - i, l]
                v2 = zeros(3, dtype=float64)
                for j in range(1, l + 1):
                    v2 += bin_coeff(l, j) * wders[i, j] * skl[k - i, l - j]
                v1 -= bin_coeff(k, i) * v2
            skl[k, l] = v1 / wders[0, 0]
    return skl


def find_span_mult(n, p, u, uk):
    """
    Determine the knot span index and multiplicity.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param float u: Parameter.
    :param ndarray uk: Knot vector.

    :return: Knot span and multiplicity.
    :rtype: tuple
    """
    k = find_span(n, p, u, uk)
    ptol = Settings.ptol
    # if uk[k] < u < uk[k + 1]:
    if CmpFlt.gt(u, uk[k], ptol) and CmpFlt.lt(u, uk[k + 1], ptol):
        return k, 0
    # if uk[k] < u <= uk[k + 1]:
    if CmpFlt.gt(u, uk[k], ptol) and CmpFlt.le(u, uk[k + 1], ptol):
        first = k + 1
        last = n + p + 1
        step = 1
    else:
        first = k
        last = 0
        step = -1
    s = 0
    for i in range(first, last + step, step):
        # if uk[i] == u:
        if CmpFlt.eq(uk[i], u, ptol):
            s += 1
        else:
            break
    return k, s


def find_mult(n, p, uk, u):
    """
    Determine the multiplicity of the parameter in the knot vector.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param float u: Parameter.

    :return: Multiplicity of parameter.
    :rtype: int
    """
    _, s = find_span_mult(n, p, u, uk)
    return s


def find_mult_knots(n, p, uk):
    """
    Find the multiplicities and the unique knots of the knot vector. This
    includes the entire knot vector (not just the interior knots).

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.

    :return: Number of unique knots, multiplicities, and the knots
        (nu, um, uq).
    :rtype: tuple
    """
    ptol = Settings.ptol

    i = 0
    m = n + p + 1
    nu = 0
    uq = zeros(m + 1, dtype=float64)
    um = zeros(m + 1, dtype=int32)
    while i <= m:
        uq[nu] = uk[i]
        mult = 0
        # while i <= m and uk[i] == uq[nu]:
        while i <= m and CmpFlt.eq(uk[i], uq[nu], ptol):
            i += 1
            mult += 1
        um[nu] = mult
        nu += 1
    return nu, um[:nu], uq[:nu]


def check_param(n, p, uk, u, is_closed=False):
    """
    Check that the parameter is within the global domain of the knot vector or
    is within tolerance of a unique knot value. Use :class:`.CompareFloats` to
    compare floats.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param float u: Global parameter.
    :param bool is_closed: Option to specify closed geometry so the knot
        value will wrap around to start/end value.

    :return: Parameter within global domain or near interior knot value.
    :rtype: float
    """
    ptol = Settings.ptol
    if is_closed:
        if u > uk[n + 1]:
            u = u - uk[n + 1] + uk[p]
        elif u < uk[p]:
            u = uk[n + 1] + u - uk[p]

    k = find_span(n, p, u, uk)
    if CmpFlt.le(u, uk[k], ptol):
        return uk[k]
    if CmpFlt.ge(u, uk[k + 1], ptol):
        return uk[k + 1]
    return u


def bezier_curve_points(n, cpw, ulist):
    """
    Compute Bezier curve points.

    :param int n: Number of control points - 1.
    :param ndarray cpw: Control points.
    :param array_like ulist: Parametric points.

    :return: Points Bezier on curve.
    :rtype: ndarray
    """
    arr_u = array(ulist, dtype=float64)
    nu = arr_u.shape[0]

    pnts = zeros((nu, 3), dtype=float64)
    for i in range(0, nu):
        pi = deCasteljau1(cpw, n, arr_u[i])
        pnts[i, :] = pi[:-1] / pi[-1]

    return pnts


def curve_points(n, p, uk, cpw, ulist):
    """
    Compute curve points.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: knot sequence.
    :param ndarray cpw: Control points.
    :param array_like ulist: Parameters.

    :return: Points on curve.
    :rtype: ndarray
    """
    arr_u = array(ulist, dtype=float64)
    nu = arr_u.shape[0]

    pnts = zeros((nu, 3), dtype=float64)
    for i in range(0, nu):
        pi = curve_point(n, p, uk, cpw, arr_u[i])
        pnts[i, :] = pi

    return pnts


def bezier_surface_points(n, m, cpw, ulist, vlist):
    """
    Compute Bezier surface points.

    :param int n: Number of control points - 1 for u-direction.
    :param int m: Number of control points - 1 for v-direction.
    :param ndarray cpw: Control points.
    :param array_like ulist: Parameters in u-direction.
    :param array_like vlist: Parameters in v-direction.

    :return: Points on surface.
    :rtype: ndarray
    """
    arr_u = array(ulist, dtype=float64)
    nu = arr_u.shape[0]
    arr_v = array(vlist, dtype=float64)
    nv = arr_v.shape[0]
    np = min(nu, nv)

    pnts = zeros((np, 3), dtype=float64)
    for i in range(0, np):
        pi = deCasteljau2(cpw, n, m, arr_u[i], arr_v[i])
        pnts[i, :] = pi[:-1] / pi[-1]

    return pnts


def surface_points(n, p, uk, m, q, vk, cpw, ulist, vlist):
    """
    Compute surface points.

    :param int n: Number of control points - 1 for u-direction.
    :param int p: Degree for u-direction.
    :param ndarray uk: Knot vector for u-direction.
    :param int m: Number of control points - 1 for v-direction.
    :param int q: Degree for v-direction.
    :param ndarray vk: Knot vector for v-direction.
    :param ndarray cpw: Control points.
    :param array_like ulist: Parameters in u-direction.
    :param array_like vlist: Parameters in v-direction.

    :return: Points on surface.
    :rtype: ndarray
    """
    arr_u = array(ulist, dtype=float64)
    nu = arr_u.shape[0]
    arr_v = array(vlist, dtype=float64)
    nv = arr_v.shape[0]
    np = min(nu, nv)

    pnts = zeros((np, 3), dtype=float64)
    for i in range(0, np):
        pi = surface_point(n, p, uk, m, q, vk, cpw, arr_u[i], arr_v[i])
        pnts[i, :] = pi

    return pnts
