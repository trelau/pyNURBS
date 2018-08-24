from __future__ import division

from numpy import zeros, float64, array, insert, floor

from pynurbs.config import Settings
from pynurbs.geometry.methods.compare_floats import CompareFloats as CmpFlt
from pynurbs.geometry.methods.evaluate import (find_span_mult, find_span,
                                               bin_coeff, find_mult_knots)


def adjust_knot_vector(n, p, uk, a=0., b=1.):
    """
    Adjust the knot vector *uk* to be between [a. b], where a < b. The
    interior knots are adjusted to maintain original parametrization.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param float a: Lower bound.
    :param float b: Upper bound.

    :return: New knot vector *qk*.
    :rtype: ndarray
    """
    if a > b:
        return uk
    a0 = uk[p]
    b0 = uk[n + 1]
    ab0 = b0 - a0
    ab = b - a
    qk = zeros(n + p + 2, dtype=float64)
    qk[:p + 1] = a
    qk[n + 1:] = b
    for i in range(p + 1, n + 1):
        u = uk[i]
        qk[i] = a + (u - a0) * ab / ab0
    return qk


def elevate_bezier_curve_degree(n, deg, cpw):
    """
    Elevate degree of Bezier curve without changing its shape.

    :param int n: Number of control points - 1 (degree).
    :param int deg: Degree to elevate curve to.
    :param ndarray cpw: Control points.

    :return: Control points of Bezier curve with elevated degree.
    :rtype ndarray
    """
    qw = array(cpw)
    for ni in range(n, deg - n + 2):
        qpw_p1 = insert(qw, -1, [0., 0., 0., 1.], axis=0)
        for j in range(1, ni + 1):
            pj = qw[j]
            p_jm1 = qw[j - 1]
            qpw_p1[j] = (j / (ni + 1)) * p_jm1 + (1 - j / (ni + 1)) * pj
        qw = array(qpw_p1)
    return qw


def curve_knot_ins(n, p, uk, cpw, u, r):
    """
    Curve knot insertion.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.
    :param float u: Knot value to insert.
    :param int r: Number of times to insert knot.

    :return: Knot vector and control points after knot insertion
        (knots, cpw).
    :rtype: tuple

    *Reference:* Algorithm A5.1 from "The NURBS Book."
    """
    # Find span and multiplicity.
    k, s = find_span_mult(n, p, u, uk)
    if r + s > p:
        r = p - s

    # Generate knot vector
    nq = n + r
    uq = zeros(nq + p + 2, dtype=float64)
    for i in range(0, k + 1):
        uq[i] = uk[i]
    for i in range(1, r + 1):
        uq[k + i] = u
    for i in range(k + 1, n + p + 2):
        uq[i + r] = uk[i]
    # Save unchanged control points.
    qw = zeros((nq + 1, 4), dtype=float64)
    for i in range(0, k - p + 1):
        qw[i] = cpw[i]
    for i in range(k - s, n + 1):
        qw[i + r] = cpw[i]
    rw = zeros((p + 1, 4), dtype=float64)
    for i in range(0, p - s + 1):
        rw[i] = cpw[k - p + i]
    # Insert knot r times.
    t = k - p
    for j in range(1, r + 1):
        t = k - p + j
        for i in range(0, p - j - s + 1):
            alpha = (u - uk[t + i]) / (uk[i + k + 1] - uk[t + i])
            rw[i] = alpha * rw[i + 1] + (1.0 - alpha) * rw[i]
        qw[t] = rw[0]
        qw[k + r - j - s] = rw[p - j - s]
    # Load remaining control points.
    for i in range(t + 1, k - s):
        qw[i] = rw[i - t]
    return uq, qw


# def curve_pnt_by_corner_cut(n, p, uk, cpw, u):
#     """
#     Compute point on rational B-spline curve by corner cutting.
#
#     :param int n: Number of control points - 1.
#     :param int p: Degree.
#     :param ndarray uk: Knot vector.
#     :param ndarray cpw: Control points.
#     :param float u: Parametric point.
#
#     :return: Point on curve.
#     :rtype: ndarray
#
#     *Reference:* Algorithm A5.2 "The NURBS Book."
#     """
#     # Endpoints are special cases.
#     if u == uk[0]:
#         return cpw[0, :-1] / cpw[0, -1]
#     if u == uk[-1]:
#         return cpw[-1, :-1] / cpw[-1, -1]
#     # General case.
#     k, s = find_span_mult(n, p, u, uk)
#     r = p - s
#     rw = zeros((p + 1, 4), dtype=float64)
#     for i in range(0, r + 1):
#         rw[i] = cpw[k - p + i]
#     for j in range(1, r + 1):
#         for i in range(0, r - j + 1):
#             alpha = ((u - uk[k - p + j + i]) /
#                      (uk[i + k + 1] - uk[k - p + j + i]))
#             rw[i] = alpha * rw[i + 1] + (1.0 - alpha) * rw[i]
#     return rw[0, :-1] / rw[0, -1]


def surface_knot_ins(n, p, uk, m, q, vk, cpw, d, uv, r):
    """
    Surface knot insertion.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param str d: Direction to insert knot *uv* ("u": u-direction,
        "v": v-direction.
    :param float uv: Knot to insert.
    :param int r: Number of times to insert knot.

    :return: Knot vector and control points after knot insertion (uq, vq, cpw).
    :rtype: tuple

    *Reference:* Algorithm A5.3 from "The NURBS Book."
    """
    if d.lower() in ['u']:
        # Generate u-direction knot vector.
        k, s = find_span_mult(n, p, uv, uk)
        if r + s > p:
            r = p - s
        nq = n + r
        uq = zeros(nq + p + 2, dtype=float64)
        for i in range(0, k + 1):
            uq[i] = uk[i]
        for i in range(1, r + 1):
            uq[k + i] = uv
        for i in range(k + 1, n + p + 2):
            uq[i + r] = uk[i]
        # Copy v-direction knot vector.
        vq = array(vk)
        # Save alphas.
        alpha = zeros((p + 1, r + 1), dtype=float64)
        for j in range(1, r + 1):
            t = k - p + j
            for i in range(0, p - j - s + 1):
                alpha[i, j] = (uv - uk[t + i]) / (uk[i + k + 1] - uk[t + i])
        # Insert knot for each row at each column.
        qw = zeros((nq + 1, m + 1, 4), dtype=float64)
        for col in range(0, m + 1):
            # Save unchanged control points.
            for i in range(0, k - p + 1):
                qw[i, col] = cpw[i, col]
            for i in range(k - s, n + 1):
                qw[i + r, col] = cpw[i, col]
            # Load auxiliary control points.
            rw = zeros((p + 1, 4), dtype=float64)
            for i in range(0, p - s + 1):
                rw[i] = cpw[k - p + i, col]
            # Insert the knot r times.
            t = k - p
            for j in range(1, r + 1):
                t = k - p + j
                for i in range(0, p - j - s + 1):
                    rw[i] = (alpha[i, j] * rw[i + 1] +
                             (1.0 - alpha[i, j]) * rw[i])
                qw[t, col] = rw[0]
                qw[k + r - j - s, col] = rw[p - j - s]
            # Load remaining control points.
            for i in range(t + 1, k - s):
                qw[i, col] = rw[i - t]
        return uq, vq, qw
    if d.lower() in ['v']:
        # Generate v-direction knot vector.
        k, s = find_span_mult(m, q, uv, vk)
        if r + s > q:
            r = q - s
        mq = m + r
        vq = zeros(mq + q + 2, dtype=float64)
        for i in range(0, k + 1):
            vq[i] = vk[i]
        for i in range(1, r + 1):
            vq[k + i] = uv
        for i in range(k + 1, m + q + 2):
            vq[i + r] = vk[i]
        # Copy u-direction knot vector.
        uq = array(uk)
        # Save alphas.
        alpha = zeros((q + 1, r + 1), dtype=float64)
        for j in range(1, r + 1):
            t = k - q + j
            for i in range(0, q - j - s + 1):
                alpha[i, j] = (uv - vk[t + i]) / (vk[i + k + 1] - vk[t + i])
        # Insert knot for each column at each row.
        qw = zeros((n + 1, mq + 1, 4), dtype=float64)
        for row in range(0, n + 1):
            # Save unchanged control points.
            for i in range(0, k - q + 1):
                qw[row, i] = cpw[row, i]
            for i in range(k - s, m + 1):
                qw[row, i + r] = cpw[row, i]
            # Load auxiliary control points.
            rw = zeros((q + 1, 4), dtype=float64)
            for i in range(0, q - s + 1):
                rw[i] = cpw[row, k - q + i]
            # Insert the knot r times.
            t = k - q
            for j in range(1, r + 1):
                t = k - q + j
                for i in range(0, q - j - s + 1):
                    rw[i] = (alpha[i, j] * rw[i + 1] +
                             (1.0 - alpha[i, j]) * rw[i])
                qw[row, t] = rw[0]
                qw[row, k + r - j - s] = rw[q - j - s]
            # Load remaining control points.
            for i in range(t + 1, k - s):
                qw[row, i] = rw[i - t]
        return uq, vq, qw


def filter_knot_vect(n, p, uk, x):
    """
    Filter the knot vector *x* so that the knot multiplicity is less than or
    equal to the degree *p* when inserted into *uk*.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector (sorted).
    :param ndarray x: Knots to add to existing knot vector (sorted).

    :return: Filtered knot vector *x*.
    :rtype: ndarray
    """
    temp = array(uk)
    xnew = []
    for i in range(0, len(x)):
        span, s = find_span_mult(n, p, x[i], temp)
        if s <= p:
            temp = insert(temp, span, x[i])
            n += 1
            xnew.append(x[i])
    return array(xnew, dtype=float64)


def refine_knot_vect_curve(n, p, uk, cpw, x):
    """
    Refine curve knot vector.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.
    :param ndarray x: Knots to add to existing knot vector.

    :return: New knot vector and control points (ubar, qw).
    :rtypeL tuple

    *Reference:* Algorithm A5.4 from "The NURBS Book."
    """
    x = filter_knot_vect(n, p, uk, x)
    m = n + p + 1
    r = len(x) - 1
    a = find_span(n, p, x[0], uk)
    b = find_span(n, p, x[r], uk) + 1
    nq = m + (r + 1) - p - 1
    qw = zeros((nq + 1, 4), dtype=float64)
    # Load unchanged control points.
    for j in range(0, a - p + 1):
        qw[j] = cpw[j]
    for j in range(b - 1, n + 1):
        qw[j + r + 1] = cpw[j]
    # Load unchanged knot vector.
    ubar = zeros(m + (r + 1) + 1, dtype=float64)
    for j in range(0, a + 1):
        ubar[j] = uk[j]
    for j in range(b + p, m + 1):
        ubar[j + r + 1] = uk[j]
    # Find new control points.
    i = b + p - 1
    k = b + p + r
    for j in range(r, -1, -1):
        while x[j] <= uk[i] and i > a:
            qw[k - p - 1] = cpw[i - p - 1]
            ubar[k] = uk[i]
            k -= 1
            i -= 1
        qw[k - p - 1] = qw[k - p]
        for t in range(1, p + 1):
            ind = k - p + t
            alpha = ubar[k + t] - x[j]
            if abs(alpha) == 0.0:
                qw[ind - 1] = qw[ind]
            else:
                alpha /= (ubar[k + t] - uk[i - p + t])
                qw[ind - 1] = alpha * qw[ind - 1] + (1.0 - alpha) * qw[ind]
        ubar[k] = x[j]
        k -= 1
    return ubar, qw


def refine_knot_vect_surface(n, p, uk, m, q, vk, cpw, x, d):
    """
    Refine surface knot vector.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param ndarray x: Knots to add to existing knot vector.
    :param str d: Direction to insert knot *uv* ("u": u-direction,
        "v": v-direction.

    :return: New knot vectors in both u- and v-direction and control points
        (ubar, vbar, qw).

    *Reference:* Algorithm A5.5 from "The NURBS Book."
    """
    if d.lower() in ['u']:
        x = filter_knot_vect(n, p, uk, x)
        a = find_span(n, p, x[0], uk)
        b = find_span(n, p, x[-1], uk) + 1
        r = len(x) - 1
        nn = n + p + 1
        # Save unchanged control points.
        qw = zeros((r + p + b - a + 1, m + 1, 4), dtype=float64)
        for col in range(0, m + 1):
            for i in range(0, a - p + 1):
                qw[i, col] = cpw[i, col]
            for i in range(b - 1, n + 1):
                qw[i + r + 1, col] = cpw[i, col]
        # Load unchanged knot vector.
        ubar = zeros(nn + r + 2, dtype=float64)
        vbar = array(vk)
        for j in range(0, a + 1):
            ubar[j] = uk[j]
        for j in range(b + p, nn + 1):
            ubar[j + r + 1] = uk[j]
        # Find new control points.
        i = b + p - 1
        k = b + p + r
        for j in range(r, -1, -1):
            while x[j] <= uk[i] and i > a:
                ubar[k] = uk[i]
                for col in range(0, m + 1):
                    qw[k - p - 1, col] = cpw[i - p - 1, col]
                k -= 1
                i -= 1
            for col in range(0, m + 1):
                qw[k - p - 1, col] = qw[k - p, col]
            for t in range(1, p + 1):
                ind = k - p + t
                alpha = ubar[k + t] - x[j]
                if abs(alpha) == 0.0:
                    for col in range(0, m + 1):
                        qw[ind - 1, col] = qw[ind, col]
                else:
                    alpha /= (ubar[k + t] - uk[i - p + t])
                    for col in range(0, m + 1):
                        qw[ind - 1, col] = (alpha * qw[ind - 1, col] +
                                            (1.0 - alpha) * qw[ind, col])
            ubar[k] = x[j]
            k -= 1
        return ubar, vbar, qw
    if d.lower() in ['v']:
        x = filter_knot_vect(m, q, vk, x)
        a = find_span(m, q, x[0], vk)
        b = find_span(m, q, x[-1], vk) + 1
        r = len(x) - 1
        mm = m + q + 1
        # Save unchanged control points.
        qw = zeros((n + 1, r + q + b - a + 1, 4), dtype=float64)
        for row in range(0, n + 1):
            for i in range(0, a - q + 1):
                qw[row, i] = cpw[row, i]
            for i in range(b - 1, m + 1):
                qw[row, i + r + 1] = cpw[row, i]
        # Load unchanged knot vector.
        vbar = zeros(mm + r + 2, dtype=float64)
        ubar = array(uk)
        for j in range(0, a + 1):
            vbar[j] = vk[j]
        for j in range(b + p, mm + 1):
            vbar[j + r + 1] = vk[j]
        # Find new control points.
        i = b + q - 1
        k = b + q + r
        for j in range(r, -1, -1):
            while x[j] <= vk[i] and i > a:
                vbar[k] = vk[i]
                for row in range(0, n + 1):
                    qw[row, k - q - 1] = cpw[row, i - q - 1]
                k -= 1
                i -= 1
            for row in range(0, n + 1):
                qw[row, k - q - 1] = qw[row, k - q]
            for t in range(1, q + 1):
                ind = k - q + t
                alpha = vbar[k + t] - x[j]
                if abs(alpha) == 0.0:
                    for row in range(0, n + 1):
                        qw[row, ind - 1] = qw[row, ind]
                else:
                    alpha /= (vbar[k + t] - vk[i - q + t])
                    for row in range(0, n + 1):
                        qw[row, ind - 1] = (alpha * qw[row, ind - 1] +
                                            (1.0 - alpha) * qw[row, ind])
            vbar[k] = x[j]
            k -= 1
        return ubar, vbar, qw


def decompose_curve(n, p, uk, cpw):
    """
    Decompose NURBS curve into Bezier segments.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.

    :return: Number and Bezier segments in qw[j, k],
        where *k* is the control point of the *j*-th segment (nb, qw).
    :rtype: tuple

    *Reference:* Algorithm A5.6 from "The NURBS Book."
    """
    m = n + p + 1
    a = p
    b = p + 1
    nb = 0
    qw = zeros((m + 1, p + 1, 4), dtype=float64)
    # Track domain of each decomposed curve with ab.
    ab = zeros((m + 1, 2), dtype=float64)
    for i in range(0, p + 1):
        qw[nb, i] = cpw[i]
    while b < m:
        i = b
        ab[nb, 0] = uk[b - 1]
        # while b < m and uk[b + 1] == uk[b]:
        while b < m and CmpFlt.eq(uk[b + 1], uk[b], Settings.ptol):
            b += 1
        ab[nb, 1] = uk[b]
        mult = b - i + 1
        if mult < p:
            numer = uk[b] - uk[a]
            alphas = zeros(p + 1, dtype=float64)
            for j in range(p, mult, -1):
                alphas[j - mult - 1] = numer / (uk[a + j] - uk[a])
            for j in range(1, p - mult + 1):
                save = p - mult - j
                s = mult + j
                for k in range(p, s - 1, -1):
                    alpha = alphas[k - s]
                    qw[nb, k] = (alpha * qw[nb, k] +
                                 (1.0 - alpha) * qw[nb, k - 1])
                if b < m:
                    qw[nb + 1, save] = qw[nb, p]
        nb += 1
        if b < m:
            for i in range(p - mult, p + 1):
                qw[nb, i] = cpw[b - p + i]
            a = b
            b += 1
    return nb, qw, ab


def decompose_surface(n, p, uk, m, q, vk, cpw, d):
    """
    Decompose NURBS surface into Bezier patches.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param str d: Direction to decompose surface ("u": u-direction,
        "v": v-direction.

    :return: Computes a Bezier strip, i.e., a NURBS surface that is Bezier in
        one direction and B-spline in the other. The return must be called
        twice, once in the u-direction to get the Bezier strips, and then each
        strip must be fed into the routine in the v-direction to get the Bezier
        patches. The routine returns (nb, qw, ab) where *nb* is the number of
        Bezier strips, *qw* is the control points for each strip, and *ab* is
        the domain of each patch in the specified direction.
    :rtype: tuple

    *Reference:* Algorithm A5.7 from "The NURBS Book."
    """
    if d.lower() in ['u']:
        a = p
        b = p + 1
        nb = 0
        r = n + p + 1
        qw = zeros((r + 1, p + 1, m + 1, 4), dtype=float64)
        # Track domain of each decomposed surface with ab.
        ab = zeros((r + 1, 2), dtype=float64)
        for i in range(0, p + 1):
            for col in range(0, m + 1):
                qw[nb, i, col] = cpw[i, col]
        while b < r:
            i = b
            ab[nb, 0] = uk[b - 1]
            # while b < r and uk[b + 1] == uk[b]:
            while b < r and CmpFlt.eq(uk[b + 1], uk[b], Settings.ptol):
                b += 1
            ab[nb, 1] = uk[b]
            mult = b - i + 1
            if mult < p:
                numer = uk[b] - uk[a]
                alphas = zeros(p + 1, dtype=float64)
                for j in range(p, mult, -1):
                    alphas[j - mult - 1] = numer / (uk[a + j] - uk[a])
                for j in range(1, p - mult + 1):
                    save = p - mult - j
                    s = mult + j
                    for k in range(p, s - 1, -1):
                        alpha = alphas[k - s]
                        for col in range(0, m + 1):
                            qw[nb, k, col] = (alpha * qw[nb, k, col] +
                                              (1.0 - alpha) *
                                              qw[nb, k - 1, col])
                    if b < r:
                        for col in range(0, m + 1):
                            qw[nb + 1, save, col] = qw[nb, p, col]
            nb += 1
            if b < r:
                for i in range(p - mult, p + 1):
                    for col in range(0, m + 1):
                        qw[nb, i, col] = cpw[b - p + i, col]
                a = b
                b += 1
        return nb, qw, ab
    if d.lower() in ['v']:
        a = q
        b = q + 1
        nb = 0
        # Note r and s are not switched for u/v directions. This was only done
        # to allow for minor changes for v-direction code and since s is
        # already used.
        r = m + q + 1
        qw = zeros((r + 1, n + 1, q + 1, 4), dtype=float64)
        # Track domain of each decomposed surface with ab.
        ab = zeros((r + 1, 2), dtype=float64)
        for i in range(0, q + 1):
            for row in range(0, n + 1):
                qw[nb, row, i] = cpw[row, i]
        while b < r:
            i = b
            ab[nb, 0] = vk[b - 1]
            # while b < r and vk[b + 1] == vk[b]:
            while b < r and CmpFlt.eq(vk[b + 1], vk[b], Settings.ptol):
                b += 1
            ab[nb, 1] = vk[b]
            mult = b - i + 1
            if mult < q:
                numer = vk[b] - vk[a]
                alphas = zeros(q + 1, dtype=float64)
                for j in range(q, mult, -1):
                    alphas[j - mult - 1] = numer / (vk[a + j] - vk[a])
                for j in range(1, q - mult + 1):
                    save = q - mult - j
                    s = mult + j
                    for k in range(q, s - 1, -1):
                        alpha = alphas[k - s]
                        for row in range(0, n + 1):
                            qw[nb, row, k] = (alpha * qw[nb, row, k] +
                                              (1.0 - alpha) *
                                              qw[nb, row, k - 1])
                    if b < r:
                        for row in range(0, n + 1):
                            qw[nb + 1, row, save] = qw[nb, row, q]
            nb += 1
            if b < r:
                for i in range(q - mult, q + 1):
                    for row in range(0, n + 1):
                        qw[nb, row, i] = cpw[row, b - q + i]
                a = b
                b += 1
        return nb, qw, ab


# Implement curve knot removal when needed (Algorithm A5.8).
# For advanced curve/surface knot removal algorithms see:
# Tiller, W., Knot-removal algorithms for NURBS curves and surfaces,
# CAD, Vol. 24, No. 8, pp. 445-453, 1992.

# Implement curve/surface degree elevation when needed (A5.9 and A5.10).


def elevate_nurbs_curve_degree(n, p, uk, cpw, t=1):
    """
    Elevate the degree of a NURBS curve *t* times.

    :param int n:
    :param int p:
    :param ndarray uk:
    :param ndarray cpw:
    :param int t:

    :return: New number of control points - 1, new knot vector, and the new
        control points of NURBS curve with elevated degree (nq, uq, qw).
    :rtype: tuple

    *Reference:* Algorithm A5.9 from "The NURBS Book."
    """
    m = n + p + 1
    ph = p + t
    ph2 = int(floor(ph / 2))
    bezalfs = zeros((p + t + 1, p + 1), dtype=float64)
    bezalfs[0, 0], bezalfs[ph, p] = 1.0, 1.0
    for i in range(1, ph2 + 1):
        inv = 1.0 / bin_coeff(ph, i)
        mpi = min(p, i)
        for j in range(max(0, i - t), mpi + 1):
            bezalfs[i, j] = inv * bin_coeff(p, j) * bin_coeff(t, i - j)
    for i in range(ph2 + 1, ph):
        mpi = min(p, i)
        for j in range(max(0, i - t), mpi + 1):
            bezalfs[i, j] = bezalfs[ph - i, p - j]
    mh = ph
    kind = ph + 1
    r = -1
    a = p
    b = p + 1
    cind = 1
    ua = uk[0]
    _, mult, _ = find_mult_knots(n, p, uk)
    s = mh
    for mi in mult[1:-1]:
        s += mi + t
    qw = zeros((s + 1, 4), dtype=float64)
    uq = zeros(s + ph + 2, dtype=float64)
    bpts = zeros((p + 1, 4), dtype=float64)
    nextbpts = zeros((p - 1, 4), dtype=float64)
    ebpts = zeros((p + t + 1, 4), dtype=float64)
    qw[0] = cpw[0]
    for i in range(0, ph + 1):
        uq[i] = ua
    for i in range(0, p + 1):
        bpts[i] = cpw[i]
    while b < m:
        i = b
        # while b < m and uk[b] == uk[b + 1]:
        while b < m and CmpFlt.eq(uk[b], uk[b + 1], Settings.ptol):
            b += 1
        mul = b - i + 1
        mh += mul + t
        ub = uk[b]
        oldr = r
        r = p - mul
        # Insert knot u[b] r times.
        if oldr > 0:
            lbz = int(floor((oldr + 2) / 2))
        else:
            lbz = 1
        if r > 0:
            rbz = ph - int(floor((r + 1) / 2))
        else:
            rbz = ph
        if r > 0:
            numer = ub - ua
            alfs = zeros(p - 1, dtype=float64)
            for k in range(p, mul, -1):
                alfs[k - mul - 1] = numer / (uk[a + k] - ua)
            for j in range(1, r + 1):
                sav = r - j
                s = mul + j
                for k in range(p, s - 1, -1):
                    bpts[k] = (alfs[k - s] * bpts[k] +
                               (1.0 - alfs[k - s]) * bpts[k - 1])
                nextbpts[sav] = bpts[p]
        for i in range(lbz, ph + 1):
            ebpts[i, :] = 0.0
            mpi = min(p, i)
            for j in range(max(0, i - t), mpi + 1):
                ebpts[i] += bezalfs[i, j] * bpts[j]
        if oldr > 1:
            first = kind - 2
            last = kind
            den = ub - ua
            bet = (ub - uq[kind - 1]) / den
            for tr in range(1, oldr):
                i = first
                j = last
                kj = j - kind + 1
                while j - i > tr:
                    if i < cind:
                        alf = (ub - uq[i]) / (ua - uq[i])
                        qw[i] = alf * qw[i] + (1.0 - alf) * qw[i - 1]
                    if j >= lbz:
                        if j - tr <= kind - ph + oldr:
                            gam = (ub - uq[j - tr]) / den
                            ebpts[kj] = (gam * ebpts[kj] +
                                         (1.0 - gam) * ebpts[kj + 1])
                        else:
                            ebpts[kj] = (bet * ebpts[kj] +
                                         (1.0 - bet) * ebpts[kj + 1])
                        i += 1
                        j -= 1
                        kj -= 1
                    first -= 1
                    last += 1
        if a != p:
            for i in range(0, ph - oldr):
                uq[kind] = ua
                kind += 1
        for j in range(lbz, rbz + 1):
            qw[cind] = ebpts[j]
            cind += 1
        if b < m:
            for j in range(0, r):
                bpts[j] = nextbpts[j]
            for j in range(r, p + 1):
                bpts[j] = cpw[b - p + j]
            a = b
            b += 1
            ua = ub
        else:
            for i in range(0, ph + 1):
                uq[kind + i] = ub
    nq = mh - ph - 1
    return nq, uq, qw


def reverse_nurbs_curve(n, p, uk, cpw):
    """
    Reverse the direction of a NUBRS curve without changing its
    parameterization.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.

    :return: New knot vector and control points of curve (uq, qw).
    :rtype: tuple
    """
    m = n + p + 1
    uq = zeros(m + 1, dtype=float64)
    qw = zeros((n + 1, 4), dtype=float64)
    # Knot vector.
    a = uk[p]
    b = uk[n + 1]
    uq[:p + 1] = a
    uq[n + 1:] = b
    for i in range(1, m - 2 * p):
        uq[m - p - i] = -uk[p + i] + a + b
    # Control points.
    for i in range(0, n + 1):
        qw[i] = cpw[n - i]
    return uq, qw


def move_nurbs_surface_seam(n, p, uk, m, q, vk, cpw, uv, d):
    """
    Move the seam of a closed surface to a new parameter value.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param float uv: New seam parameter.
    :param str d: Direction to move seam in (surface must be closed in that
        direction).

    :return: Knot vector and control points of new surface.
    :rtype: tuple
    """
    ptol = Settings.ptol
    if d.lower() in ['u']:
        if CmpFlt.le(uv, uk[p], ptol) or CmpFlt.ge(uv, uk[n + 1], ptol):
            return uk, vk, cpw

        # Check multiplicity of seam parameter.
        _, s = find_span_mult(n, p, uv, uk)

        # Insert seam parameter s - p times.
        uq, vq, qw = surface_knot_ins(n, p, uk, m, q, vk, cpw, 'u', uv, p - s)

        # Find span of seam parameter in new knot vector.
        nw = qw.shape[0] - 1
        k = find_span(nw, p, uv, uq)

        # Get index for first row of control points.
        n0 = k - p

        # Generate control points for new seam.
        qw_seam = array(qw)
        indx = 0
        # For new seam to end of old seam
        for i in range(n0, nw + 1):
            qw_seam[indx, :, :] = qw[i, :, :]
            indx += 1
        # From start of old seam to new seam.
        for i in range(1, n0):
            qw_seam[indx, :, :] = qw[i, :, :]
            indx += 1
        # Set last row.
        qw_seam[-1, :, :] = qw_seam[0, :, :]

        # Generate new knot vector.
        uq_seam = array(uq)
        # From seam to end
        indx = p + 1
        for i in range(k + 1, nw + p + 1):
            uq_seam[indx] = uq[i] - uv
            indx += 1
        # From start to seam
        for i in range(p + 1, k + 1 - p):
            uq_seam[indx] = uq_seam[indx - 1] + (uq[i] - uq[i - 1])
            indx += 1
        uq_seam[-p - 1:] = uq[-1]
        return uq_seam, vk, qw_seam

    if d.lower() in ['v']:
        if CmpFlt.le(uv, vk[q], ptol) or CmpFlt.ge(uv, vk[m + 1], ptol):
            return uk, vk, cpw

        # Check multiplicity of seam parameter.
        _, s = find_span_mult(m, q, uv, vk)

        # Insert seam parameter s - q times.
        uq, vq, qw = surface_knot_ins(n, p, uk, m, q, vk, cpw, 'v', uv, q - s)

        # Find span of seam parameter in new knot vector.
        mw = qw.shape[1] - 1
        k = find_span(mw, q, uv, vq)

        # Get index for first row of control points.
        m0 = k - q

        # Generate control points for new seam.
        qw_seam = array(qw)
        indx = 0
        # For new seam to end of old seam
        for i in range(m0, mw + 1):
            qw_seam[:, indx, :] = qw[:, i, :]
            indx += 1
        # From start of old seam to new seam.
        for i in range(1, m0):
            qw_seam[:, indx, :] = qw[:, i, :]
            indx += 1
        # Set last row.
        qw_seam[:, -1, :] = qw_seam[:, 0, :]

        # Generate new knot vector.
        vq_seam = array(vq)
        # From seam to end
        indx = q + 1
        for i in range(k + 1, mw + q + 1):
            vq_seam[indx] = vq[i] - uv
            indx += 1
        # From start to seam
        for i in range(q + 1, k + 1 - q):
            vq_seam[indx] = vq_seam[indx - 1] + (vq[i] - vq[i - 1])
            indx += 1
        vq_seam[-q - 1:] = vq[-1]
        return uk, vq_seam, qw_seam
