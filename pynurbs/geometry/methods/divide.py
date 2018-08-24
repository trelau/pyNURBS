from __future__ import division

from numpy import zeros, float64, array

from pynurbs.config import Settings
from pynurbs.geometry.methods.compare_floats import CompareFloats as CmpFlt
from pynurbs.geometry.methods.evaluate import find_span_mult, deCasteljau1
from pynurbs.geometry.methods.modify import curve_knot_ins, surface_knot_ins


def split_bezier_curve(cpw, n, u, a=0., b=1.):
    """
    Split Bezier curve into two segments at *u*.

    :param int n: Number of control points - 1.
    :param ndarray cpw: Control points.
    :param float u: Parametric point (0 <= u <= 1).
    :param float a: Lower domain of curve (0 <= a < b).
    :param float b: Upper domain of curve (a < b <= 1).

    :return: New control points and domain for two curves
        (qw1, a1, b1, qw2, a2, b2).
    :rtype: tuple
    """
    # Split using deCasteljau algorithm, tracking first and last points in
    # each column.
    qw1 = zeros((n + 1, 4), dtype=float64)
    qw2 = array(cpw, copy=True)
    qw1[0] = qw2[0]
    for k in range(1, n + 1):
        for i in range(0, n - k + 1):
            qw2[i] = (1. - u) * qw2[i] + u * qw2[i + 1]
        qw1[k] = qw2[0]
    # Adjust domain.
    a1 = a
    b1 = a + u * (b - a)
    a2 = b1
    b2 = b
    return qw1, a1, b1, qw2, a2, b2


def extract_bezier_curve(cpw, n, u0, u1, a=0., b=1.):
    """
    Extract a Bezier curve between parameters *u0* and *u1*.

    :param ndarray cpw: Control points.
    :param int n: Number of control points - 1 (degree).
    :param float u0: Starting parameter (0 <= u0 < u1).
    :param float u1: Ending parameter (u0 < u1 <= 1).
    :param float a: Lower domain of curve (0 <= a < b).
    :param float b: Upper domain of curve (a < b <= 1).

    :return: Control points of Bezier curve between *u0* and *u1* and new
        domain (qw, a1, b1).
    :rtype: ndarray
    """
    # Split at u0 keeping the right part.
    _, _, _, qw2, a2, b2 = split_bezier_curve(cpw, n, u0, a, b)
    # Split at u1 keeping the left part. Adjust u1 so it considers new curve.
    u1b = (u1 - u0) / (1. - u0)
    qw1, a1, b1, _, _, _ = split_bezier_curve(qw2, n, u1b, a2, b2)
    return qw1, a1, b1


def split_nurbs_curve(n, p, uk, cpw, u):
    """
    Split NURBS curve into two segments at *u*.

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.
    :param float u: Parametric point to split curve at.

    :return: New knot vector and control points for two curves
        (uk1, qw1, uk2, qw2).
    :rtype: tuple
    """
    # Special case if a >= u >= b.
    # if u <= uk[p]:
    ptol = Settings.ptol
    if CmpFlt.le(u, uk[p], ptol):
        qw1 = zeros((n + 1, 4), dtype=float64)
        qw1[:] = cpw[0]
        uk1 = array(uk, dtype=float64)
        uk1[:] = uk[p]
        qw2 = array(cpw, dtype=float64)
        uk2 = array(uk, dtype=float64)
        return uk1, qw1, uk2, qw2
    # elif u >= uk[n + 1]:
    elif CmpFlt.ge(u, uk[n + 1], ptol):
        qw1 = array(cpw, dtype=float64)
        uk1 = array(uk, dtype=float64)
        qw2 = zeros((n + 1, 4), dtype=float64)
        qw2[:] = cpw[n]
        uk2 = array(uk, dtype=float64)
        uk2[:] = uk[n + 1]
        return uk1, qw1, uk2, qw2
    # Find multiplicity of knot.
    k, s = find_span_mult(n, p, u, uk)
    # Insert knot p - s times and update knot vector and control points.
    if s >= p:
        uq, qw = array(uk, dtype=float64), array(cpw, dtype=float64)
    else:
        uq, qw = curve_knot_ins(n, p, uk, cpw, u, p - s)
    # Control points and knot vectors.
    qw1 = qw[:k - s + 1]
    qw2 = qw[k - s:]
    m = n + p + 1
    uk1 = zeros((k + 1) + (p - s + 1), dtype=float64)
    uk2 = zeros((m - k) + p + 1, dtype=float64)
    uk1[:-1] = uq[:k + (p - s) + 1]
    uk1[-1] = u
    uk2[1:] = uq[k - s + 1:]
    uk2[0] = u
    return uk1, qw1, uk2, qw2


def extract_nurbs_curve(n, p, uk, cpw, u0, u1):
    """
    Extract a NURBS curve between parameters *u0* and *u1* where (u0 < u1).

    :param int n: Number of control points - 1.
    :param int p: Degree.
    :param ndarray uk: Knot vector.
    :param ndarray cpw: Control points.
    :param float u0: Starting parameter.
    :param float u1: Ending parameter.

    :return: Knot vector and control points of NURBS curve between *u0* and
        *u1* (uq, qw).
    :rtype: tuple
    """
    if u0 > u1:
        u0, u1 = u1, u0
    # Split at u0 keeping the right part.
    _, _, bk, bw = split_nurbs_curve(n, p, uk, cpw, u0)
    # Split at u1 keeping the left part.
    n = bw.shape[0] - 1
    uq, qw, _, _ = split_nurbs_curve(n, p, bk, bw, u1)
    return uq, qw


def split_bezier_surface(cpw, n, m, u=None, v=None,
                         au=0., bu=1., av=0., bv=1.):
    """
    Split Bezier surface.

    :param ndarray cpw: Control points.
    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param u: Location of *u* split (in v-direction).
    :type u: float or None
    :param float v: Location of *v* split (in u-direction).
    :type v: float or None
    :param float au: Lower domain in u-direction.
    :param float bu: Upper domain in u-direction.
    :param float av: Lower domain in v-direction.
    :param float bv: Upper domain in v-direction.

    :return: Control points of two Bezier surfaces after split and domain
        data after split (qw1, qw2, ab1, ab2).
    :rtype: tuple

    **Note**: The domain after split for each surface are contained in 2 x 2
    matrices in the form:

         Data         Index
        au | bu     0,0 | 0,1
        ------- ---> ---------
        av | bv     1,0 | 1,1

    """
    # Use deCasteljau algorithm.
    qw2 = array(cpw, copy=True)
    if u is not None:
        # Split at u in v-direction.
        qw1 = zeros((n + 1, m + 1, 4), dtype=float64)
        qw1[0, :] = qw2[0, :]
        for j in range(0, m + 1):
            for k in range(1, n + 1):
                for i in range(0, n - k + 1):
                    qw2[i, j] = (1. - u) * qw2[i, j] + u * qw2[i + 1, j]
                qw1[k, j] = qw2[0, j]
        # Adjust domain.
        au1 = au
        bu1 = au + u * (bu - au)
        av1, bv1 = av, bv
        au2, bu2 = bu1, bu
        av2, bv2 = av, bv
        ab1 = array([[au1, bu1], [av1, bv1]], dtype=float64)
        ab2 = array([[au2, bu2], [av2, bv2]], dtype=float64)
        return qw1, qw2, ab1, ab2
    elif v is not None:
        # Split at v in u-direction.
        qw1 = zeros((n + 1, m + 1, 4), dtype=float64)
        qw1[:, 0] = qw2[:, 0]
        for i in range(0, n + 1):
            for k in range(1, m + 1):
                for j in range(0, m - k + 1):
                    qw2[i, j] = (1. - v) * qw2[i, j] + v * qw2[i, j + 1]
                qw1[i, k] = qw2[i, 0]
        # Adjust domain.
        au1, bu1 = au, bu
        av1 = av
        bv1 = av + v * (bv - av)
        au2, bu2 = au, bu
        av2, bv2 = bv1, bv
        ab1 = array([[au1, bu1], [av1, bv1]], dtype=float64)
        ab2 = array([[au2, bu2], [av2, bv2]], dtype=float64)
        return qw1, qw2, ab1, ab2


def extract_bezier_isocurve(n, m, cpw, u=None, v=None):
    """
    Extract iso-curve from Bezier surface.

    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param ndarray cpw: Control points.
    :param float u: Parameter to extract curve at (in v-direction).
    :param float v: Parameter to extract curve at (in u-direction).

    :return: Control points of Bezier iso-curve.
    :rtype: ndarray
    """
    if u is not None:
        qw = zeros((m + 1, 4), dtype=float64)
        for j in range(0, m + 1):
            qw[j] = deCasteljau1(cpw[:, j], n, u)
    else:
        qw = zeros((n + 1, 4), dtype=float64)
        for i in range(0, n + 1):
            qw[i] = deCasteljau1(cpw[i, :], m, v)
    return qw


def split_nurbs_surface(n, p, uk, m, q, vk, cpw, u=None, v=None):
    """
    Split NURBS surface into two segments in specified direction.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param u: Location of *u* split (in v-direction).
    :type u: float or None
    :param v: Location of *v* split (in u-direction).
    :type v: float or None

    :return: New knot vectors and control points for two surfaces
        (uk1, vk1, qw1, uk2, vk2, qw2).
    :rtype: tuple
    """
    if u is not None:
        # Split at u in v-direction.
        vk12 = array(vk, dtype=float64)
        if u <= uk[p]:
            qw1 = zeros((n + 1, m + 1, 4), dtype=float64)
            for i in range(0, n + 1):
                qw1[i, :] = cpw[0, :]
            qw2 = array(cpw, dtype=float64)
            uk1 = array(uk, dtype=float64)
            uk1[:] = uk[p]
            uk2 = array(uk, dtype=float64)
            return uk1, vk12, qw1, uk2, vk12, qw2
        if u >= uk[n + 1]:
            qw1 = array(cpw, dtype=float64)
            qw2 = zeros((n + 1, m + 1, 4), dtype=float64)
            for i in range(0, n + 1):
                qw2[i, :] = cpw[n, :]
            uk1 = array(uk, dtype=float64)
            uk2 = array(uk, dtype=float64)
            uk2[:] = uk[n + 1]
            return uk1, vk12, qw1, uk2, vk12, qw2
        # Find multiplicity of knot.
        k, s = find_span_mult(n, p, u, uk)
        # Insert knot p - s times and update knot vector and control points.
        if s >= p:
            uq, qw = array(uk, dtype=float64), array(cpw, dtype=float64)
        else:
            uq, _, qw = surface_knot_ins(n, p, uk, m, q, vk, cpw,
                                         'u', u, p - s)
        qw1 = qw[:k - s + 1, :]
        qw2 = qw[k - s:, :]
        r = n + p + 1
        uk1 = zeros((k + 1) + (p - s + 1), dtype=float64)
        uk2 = zeros((r - k) + p + 1, dtype=float64)
        uk1[:-1] = uq[:k + (p - s) + 1]
        uk1[-1] = u
        uk2[1:] = uq[k - s + 1:]
        uk2[0] = u
        return uk1, vk12, qw1, uk2, vk12, qw2
    elif v is not None:
        # Split at v in u-direction.
        uk12 = array(uk, dtype=float64)
        if v <= vk[q]:
            qw1 = zeros((n + 1, m + 1, 4), dtype=float64)
            for j in range(0, m + 1):
                qw1[:, j] = cpw[:, 0]
            qw2 = array(cpw, dtype=float64)
            vk1 = array(vk, dtype=float64)
            vk1[:] = vk[q]
            vk2 = array(vk, dtype=float64)
            return uk12, vk1, qw1, uk12, vk2, qw2
        if v >= vk[m + 1]:
            qw1 = array(cpw, dtype=float64)
            qw2 = zeros((n + 1, m + 1, 4), dtype=float64)
            for j in range(0, m + 1):
                qw2[:, j] = cpw[:, m]
            vk1 = array(vk, dtype=float64)
            vk2 = array(vk, dtype=float64)
            vk2[:] = vk[m + 1]
            return uk12, vk1, qw1, uk12, vk2, qw2
        # Find multiplicity of knot.
        k, s = find_span_mult(m, q, v, vk)
        # Insert knot q - s times and update knot vector and control points.
        if s >= q:
            vq, qw = array(vk, dtype=float64), array(cpw, dtype=float64)
        else:
            _, vq, qw = surface_knot_ins(n, p, uk, m, q, vk, cpw,
                                         'v', v, q - s)
        qw1 = qw[:, :k - s + 1]
        qw2 = qw[:, k - s:]
        r = m + q + 1
        vk1 = zeros((k + 1) + (q - s + 1), dtype=float64)
        vk2 = zeros((r - k) + q + 1, dtype=float64)
        vk1[:-1] = vq[:k + (q - s) + 1]
        vk1[-1] = v
        vk2[1:] = vq[k - s + 1:]
        vk2[0] = v
        return uk12, vk1, qw1, uk12, vk2, qw2


def extract_nurbs_isocurve(n, p, uk, m, q, vk, cpw, u=None, v=None):
    """
    Extract iso-curve from NURBS surface.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param float u: Location of *u* split (in v-direction).
    :param float v: Location of *v* split (in u-direction).

    :return: Knot vector and control points for NURBS iso-curve in specified
        direction (uq, qw).
    :rtype: tuple
    """
    if u is not None:
        # Split at u in v-direction.
        uq = array(vk, dtype=float64)
        # Special case if a >= u >= b.
        if u <= uk[p]:
            return uq, cpw[0, :]
        if u >= uk[n + 1]:
            return uq, cpw[n, :]
        # Determine multiplicity.
        k, s = find_span_mult(n, p, u, uk)
        qw = zeros((m + 1, 4), dtype=float64)
        for j in range(0, m + 1):
            _, qwi = curve_knot_ins(n, p, uk, cpw[:, j], u, p - s)
            qw[j] = qwi[k - s]
        return uq, qw
    elif v is not None:
        # Split at v in u-direction.
        uq = array(uk, dtype=float64)
        # Special case if a >= v >= b.
        if v <= vk[q]:
            return uq, cpw[:, 0]
        if v >= vk[m + 1]:
            return uq, cpw[:, m]
        # Determine multiplicity.
        k, s = find_span_mult(m, q, v, vk)
        qw = zeros((n + 1, 4), dtype=float64)
        for i in range(0, n + 1):
            _, qwi = curve_knot_ins(m, q, vk, cpw[i, :], v, q - s)
            qw[i] = qwi[k - s]
        return uq, qw


def extract_bezier_surface(n, m, cpw, u0, u1, v0, v1, au, bu, av, bv):
    """
    Extract Bezier surface bounded by parameters.

    :param int n: Number of control points - 1 in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param ndarray cpw: Control points.
    :param float u0: Starting parameter in u-direction.
    :param float u1: Ending parameter in u-direction.
    :param float v0: Starting parameter in v-direction.
    :param float v1: Ending parameter in v-direction.
    :param float au: Lower domain in u-direction.
    :param float bu: Upper domain in u-direction.
    :param float av: Lower domain in v-direction.
    :param float bv: Upper domain in v-direction.

    :return: Extracted Bezier surface.
    :rtype: :class:`.BezierSurface`
    """
    # Split at u0 keeping upper part.
    _, qw1, _, ab1 = split_bezier_surface(cpw, n, m, u0, None, au, bu, av, bv)
    # Split at u1 keeping lower part.
    u1b = (u1 - u0) / (1. - u0)
    qw2, _, ab2, _ = split_bezier_surface(qw1, n, m, u1b, None,
                                          ab1[0, 0], ab1[0, 1],
                                          ab1[1, 0], ab1[1, 1])
    # Split at v0 keeping upper part.
    _, qw3, _, ab3 = split_bezier_surface(qw2, n, m, None, v0,
                                          ab2[0, 0], ab2[0, 1],
                                          ab2[1, 0], ab2[1, 1])
    # Split at v1 keeping lower part.
    v1b = (v1 - v0) / (1. - v0)
    qw4, _, ab4, _ = split_bezier_surface(qw3, n, m, None, v1b,
                                          ab3[0, 0], ab3[0, 1],
                                          ab3[1, 0], ab3[1, 1])
    # Return control points and domain.
    return qw4, ab4


def extract_nurbs_surface(n, p, uk, m, q, vk, cpw, u0, u1, v0, v1):
    """
    Extract NURBS surface bounded by parameters.

    :param int n: Number of control points - 1 in u-direction.
    :param int p: Degree in u-direction.
    :param ndarray uk: Knot vector in u-direction.
    :param int m: Number of control points - 1 in v-direction.
    :param int q: Degree in v-direction.
    :param ndarray vk: Knot vector in v-direction.
    :param ndarray cpw: Control points.
    :param float u0: Starting parameter in u-direction.
    :param float u1: Ending parameter in u-direction.
    :param float v0: Starting parameter in v-direction.
    :param float v1: Ending parameter in v-direction.

    :return: Extracted NURBS surface (nq, mq, uq, vq, qw).
    :rtype: :class:`.NurbsSurface`
    """
    # Split at u0 keeping upper part.
    _, _, _, uk1, vk1, qw1 = split_nurbs_surface(n, p, uk, m, q, vk, cpw,
                                                 u0, None)
    # Split at u1 keeping lower part.
    nm = qw1.shape[:-1]
    n1, m1 = nm[0] - 1, nm[1] - 1
    uk2, vk2, qw2, _, _, _ = split_nurbs_surface(n1, p, uk1, m1, q, vk1, qw1,
                                                 u1, None)
    # Split at v0 keeping upper part.
    nm = qw2.shape[:-1]
    n2, m2 = nm[0] - 1, nm[1] - 1
    _, _, _, uk3, vk3, qw3 = split_nurbs_surface(n2, p, uk2, m2, q, vk2, qw2,
                                                 None, v0)
    # Split at v1 keeping lower part.
    nm = qw3.shape[:-1]
    n3, m3 = nm[0] - 1, nm[1] - 1
    uk4, vk4, qw4, _, _, _ = split_nurbs_surface(n3, p, uk3, m3, q, vk3, qw3,
                                                 None, v1)
    nm = qw4.shape[:-1]
    n4, m4 = nm[0] - 1, nm[1] - 1
    return n4, m4, uk4, vk4, qw4
