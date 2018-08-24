from numpy import array, float64, inf
from numpy.linalg import norm


def si_curve_nearest_point(si, pref, tol=None):
    """
    Find the index of an intersection curve nearest the given point.
    """
    pref = array(pref, dtype=float64)

    if not si.success:
        return None
    if si.ncrvs == 1:
        return 0

    # For each curve, calculate the distance to the segments between
    # intersections points and find the minimum distance and index.
    dmin = inf
    crv_indx = None
    for i in range(0, si.ncrvs):
        pnts = array(si.icurves[i].cp, dtype=float64)
        n = pnts.shape[0]
        for j in range(0, n):
            di = norm(pref - pnts[j])
            if di < dmin:
                dmin = di
                crv_indx = i

    # Check tolerance and return.
    if tol is None:
        return crv_indx

    if dmin <= tol:
        return crv_indx
    return None
