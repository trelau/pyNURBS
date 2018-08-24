from numpy import float64, linspace, zeros, hstack, unique
from numpy.linalg import norm
from scipy.interpolate import LinearNDInterpolator, RectBivariateSpline


def surface_distance_map(surface, n=None, h=None):
    """
    Build a distance map of the
    """
    if n is None and h is None:
        return None, None, None, None

    # Use length of bounding box diagonal to estimate grid steps.
    if h is not None:
        bbox = surface.bbox
        diag = bbox.diagonal_length()
        ntemp = int(diag / h) + 1
        if ntemp > n:
            n = ntemp

    # Parameter domain.
    umin = surface.au
    umax = surface.bu
    vmin = surface.av
    vmax = surface.bv

    # Get knots of surfaces.
    uknots, vknots = surface.uknots, surface.vknots

    # Grid of uv values.
    ugrid = linspace(umin, umax, n)
    vgrid = linspace(vmin, vmax, n)

    # Combine and make unique.
    ugrid = hstack((ugrid, uknots))
    vgrid = hstack((vgrid, vknots))
    ugrid = unique(ugrid)
    vgrid = unique(vgrid)

    # Number of grids.
    nu = ugrid.size
    nv = vgrid.size

    # Surface evaluator.
    seval = surface.eval

    # Map u-direction.
    udist = zeros((nu, nv), dtype=float64)
    for j in range(0, nv):
        dv = vgrid[j]
        d = 0.
        for i in range(1, nu):
            du_0 = ugrid[i - 1]
            du_1 = ugrid[i]
            pi = seval(du_0, dv, domain='global', rtype='ndarray')
            pi1 = seval(du_1, dv, domain='global', rtype='ndarray')
            d += norm(pi1 - pi)
            udist[i, j] = d

    # Map v-direction.
    vdist = zeros((nu, nv), dtype=float64)
    for i in range(0, nu):
        du = ugrid[i]
        d = 0.
        for j in range(1, nv):
            dv_0 = vgrid[j - 1]
            dv_1 = vgrid[j]
            pi = seval(du, dv_0, domain='global', rtype='ndarray')
            pi1 = seval(du, dv_1, domain='global', rtype='ndarray')
            d += norm(pi1 - pi)
            vdist[i, j] = d

    # Build interpolation functions for parameters.
    xy = []
    uparam = []
    vparam = []
    for j in range(0, nv):
        for i in range(0, nu):
            xy.append([udist[i, j], vdist[i, j]])
            uparam.append(ugrid[i])
            vparam.append(vgrid[j])

    # Build interpolation functions for distance.
    interp_udist = RectBivariateSpline(ugrid, vgrid, udist, kx=1, ky=1)
    interp_vdist = RectBivariateSpline(ugrid, vgrid, vdist, kx=1, ky=1)
    interp_uparam = LinearNDInterpolator(xy, uparam, rescale=True)
    interp_vparam = LinearNDInterpolator(xy, vparam, rescale=True)

    return interp_udist, interp_vdist, interp_uparam, interp_vparam
