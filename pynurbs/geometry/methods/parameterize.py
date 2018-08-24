from __future__ import division

from numpy import zeros, float64, diff, sum, array
from numpy.linalg import norm

from pynurbs.geometry.point import Point


def uniform(pnts, a=0., b=1.):
    """
    Generate a uniform parameters.

    :param pnts: 1D set of ordered points.
    :etype pnts: Points or array_like
    :param float a: Beginning domain if not 0.
    :param float b: Ending domain if not 1.

    :return: Uniformly spaced parameters between [a, b].
    :rtype: ndarray
    """
    if isinstance(pnts[0], Point):
        n = len(pnts)
    else:
        n = pnts.shape[0]

    u = zeros(n, dtype=float64)
    u[0] = a
    u[-1] = b
    for i in range(1, n - 1):
        u[i] = a + i * (b - a) / n
    return u


def chord_length(pnts, a=0., b=1.):
    """
    Generate parameters using chord length method.

    :param pnts: List or array of ordered points.
    :etype pnts: Points or array_like
    :param float a: Beginning domain if not 0.
    :param float b: Ending domain if not 1.

    :return: Parameters between [a, b].
    :rtype: ndarray
    """
    pnts = array(pnts, dtype=float64)
    n = len(pnts)
    dtotal = sum(norm(diff(pnts, axis=0), axis=1))
    u = zeros(n, dtype=float64)
    u[0] = a
    u[-1] = b
    if dtotal <= 0.:
        return u
    for i in range(1, n - 1):
        di = norm(pnts[i] - pnts[i - 1]) / dtotal
        u[i] = u[i - 1] + di
    return u


def centripetal(pnts, a=0., b=1.):
    """
    Generate parameters using centripetal method.

    :param pnts: List or array of ordered points.
    :etype pnts: Points or array_like
    :param float a: Beginning domain if not 0.
    :param float b: Ending domain if not 1.

    :return: Parameters between [a, b].
    :rtype: ndarray
    """
    pnts = array(pnts, dtype=float64)
    n = len(pnts)
    dtotal = sum(norm(diff(pnts, axis=0), axis=1))
    u = zeros(n, dtype=float64)
    u[0] = a
    u[-1] = b
    if dtotal <= 0.:
        return u
    for i in range(1, n - 1):
        di = (norm(pnts[i] - pnts[i - 1]) / dtotal) ** 0.5
        u[i] = u[i - 1] + di
    return u
