from numpy import identity, float64, array, mat


def _get_tm(dx=0., dy=0., dz=0.):
    """
    Generate a homogeneous translation matrix.
    """
    tm = identity(4, dtype=float64)
    tm[:-1, -1] = [dx, dy, dz]
    return mat(tm)


def translate_curve(c, v, d, inplace=False):
    """
    Translate a curve along the vector by a specified distance.

    :param c:
    :param v:
    :param d:
    :param inplace:

    :return:
    """
    v = array(v, dtype=float64)
    cpw = c.cpw
    tm = _get_tm(*v * d)
    cpw_new = array(cpw, dtype=float64)
    for i in range(0, cpw_new.shape[0]):
        cpw_new[i, :] = (tm * cpw_new[i, :].reshape(-1, 1)).flatten()
    if inplace:
        c.set_cpw(cpw_new)
        return True
    elif not inplace:
        cnew = c.copy()
        cnew.set_cpw(cpw_new)
        return cnew
    return False
