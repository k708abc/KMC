import math


def rate(pre: float, kbt: float, energy: float) -> float:
    """Return reaction the rate.

    Parameters
    ----------
    pre : float
        Prefactor
    kbt : float
        temperature in eV unit. (k_B T)
    energy : float
        Energy in eV unit

    Returns
    -------
    float
        reaction rate
    """
    return pre * math.exp(energy / kbt)
