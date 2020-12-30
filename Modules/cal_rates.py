import math
from InputParameter import Params


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


if __name__ == "__main__":
    init_values = Params()
    pre = float(init_values.prefactor)
    kbt = init_values.temperature_eV
    print("Test : caluculate rate constant")
    energy = float(input("Input energy: "))
    print("rate = " + str(rate(pre, kbt, energy)))
