# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False


import math

cdef double rate(double pre, kbt, energy):
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
    return pre * math.exp((-1) * energy / kbt)


"""
if __name__ == "__main__":
    from InputParameter import Params

    init_values = Params()
    pre = float(init_values.prefactor)
    kbt = init_values.temperature_eV
    print("Test : caluculate rate constant")
    energy = float(input("Input energy: "))
    print("prefactor: " + str(pre))
    print("kbt: " + str(kbt))
    print("energy: " + str(energy))
    print("rate = " + str(rate(pre, kbt, energy)))
"""
