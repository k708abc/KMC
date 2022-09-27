# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from libc.math cimport exp

cdef double rate(double pre, double kbt, double energy):
    return pre * exp((-1)* energy / kbt)


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
