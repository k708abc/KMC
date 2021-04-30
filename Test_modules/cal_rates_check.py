import sys

sys.path.append("../../KMC/")

from Modules.InputParameter import Params
from Modules.cal_rates import rate

if __name__ == "__main__":
    print("Chek cal_rates.py")
    init_values = Params("../kmc_input.yml")
    pre = float(init_values.prefactor)
    kbt = init_values.temperature_eV
    energy = float(input("Input energy: "))
    print("prefactor: " + str(pre))
    print("kbt: " + str(kbt))
    print("energy: " + str(energy))
    print("rate = " + str(rate(pre, kbt, energy)))