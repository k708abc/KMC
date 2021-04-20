from InputParameter import Params
from Modules.lattice_form import lattice_form
import os
from typing import List


def form_poscar(target, lattice, candidate, unit_length, maxz):
    xp: List[float] = [lattice[target][0] / unit_length]
    yp: List[float] = [lattice[target][1] / unit_length]
    zp: List[float] = [lattice[target][2] / maxz / 2.448]
    atom_cand = 0
    atom_lattice = 0
    for index in candidate[target]:
        if index != target:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_cand += 1
    for index in lattice:
        if (index not in candidate[target]) and (index != target):
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_lattice += 1

    file_name = "Test_modules/Possible_diffusion/pos_" + str(target) + ".vasp"
    file_data = open(file_name, "w")
    file_data.write("check possible diffusion" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Bi" + "\t" + "Si" + "\t" + "O" + "\n")

    file_data.write(str(1) + "\t" + str(atom_cand) + "\t" + str(atom_lattice) + "\n")

    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()


if __name__ == "__main__":
    init_value = Params()

    dir_name = "Test_modules/Possible_diffusion/"
    os.makedirs(dir_name, exist_ok=True)

    (
        lattice,
        bonds,
        atom_set,
        event,
        event_time,
        event_time_tot,
        _,
        candidate,
    ) = lattice_form(init_value)

    unit_length = init_value.n_cell_init
    maxz = init_value.z_unit_init
    for key in lattice:
        if key[2] > maxz:
            pass
        else:
            form_poscar(key, lattice, candidate, unit_length, maxz)
    print("diffusion check file is fored")