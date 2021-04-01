from typing import List, Dict, Tuple


def dep_check_poscar(
    atom_set: Dict[Tuple[int, int, int], int],
    candidate: List[Tuple[int, int, int]],
    lattice: Dict,
    unit_length: int,
    maxz: int,
):
    xp: List[float] = []
    yp: List[float] = []
    zp: List[float] = []
    atom_i = 0
    cand_i = 0
    for index, atom_state in atom_set.items():
        if atom_state != 0:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_i += 1
    for cand in candidate:
        xp.append(lattice[cand][0] / unit_length)
        yp.append(lattice[cand][1] / unit_length)
        zp.append(lattice[cand][2] / maxz / 2.448)
        cand_i += 1
    #
    file_data = open(".\Test_modules\poscar_for_check_depsotion.vasp", "w")
    file_data.write("check depsotion" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Si" + "\t" + "O" + "\n")
    file_data.write(str(atom_i) + "\t" + str(cand_i) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()
