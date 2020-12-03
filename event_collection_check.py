from typing import List, Dict, Tuple


def random_target(atom_set: Dict):
    candidate = []
    for index, state in atom_set.items():
        if state != 0:
            candidate.append(index)
    return candidate


def event_check_poscar(
    atom_set: Dict,
    events: List,
    lattice: Dict,
    unit_length: int,
    maxz: int,
    target: Tuple[int, int, int],
):
    xp: List[float] = [lattice[target][0] / unit_length]
    yp: List[float] = [lattice[target][1] / unit_length]
    zp: List[float] = [lattice[target][2] / maxz / 2.448]
    atom_i = 0
    eve_i = 0
    for index, atom_state in atom_set.items():
        if atom_state != 0:
            xp.append(lattice[index][0] / unit_length)
            yp.append(lattice[index][1] / unit_length)
            zp.append(lattice[index][2] / maxz / 2.448)
            atom_i += 1
    for cand in events:
        xp.append(lattice[cand][0] / unit_length)
        yp.append(lattice[cand][1] / unit_length)
        zp.append(lattice[cand][2] / maxz / 2.448)
        eve_i += 1
    #
    file_name = "Event_check/Check_events_" + str(target) + ".vasp"
    file_data = open(file_name, "w")
    file_data.write("check events" + "\n")
    file_data.write("10" + "\n")
    file_data.write(str(unit_length) + "\t" + "0" + "\t" + "0" + "\n")
    file_data.write(
        str(unit_length / 2) + "\t" + str(unit_length / 2 * 1.732) + "\t" + "0" + "\n"
    )
    file_data.write("0" + "\t" + "0" + "\t" + str(maxz * 2.448) + "\n")
    file_data.write("Bi" + "\t" + "Si" + "\t" + "O" + "\n")
    file_data.write(str(1) + "\t" + str(atom_i) + "\t" + str(eve_i) + "\n")
    file_data.write("direct" + "\n")
    for i in range(len(xp)):
        file_data.write(str(xp[i]) + "\t" + str(yp[i]) + "\t" + str(zp[i]) + "\n")
    file_data.close()
