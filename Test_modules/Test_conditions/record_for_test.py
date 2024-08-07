from typing import Dict


def rec_for_test(
    atom_set: Dict, bonds: Dict, lattice: Dict, candidate: Dict, path
) -> None:
    file_data = open(path + "atm_set_example.txt", "w")
    for atom_index, state in atom_set.items():
        for component_index in atom_index:
            file_data.write(str(component_index) + " ")
        file_data.write("\t" + str(state) + "\n")
    file_data.close()
    #
    file_data = open(path + "bonds_example.txt", "w")
    for atom_index, bond_all in bonds.items():
        for component_index in atom_index:
            file_data.write(str(component_index) + " ")
        file_data.write("\t")
        for bond in bond_all:
            for bond_components in bond:
                file_data.write(str(bond_components) + " ")
            file_data.write("\t")
        file_data.write("\n")
    file_data.close()
    #
    file_data = open(path + "lattice_example.txt", "w")
    for atom_index, position in lattice.items():
        for component in atom_index:
            file_data.write(str(component) + " ")
        file_data.write("\t")
        for pos in position:
            file_data.write(str(pos) + " ")
        file_data.write("\n")
    file_data.close()
    #
    file_data = open(path + "candidate_example.txt", "w")
    for atom_index, cands in candidate.items():
        for component in atom_index:
            file_data.write(str(component) + " ")
        file_data.write("\t")
        for cand in cands:
            for cand_component in cand:
                file_data.write(str(cand_component) + " ")
            file_data.write("\t")
        file_data.write("\n")
    file_data.close()
