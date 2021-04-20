from typing import Dict, Tuple, List


def read_atom_set() -> Dict[Tuple, int]:
    atom_set: Dict[Tuple, int] = {}
    with open("./Test_modules/atm_set_example.txt", "r") as f:
        for line in f:
            line_s = line.split("\t")
            indexes = line_s[0].split(" ")
            index = (int(indexes[0]), int(indexes[1]), int(indexes[2]))
            atom_set[index] = int(line_s[1])
    return atom_set


def read_bonds() -> Dict[Tuple, List]:
    bonds: Dict[Tuple, List] = {}
    with open("./Test_modules/bonds_example.txt", "r") as f:
        for line in f:
            line_s = line.split("\t")
            indexes = line_s[0].split(" ")
            index = (int(indexes[0]), int(indexes[1]), int(indexes[2]))
            values = []
            for i in range(1, len(line_s)):
                bond_indexes = line_s[i].split(" ")
                if bond_indexes == ["\n"]:
                    pass
                else:
                    bond_index = (
                        int(bond_indexes[0]),
                        int(bond_indexes[1]),
                        int(bond_indexes[2]),
                    )
                    values.append(bond_index)

            bonds[index] = values
    return bonds


def read_lattice() -> Dict:
    lattice: Dict[Tuple, List] = {}
    with open("./Test_modules/lattice_example.txt", "r") as f:
        for line in f:
            line_s = line.split("\t")
            indexes = line_s[0].split(" ")
            index = (int(indexes[0]), int(indexes[1]), int(indexes[2]))
            cood = line_s[1].split(" ")
            value = [float(cood[0]), float(cood[1]), float(cood[2])]
            lattice[index] = value
    return lattice


def read_candidate() -> Dict:
    candidate: Dict[Tuple, List] = {}
    with open("./Test_modules/candidate_example.txt", "r") as f:
        for line in f:
            line_s = line.split("\t")
            indexes = line_s[0].split(" ")
            index = (int(indexes[0]), int(indexes[1]), int(indexes[2]))
            values = []

            for i in range(1, len(line_s)):
                cand_indexes = line_s[i].split(" ")
                if cand_indexes == ["\n"]:
                    pass
                else:
                    cand_index = (
                        int(cand_indexes[0]),
                        int(cand_indexes[1]),
                        int(cand_indexes[2]),
                    )
                    values.append(cand_index)

            candidate[index] = values
    return candidate


if __name__ == "__main__":
    atom_set = read_atom_set()
    # print(atom_set[(0, 1, 0)])
    bonds = read_bonds()
    # print(bonds[(0, 0, 10)])
    lattice = read_lattice()
    # print(lattice[(1, 1, 7)])
