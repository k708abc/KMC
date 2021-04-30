from typing import List, Tuple, Dict
import random
import math


# 20210401 checked


def find_candidates(atom_set: Dict, bonds: Dict) -> List[Tuple[int, int, int]]:
    candidate: List[Tuple[int, int, int]] = []
    for atom_index, condition in atom_set.items():
        if condition == 1:
            pass
        elif atom_index[2] == 0:
            candidate.append(atom_index)
        else:
            for bond in bonds[atom_index]:
                if atom_set[bond] == 1:
                    candidate.append(atom_index)
                    break
    return list(set(candidate))


def dep_position(candidate: List) -> Tuple[int, int, int]:
    return candidate[math.floor(random.random() * len(candidate))]


def remove_first(candidate) -> List[Tuple]:
    return [site for site in candidate if site[2] not in (0, 1)]


def deposit_an_atom(atom_set: Dict, bonds: Dict, params, empty_first: int) -> Tuple:
    candidate = find_candidates(atom_set, bonds)
    if (params.keep_defect_check is True) and (empty_first == int(params.num_defect)):
        candidate = remove_first(candidate)
    return dep_position(candidate)
