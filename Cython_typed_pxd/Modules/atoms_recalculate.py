from typing import List, Dict, Tuple
import cython


def recalculate(
    target: Tuple[int, int, int],
    atom_set: Dict[Tuple[int, int, int], int],
    bonds,
    diffuse_candidates: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]],
    height_change,
    trans_val,
):
    candidate = [fill for fill in diffuse_candidates[target] if atom_set[fill] == 1] + [
        target
    ]
    for bond in bonds[target]:
        if atom_set[bond] == 0:
            candidate += [
                fill for fill in diffuse_candidates[bond] if atom_set[fill] == 1
            ]

    if height_change is True:
        candidate += [
            (target[0], target[1], i)
            for i in range(trans_val, -1, -1)
            if atom_set[(target[0], target[1], i)] == 1
        ]
    return candidate
