# from InputParameter import Params
from typing import List, Dict, Tuple


def recalculate(
    target: Tuple[int, int, int],
    atom_set: Dict[Tuple[int, int, int], int],
    diffuse_candidates: Dict[Tuple[int, int, int], List[Tuple[int, int, int]]],
    height_change,
    trans_val,
):
    candidate = [fill for fill in diffuse_candidates[target] if atom_set[fill] == 1] + [
        target
    ]
    if (height_change[0] < trans_val) and (height_change[1] >= trans_val):
        candidate += [
            (target[0], target[1], i)
            for i in range(target[0], -1, -1)
            if atom_set[(target[0], target[1], i)] == 1
        ]
    return candidate
