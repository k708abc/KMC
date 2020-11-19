from typing import List
import random
import math


def get_candidates(atom_set, bonds):
    candidate: List[tuple] = []
    for key, val in atom_set.items():
        if key[2] == 0:
            candidate.append(key)
        else:
            for bond in bonds[key]:
                if atom_set[bond] != 0:
                    candidate.append(key)
                    break
    return list(set(candidate))


def dep_position(candidate):
    random_n = random.random()
    num = math.floor(random_n * len(candidate))
    return candidate[num]


def judge_type(atom_set, bonds, dep_pos):
    # 周囲の原子の様子等に応じて原子の状態を区別
    # silicene or diamond etc...
    return 2


def deposit_an_atom(atom_set: dict, bonds: dict):
    candidate = get_candidates(atom_set, bonds)
    dep_pos = dep_position(candidate)
    atom_type = judge_type(atom_set, bonds, dep_pos)

    # rejection free ではイベント更新を行うサイトのリストが必要（あとまわし）
    # 蒸着サイトと関係する周囲のサイト

    return dep_pos, atom_type
