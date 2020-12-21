from typing import List, Tuple, Dict
import random
import math
from event_collection import state_after_move
from read_examples import read_atom_set, read_bonds, read_lattice
from InputParameter import Params
from deposition_check import dep_check_poscar


def find_candidates(atom_set: Dict, bonds: Dict) -> List[Tuple[int, int, int]]:
    candidate: List[Tuple[int, int, int]] = []  # リストの中身は into のタプルということでＯＫ？→OKです
    for atom_index, condition in atom_set.items():
        # atom_set の key , val は何を意味しているの？
        # key val では無くて、（繰り返しになりますが。）
        # 変数名に意味を持たせるようにするのが大事です。直也のコード全般にいえることですが。
        # 適宜修正します（川上）
        if condition != 0:
            pass
        elif atom_index[2] == 0:
            candidate.append(atom_index)
        else:
            for bond in bonds[atom_index]:
                if atom_set[bond] != 0:
                    candidate.append(atom_index)
                    break
    return list(set(candidate))


def dep_position(candidate: List) -> Tuple[int, int, int]:
    random_n = random.random()
    num = math.floor(random_n * len(candidate))
    return candidate[num]


def judge_type(atom_set: Dict, bonds: Dict, dep_pos: Tuple, params) -> int:
    # 周囲の原子の様子等に応じて原子の状態を区別
    # silicene or diamond etc...
    if params.trans_check is True:
        state = state_after_move(atom_set, bonds, dep_pos, params)
        return state
    else:
        return 2


def remove_first(candidate) -> List[Tuple]:
    candidate_new = []
    for site in candidate:
        if site[2] not in (0, 1):
            candidate_new.append(site)
    return candidate_new


def deposit_an_atom(atom_set: Dict, bonds: Dict, params) -> Tuple:
    defect = params.keep_defect_check
    empty_first = int(params.put_first)
    candidate = find_candidates(atom_set, bonds)
    if (defect is True) and (empty_first == int(params.num_defect)):
        candidate = remove_first(candidate)

    dep_pos = dep_position(candidate)
    atom_type = judge_type(atom_set, bonds, dep_pos, params)

    # rejection free ではイベント更新を行うサイトのリストが必要（あとまわし）
    # そこでは、蒸着サイトと関係する周囲のサイトをイベント更新リストに入れる
    return dep_pos, atom_type


def highest_z(atom_set: Dict):
    maxz = 1
    for index, state in atom_set.items():
        if (state != 0) and (index[2] + 1 > maxz):
            maxz = index[2] + 1
    return maxz


if __name__ == "__main__":
    defect = True
    empty_first = 2
    parameter = Params()
    unit_length = parameter.n_cell_init
    #
    atom_set = read_atom_set()
    bonds = read_bonds()
    lattice = read_lattice()
    maxz = highest_z(atom_set)
    candidate = find_candidates(atom_set, bonds)
    if (defect is True) and (empty_first == 1):
        candidate = remove_first(candidate)
    dep_check_poscar(atom_set, candidate, lattice, unit_length, maxz)
    print("Test candidates for depositing an atom")
    print("candidates :" + str(candidate))
    print("A poscar is formed to check the candidate of deposition")
