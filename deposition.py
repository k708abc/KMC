from typing import List, Tuple, Dict
import random
import math


def find_candidates(atom_set: Dict, bonds: Dict) -> list:
    candidate: List[Tuple[int, int, int]] = []  # リストの中身は into のタプルということでＯＫ？→OKです
    for index, condition in atom_set.items():
        # atom_set の key , val は何を意味しているの？
        # key val では無くて、（繰り返しになりますが。）
        # 変数名に意味を持たせるようにするのが大事です。直也のコード全般にいえることですが。
        # 適宜修正します（川上）
        if condition != 0:
            pass
        elif index[2] == 0:
            candidate.append(index)
        else:
            for bond in bonds[index]:
                if atom_set[bond] != 0:
                    candidate.append(index)
                    break
    return list(set(candidate))


def dep_position(candidate: List) -> tuple:
    random_n = random.random()
    num = math.floor(random_n * len(candidate))
    return candidate[num]


def judge_type(atom_set: Dict, bonds: Dict, dep_pos: tuple) -> int:
    # 周囲の原子の様子等に応じて原子の状態を区別
    # silicene or diamond etc...
    return 2


def deposit_an_atom(atom_set: Dict, bonds: Dict) -> tuple:
    candidate = find_candidates(atom_set, bonds)
    dep_pos = dep_position(candidate)
    atom_type = judge_type(atom_set, bonds, dep_pos)

    # rejection free ではイベント更新を行うサイトのリストが必要（あとまわし）
    # 蒸着サイトと関係する周囲のサイト

    return dep_pos, atom_type
