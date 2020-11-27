import random
from typing import List, Tuple


def choose_atom(candidates: List) -> Tuple:  # Tuple の中身？ Tuple[int, int, int] とか?
    return random.choice(candidates)
