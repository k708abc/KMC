import random
from typing import List, Tuple


def choose_atom(
    candidates: List,
) -> Tuple[int, int, int]:  # Tuple の中身？ Tuple[int, int, int] とか?
    return random.choice(candidates)


if __name__ == "__main__":
    num_cand = 10
    candidates = [(0, 0, i) for i in range(num_cand)]
    results = [0 for _ in range(num_cand)]
    for _ in range(1000):
        result = choose_atom(candidates)
        results[result[2]] += 1
    print("Test choose atom")
    print("candidates : " + str(candidates))
    print("results : " + str(results))
