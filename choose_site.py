import random
from typing import List, Tuple


def choose_atom(
    candidates: List[Tuple[int, int, int]],
) -> Tuple[int, int, int]:
    return random.choice(candidates)


if __name__ == "__main__":
    num_cand = 10
    repetition = 1000
    candidates = [(0, 0, i) for i in range(num_cand)]
    results = [0 for _ in range(num_cand)]
    for _ in range(repetition):
        result = choose_atom(candidates)
        results[result[2]] += 1
    print("Test : choose an atom from a list")
    print("repetition : " + str(repetition))
    print("candidates : " + str(candidates))
    print("results : " + str(results))
