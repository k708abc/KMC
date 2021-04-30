import sys

sys.path.append("../../KMC/")
from Modules.choose_site import choose_atom

if __name__ == "__main__":
    print("Chek choose_site.py")
    num_cand = 10
    repetition = 1000
    candidates = [(0, 0, i) for i in range(num_cand)]
    results = [0 for _ in range(num_cand)]
    for _ in range(repetition):
        result = choose_atom(candidates)
        results[result[2]] += 1
    print("repetition : " + str(repetition))
    print("candidates : " + str(candidates))
    print("results : " + str(results))