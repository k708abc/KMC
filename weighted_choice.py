import random
from typing import List, Tuple


def choice(events: List, rates: List, states: List) -> Tuple[Tuple, int]:
    tot = sum(rates)
    r_tot = random.random() * tot

    for rate, event, state in zip(rates, events, states):
        tot -= rate
        if r_tot >= tot:
            return event, state

    if r_tot < tot:
        print(len(events))
        print(len(rates))
        print(len(states))
        raise RuntimeError("Something wrong happens. Check weighted choice.")


if __name__ == "__main__":
    events: List = [(0, 0, 0), (0, 0, 1), (0, 0, 2)]
    rates: List = [0.1, 0.5, 10]
    states: List = [0, 1, 2]
    result: List = [0, 0, 0]
    repetition = 100
    #
    for _ in range(repetition):
        event, state = choice(events, rates, states)
        result[state] += 1
    print("Test: Weighted choice")
    print("repetition : " + str(repetition))
    print("rates = " + str(rates))
    print("results = " + str(result))
