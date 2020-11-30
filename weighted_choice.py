import random
from typing import List, Tuple


def choice(events: List, rates: List, states: List) -> Tuple:
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
