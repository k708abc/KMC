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
        raise RuntimeError("Something wrong happens. Check weighted choice.")
