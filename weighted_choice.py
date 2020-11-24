import random
from typing import List


def choice(events: List, rates: List) -> tuple:
    tot = sum(rates)
    r_tot = random.random() * tot

    for rate, event in zip(rates, events):
        tot -= rate
        if r_tot >= tot:
            return event

    if r_tot < tot:
        raise RuntimeError("Something wrong happens. Check weighted choice.)
