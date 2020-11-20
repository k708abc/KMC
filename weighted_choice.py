import random


def choice(events, rates):
    tot = sum(rates)
    n_r = random.random() * tot

    for rate, event in zip(rates, events):
        tot -= rate
        if n_r >= tot:
            return event
