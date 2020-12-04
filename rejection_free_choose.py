import random
from typing import List, Tuple


def choose_an_event(r_tot, event_rates):
    for i, event_rate in enumerate(event_rates):
        if event_rate >= r_tot:
            return i
        else:
            r_tot -= event_rate
    if r_tot > 0:
        print("value remains")


def rejection_free_choise(total_event_time, event_time, event_time_tot):
    r_tot = total_event_time * random.random()
    for rate_site, rate in event_time_tot.items():
        if rate >= r_tot:
            event_number = choose_an_event(r_tot, event_time[rate_site])
            return rate_site, event_number
        else:
            r_tot -= rate
    if r_tot >= 0:
        return (-1, -1, -1), 0
