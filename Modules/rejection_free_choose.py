import random
from typing import List, Tuple
import decimal


def choose_an_event(r_tot: float, event_rates: List[float]) -> int:
    for i, event_rate in enumerate(event_rates):
        if event_rate >= r_tot:
            return i
        else:
            r_tot -= event_rate


def rejection_free_choise(
    total_event_time: float, event_time: List[List[float]], event_time_tot: List[float]
):
    random_val = decimal.Decimal(random.random())
    r_tot = total_event_time * random_val

    for rate_site, rate in event_time_tot.items():
        if rate >= r_tot:
            if rate == 0:
                pass
            else:
                event_number = choose_an_event(r_tot, event_time[rate_site])
                return rate_site, event_number
        else:
            r_tot -= rate
    if r_tot >= 0:
        return (-1, -1, -1), 0
    else:
        print("Something wrong in rejection free")
        input()
