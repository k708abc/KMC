import random
from typing import List, Tuple, Dict

# import decimal


def choose_an_event(r_tot: float, event_rates: List[float]) -> int:
    for i, event_rate in enumerate(event_rates):
        if event_rate >= r_tot:
            return i
        else:
            r_tot -= event_rate


def rejection_free_choise(
    total_event_time: float,
    event_time: Dict[Tuple, List[float]],
    event_time_tot: List[float],
    list_site_correspondance: Dict[int, Tuple],
    dep_rate: float,
):
    # random_val = decimal.Decimal(random.random())
    random_val = random.random()
    r_tot = total_event_time * random_val

    for i, rate in enumerate(event_time_tot):
        if rate >= r_tot:
            if rate == 0:
                pass
            else:
                rate_site = list_site_correspondance[i]
                event_number = choose_an_event(r_tot, event_time[rate_site])
                return rate_site, event_number
        else:
            r_tot -= rate
    if r_tot <= dep_rate:
        return (-1, -1, -1), 0
    else:
        print("total: " + str(total_event_time))
        print("r_tot: " + str(r_tot))
        print("dep_rate: " + str(dep_rate))
        print("random: " + str(random_val))
        return (-1, -1, -1), 0
