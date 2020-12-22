import random
from typing import List, Tuple


def choose_an_event(r_tot, event_rates):
    for i, event_rate in enumerate(event_rates):
        if event_rate >= r_tot:
            return i
        else:
            r_tot -= event_rate

    print("value remains")
    print("rtot in event num = " + str(r_tot))
    print(event_rates)
    return len(event_rates) - 1


def rejection_free_choise(
    total_event_time: float, event_time: List[List[float]], event_time_tot: List[float]
):
    r_tot = total_event_time * random.random()
    """
    tot_time = 0
    for sites, times in event_time_tot.items():
        tot_time += times
    print("total_event_time : " + str(total_event_time))
    print("tot_time_calculatesd : " + str(tot_time))
    print("diff: " + str(total_event_time - tot_time))
    """

    for rate_site, rate in event_time_tot.items():
        if rate >= r_tot:
            event_number = choose_an_event(r_tot, event_time[rate_site])
            if event_number is None or event_number == -1:
                print("None happenes in rejection free choice")
                print("rate site = " + str(rate_site))
                print("r_tot = " + str(r_tot))
                print("rate = " + str(rate))
                print("event time = " + str(event_time[rate_site]))

            return rate_site, event_number
        else:
            r_tot -= rate
    if r_tot >= 0:
        return (-1, -1, -1), 0
    else:
        print("Something wrong in rejection free")
        input()
