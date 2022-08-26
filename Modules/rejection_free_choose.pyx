import random

# import decimal


cdef int choose_an_event(double r_tot, list event_rates):
    cdef int i
    cdef double event_rate 
    for i, event_rate in enumerate(event_rates):
        if event_rate >= r_tot:
            return i
        else:
            r_tot -= event_rate


cpdef tuple rejection_free_choise(
    double total_event_time,
    dict event_time,
    list event_time_tot,
    dict list_site_correspondance,
    double dep_rate,
):
    cdef double random_val, r_tot, rate
    cdef int i, event_number
    cdef tuple rate_site

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
