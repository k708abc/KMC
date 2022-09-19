# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from libcpp.vector cimport vector
from libcpp.pair cimport pair
import random

cdef int choose_an_event(double r_tot, vector[double] event_rates):
    cdef int i
    cdef double event_rate 
    for i, event_rate in enumerate(event_rates):
        if event_rate >= r_tot:
            return i
        else:
            r_tot -= event_rate


cdef pair[int, int] rejection_free_choise(
    double total_event_time,
    vector[vector[double]] event_time,
    vector[double] event_time_tot,
    double dep_rate,
):
    cdef double random_val, r_tot, rate
    cdef int i, event_number

    random_val = random.random()
    r_tot = total_event_time * random_val

    for i, rate in enumerate(event_time_tot):
        if rate >= r_tot:
            if rate == 0:
                pass
            else:
                event_number = choose_an_event(r_tot, event_time[i])
                return i, event_number
        else:
            r_tot -= rate
    if r_tot <= dep_rate:
        return -1, 0
    else:
        print("total: " + str(total_event_time))
        print("r_tot: " + str(r_tot))
        print("dep_rate: " + str(dep_rate))
        print("random: " + str(random_val))
        return -1, 0
