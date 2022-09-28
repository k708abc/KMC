# cython: language_level=3, boundscheck=False, wraparound=False

cdef int choose_an_event(double r_tot, list event_rates)

cdef tuple rejection_free_choise(
    double total_event_time,
    list event_time,
    list event_time_tot,
    double dep_rate,
)
