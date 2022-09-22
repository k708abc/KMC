# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False


from libcpp.vector cimport vector

cdef int choose_an_event(double r_tot, vector[double] event_rates)

cdef (int, int) rejection_free_choise(
    double total_event_time,
    vector[vector[double]] event_time,
    vector[double] event_time_tot,
    double dep_rate,
)
