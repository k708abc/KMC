# cython: language_level=3, boundscheck=False, wraparound=False

from InputParameter cimport Params

cdef rec_ppt(
    Params params,
    int minute,
    int second,
    list img_names,
    list hist_names,
    list time,
    list coverage,
    str dir_name,
    double time_per_dep,
    list growth_mode,
    list mode_val_or,
    list other_modes,
)
