# cython: language_level=3, boundscheck=False, wraparound=False

from InputParameter cimport Params

cdef int highest_z(list pos_all, list index_list)

cdef dir_formarion(str name)

cdef list[tuple] triangle(double xp, double yp, int z)

cdef list color_determinate(int z, int maxz)

cdef image_formaiton(list pos, list lattice, int length, int maxz, str img_name)

cdef list occupation_of_layers(int maxz, pos, int n_BL, list index_list)

cdef hist_formation(list pos, int maxz, int n_BL, str img_name, list index_list)

cdef rec_poscar(list pos, list lattice, int unit_length, int maxz, str rec_name, list index_list)

cdef rec_events_per_dep_fig(list atoms, list num_events, str fig_name)

cdef rec_events_per_dep(list num_events, list atoms, str rec_name)

cdef rec_growth_mode(list growth_mode, list coverage, Params params)

cdef delete_images(str dir_name)

cdef rec_yaml(Params init_value, str dir_name)

cdef record_data(
    list pos_all,
    list time,
    list coverage,
    list lattice,
    Params params,
    int minute,
    int second,
    time_per_eve,
    Params init_value,
    list index_list,
)