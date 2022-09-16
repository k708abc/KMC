# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False

from libcpp.vector cimport vector

cdef bint search(vector[int] vec_1, int target):
    cdef int i
    for i in vec_1:
        if i == target:
            return True
    return False

cdef vector[int] remove_element(vector[int] vec_2, vector[int] remove_target):
    cdef vector[int] vec_removed = []
    cdef int i
    for i in vec_2:
        if search(remove_target, i):
            pass
        else:
            vec_removed.push_back(i)
    return vec_removed

cdef vector[int] remove_duplicate(vector[int] vec_3):
    cdef vector[int] vec_removed = []
    cdef int i
    for i in vec_3:
        if search(vec_removed, i) is False:
            vec_removed.push_back(i)
    return vec_removed