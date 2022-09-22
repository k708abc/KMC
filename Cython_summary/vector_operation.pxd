# distutils: language = c++
# cython: language_level=3, boundscheck=False, wraparound=False
# cython: cdivision=True

from libcpp.vector cimport vector

cdef bint search(vector[int] vec_1, int target)
cdef vector[int] remove_element(vector[int] vec_2, vector[int] remove_target)
cdef vector[int] remove_duplicate(vector[int] vec_3)