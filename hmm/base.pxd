# base.pxd

cimport numpy


cdef class DiscreteDistribution(object):
    cdef public dict emission


cdef class State(object):
    cdef public DiscreteDistribution distribution
    cdef public str name

    #cdef public int is_silent(self)
