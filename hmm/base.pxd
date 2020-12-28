# base.pxd


cdef class DiscreteDistribution(object):
    cdef public dict emission
    cdef public dict log_emission

cdef class State(object):
    cdef public DiscreteDistribution distribution
    cdef public str name

    #cdef public int is_silent(self)
