from libc.math cimport log

cdef class DiscreteDistribution(object):

    def __init__(self, emission):
        self.emission = emission

        key_map = {'A':0, 'C':1, 'G':2, 'T':3}
        self.log_emission = {key_map[key]: log(value) for key, value in emission.items()}
        # print(self.log_emission)

    def __getitem__(self, key):
        return self.log_emission[key]
        # return self.emission[key]


cdef class State(object):

    def __init__(self, distribution=None, name=None):
        self.distribution = distribution
        self.name = name or str(id(self))

    def is_silent(self):
        if self.distribution is None:
            return 1
        else:
            return 0