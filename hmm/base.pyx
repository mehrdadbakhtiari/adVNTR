
cdef class DiscreteDistribution(object):

    def __init__(self, emission):
        self.emission = emission

    def __getitem__(self, key):
        return self.emission[chr(key)]


cdef class State(object):

    def __init__(self, distribution=None, name=None):
        self.distribution = distribution
        self.name = name or str(id(self))

    def is_silent(self):
        if self.distribution is None:
            return 1
        else:
            return 0