from collections import defaultdict
import numpy as np


class DiscreteDistribution:

    def __init__(self, emission):
        self.emission = emission

class State:

    def __init__(self, distribution=None, name=None):
        self.distribution = distribution
        self.name = name or str(id(self))

class Model:
    """ Hidden Markov Model """

    def __init__(self, name=None, start=None, end=None):
        # Save the name or make up a name.
        self.name = str(name) or str(id(self))
        self.model = "HiddenMarkovModel"

        # states
        self.states   = []
        self.n_states = 0
        self.start    = start or State(name=self.name + "-start")
        self.end      = end or State(name=self.name + "-end")

        # Put start and end in the states
        self.add_states(self.start, self.end)

        # edges
        self.edges    = []
        self.n_edges  = 0

        # 2D matrix (After conforming the topology, create one matrix for visualization)
        self.transition_matrix = None

        self.dynamic_table = []

    def add_model(self, model):
        pass

    def add_state(self, state):
        self.states.append(state)
        self.n_states += 1

    def add_states(self, *states):
        for state in states:
            self.add_state(state)

    def add_edge(self, a, b):
        self.edges.append( ( a, b ) )
        self.n_edges += 1

    def add_transition( self, a, b ):
        self.add_edge( a, b )

    def state_count(self):
        return self.n_states

    def edge_count(self):
        return self.n_edges

    def dense_transition_matrix(self):
        """
        Returns the dense transition matrix. Useful if the transitions of
        somewhat small models need to be analyzed.
        """

        m = len(self.states)
        transition_log_probabilities = numpy.zeros( (m, m) ) + NEGINF

        for i in range(m):
          for n in range( self.out_edge_count[i], self.out_edge_count[i+1] ):
            transition_log_probabilities[i, self.out_transitions[n]] = \
              self.out_transition_log_probabilities[n]

        return transition_log_probabilities

    def add_transition(self, from_state, to_state, probability, pseudocount):
        pass

    def add_transitions(self, transitions):
        pass

    def bake(self):
        pass

    def viterbi(self):
        pass

    def dense_transition_matrix(self):
        pass

    def start_index(self):
        pass

    def end_index(self):
        pass

    def from_matrix(self,mat, distributions, starts, ends, name, state_names, merge=None):
        pass

    def concatenate(self,model):
        pass
  
    def fit(self,fit_patterns, algorithm='viterbi', transition_pseudocount=1, use_pseudocount=True):
        pass

if __name__ == "__main__":
    pass
