from collections import defaultdict
import numpy


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

        # store transitions as a map
        self.transition_map = dict()
        # 2D matrix (After conforming the topology, create one matrix for visualization)
        self.transition_matrix = None

        # Put start and end in the states
        self.add_states(self.start, self.end)

        # edges
        self.edges    = []
        self.n_edges  = 0


        self.dynamic_table = []

    def add_model(self, model):
        pass

    def add_state(self, state):
        self.states.append(state)
        self.n_states += 1
        # initialize transition map
        self.transition_map[state.name] = defaultdict(lambda: 0, type=float)

    def add_states(self, *states):
        for state in states:
            self.add_state(state)

    def add_edge(self, a, b):
        '''

        :param a: state1
        :param b: state2
        :return:
        '''
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
        transition_log_probabilities = numpy.zeros( (m, m) ) + numpy.NINFNEGINF

        for i in range(m):
          for n in range( self.out_edge_count[i], self.out_edge_count[i+1] ):
            transition_log_probabilities[i, self.out_transitions[n]] = \
              self.out_transition_log_probabilities[n]

        return transition_log_probabilities

    def add_transition(self, from_state, to_state, probability, pseudocount=None):
        if from_state not in self.states:
            print ("ERROR: No such state named {}".format(from_state.name))
            raise Exception("No such state")
        elif to_state not in self.states:
            print ("ERROR: No such state named {}".format(from_state.name))
            raise Exception("No such state")
        else:
            self.transition_map[from_state.name][to_state.name] = probability

    def add_transitions(self, transitions):
        pass

    def bake(self):
        # set the topology and sort by the state name
        pass

    def viterbi(self):
        pass

    def dense_transition_matrix(self):
        pass

    def start_index(self):
        pass

    def end_index(self):
        pass

    def from_matrix(self, mat, distributions, starts, ends, name, state_names, merge=None):
        pass

    def concatenate(self, model):
        pass
  
    def fit(self, fit_patterns, algorithm='viterbi', transition_pseudocount=1, use_pseudocount=True):
        pass

if __name__ == "__main__":
    pass
