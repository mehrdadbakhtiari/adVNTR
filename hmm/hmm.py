from collections import defaultdict
import numpy as np


class DiscreteDistribution:

    def __init__(self, emission):
        self.emission = emission

class State:

    def __init__(self, name=None):
        self.name = str(name)
        self.outgoing_edges = []

    def add_outgoing_edge(self, state):
        self.outgoing_edges.append(state)

    # State has all the outgoing edges
    edges = dict()


class Model:
    """ Hidden Markov Model """

    def __init__(self, name=None, start=None, end=None):
        # Save the name or make up a name.
        self.name = str(name) or str(id(self))
        self.model = "HiddenMarkovModel"

        self.start = start or State(name=self.name + "-start")
        self.end = end or State(name=self.name + "-end")

        self.number_of_edges = 0
        self.number_of_states = 0

        # Put start and end in the state
        self.states = set()
        self.states.add(self.start)
        self.states.add(self.end)

        # 2D matrix (After conforming the topology, create one matrix for visualization)
        self.transition_matrix = None

        self.states = set()
        self.dynamic_table = []

    def add_model(self, model):
        pass

    def add_state(self, state):
        pass

    def add_states(self, states):
        pass

    def states(self):
        pass

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
  
    def fit(self,fit_patterns, algorithm='viterbi', transition_pseudocount=1, use_pseudocount=True)
        pass

if __name__ == "__main__":
    pass
