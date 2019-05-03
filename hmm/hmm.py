from collections import defaultdict
import numpy as np


class DiscreteDistribution:

    def __init__(self, emission):
        self.emission = emission

    def __getitem__(self, key):
        return self.emission[key]

class State:

    def __init__(self, distribution=None, name=None):
        self.distribution = distribution
        self.name = name or str(id(self))

class Model:
    """ Hidden Markov Model 
        start: a State representing the model start 
        end:   a State representing the model end
        states: a list of states
        edges:  a list of edges represented by tuples of (from-state, to-state)


    """

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
            if isinstance( state, list ):
                for s in state:
                    self.add_state( s )
            else:
                self.add_state( state )

    def state_count(self):
        return self.n_states

    def add_transition(self, from_state, to_state, probability, pseudocount=None):
        if from_state not in self.states:
            print ("ERROR: No such state named {}".format(from_state.name))
            raise Exception("No such state")
        elif to_state not in self.states:
            print ("ERROR: No such state named {}".format(to_state.name))
            raise Exception("No such state")
        else:
            self.transition_map[from_state.name][to_state.name] = probability

    def add_transitions(self, transitions):
        pass

    def log_probability(self, seq):
        T = len(seq)
        N = self.state_count() - 2 # exclude the first two states (model-start and model-end)
        prob_mat = np.zeros((N,2))
        for n in range(N):
            state = self.states[n+2]
            prob_mat[n,0] = self.transition_map[self.states[0].name][state.name] * state.distribution[seq[0]]
        for t in range(1,T):
            prob_mat[:,t%2] = 0.0
            for n in range(N):
                state = self.states[n+2]
                for n_prev in range(N):
                    state_prev = self.states[n_prev+2]
                    prob_mat[n,t%2] += prob_mat[n_prev,(t-1)%2] * self.transition_map[state_prev.name][state.name]
                prob_mat[n,t%2] *= state.distribution[seq[t]]
        for n in range(N):
            state = self.states[n+2]
            prob_mat[n,(T-1)%2] *= self.transition_map[state.name][self.states[1].name]
        prob = sum(prob_mat[:,(T-1)%2])
        return np.log(prob)

    def bake(self):
        # set the topology and sort by the state name
        pass

    def viterbi(self):
        pass

    def dense_transition_matrix( self ):
        """Returns the dense transition matrix.

        Parameters
        ----------
        None

        Returns
        -------
        matrix : numpy.ndarray, shape (n_states, n_states)
            A dense transition matrix, containing the probability
            of transitioning from each state to each other state.
        """

        m = len(self.states)-2
        transition_probabilities = numpy.zeros( (m, m) ) 

        for i in xrange(m):
            state1 = self.states[i+2]
            for n in xrange(m):
                state2 = self.states[n+2]
                if (self_transition_map[state1.name][state2.name] > 0.0):
                  transition_probabilities[i, n] = self_transition_map[state1.name][state2.name]

        return transition_probabilities

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

    #def add_edge(self, a, b):
    #    '''

    #    :param a: state1
    #    :param b: state2
    #    :return:
    #    '''
    #    self.edges.append( ( a, b ) )
    #    self.n_edges += 1

    #def add_transition( self, a, b ):
    #    self.add_edge( a, b )

    #def edge_count(self):
    #    return self.n_edges


if __name__ == "__main__":
    pass
