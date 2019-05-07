from collections import defaultdict
import numpy as np
from operator import attrgetter


class DiscreteDistribution:

    def __init__(self, emission):
        self.emission = emission

    def __getitem__(self, key):
        return self.emission[key]

class State:

    def __init__(self, distribution=None, name=None):
        self.distribution = distribution
        self.name = name or str(id(self))

    def is_silent(self):
        return (self.distribution == None)

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
            state = self.states[n]
            prob_mat[n,0] = self.transition_map[self.start.name][state.name] * state.distribution[seq[0]]
        for t in range(1,T):
            prob_mat[:,t%2] = 0.0
            for n in range(N):
                state = self.states[n]
                for n_prev in range(N):
                    state_prev = self.states[n_prev]
                    prob_mat[n,t%2] += prob_mat[n_prev,(t-1)%2] * self.transition_map[state_prev.name][state.name]
                prob_mat[n,t%2] *= state.distribution[seq[t]]
        for n in range(N):
            state = self.states[n]
            prob_mat[n,(T-1)%2] *= self.transition_map[state.name][self.end.name]
        prob = sum(prob_mat[:,(T-1)%2])
        return np.log(prob)

    def bake(self):
        # set the topology and sort by the state name
        """
        Right now the normal states are sorted by their name, 
        and the silent states are sorted by the order or occurence
        
        setting start_index and end_index
        setting n_states
        """
        silent_states, normal_states = [], []

        for state in self.states:
            if state.is_silent():
                silent_states.append(state)
            else:
                normal_states.append(state)

        np.random.seed(0)
        #random.seed(0)

        normal_states = list(sorted( normal_states, key=attrgetter('name')))

        self.states = normal_states + silent_states
        self.n_states = len(self.states)

        indices = { self.states[i]: i for i in range(self.n_states) }

        self.start_index = indices[self.start]
        self.end_index = indices[self.end]

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

        m = len(self.states)
        transition_probabilities = np.zeros( (m, m) ) 

        for i in range(m):
            state1 = self.states[i]
            for n in range(m):
                state2 = self.states[n]
                if (self.transition_map[state1.name][state2.name] > 0.0):
                  transition_probabilities[i, n] = self.transition_map[state1.name][state2.name]

        return transition_probabilities

    @classmethod
    def from_matrix(cls, transition_probabilities, distributions, starts, ends=None,
        state_names=None, name=None, verbose=False, merge='All' ):
        """Create a model from a more standard matrix format.

        Take in a 2D matrix of floats of size n by n, which are the transition
        probabilities to go from any state to any other state. May also take in
        a list of length n representing the names of these nodes, and a model
        name. Must provide the matrix, and a list of size n representing the
        distribution you wish to use for that state, a list of size n indicating
        the probability of starting in a state, and a list of size n indicating
        the probability of ending in a state.

        Parameters
        ----------
        transition_probabilities : array-like, shape (n_normal_states, n_normal_states)
            The probabilities of each state transitioning to each other state.

        distributions : array-like, shape (n_normal_states)
            The distributions for each state. Silent states are indicated by
            using None instead of a distribution object.

        starts : array-like, shape (n_normal_states)
            The probabilities of starting in each of the states.

        ends : array-like, shape (n_normal_states), optional
            If passed in, the probabilities of ending in each of the states.
            If ends is None, then assumes the model has no explicit end
            state. Default is None.

        state_names : array-like, shape (n_normal_states), optional
            The name of the states. If None is passed in, default names are
            generated. Default is None

        name : str, optional
            The name of the model. Default is None

        verbose : bool, optional
            The verbose parameter for the underlying bake method. Default is False.

        merge : 'None', 'Partial', 'All', optional
            The merge parameter for the underlying bake method. Default is All

        Returns
        -------
        model : Model
            The baked model ready to go.

        Examples
        --------
        matrix = [ [ 0.4, 0.5 ], [ 0.4, 0.5 ] ]
        distributions = [NormalDistribution(1, .5), NormalDistribution(5, 2)]
        starts = [ 1., 0. ]
        ends = [ .1., .1 ]
        state_names= [ "A", "B" ]

        model = Model.from_matrix( matrix, distributions, starts, ends,
            state_names, name="test_model" )
        """

        # Build the initial model
        model = Model( name=name )
        state_names = state_names or ["s{}".format(i) for i in range(len(distributions))]

        # Build state objects for every state with the appropriate distribution
        states = [ State( distribution, name=name ) for name, distribution in
            zip( state_names, distributions) ]

        n = len( states )

        # Add all the states to the model
        for state in states:
            model.add_state( state )

        # Connect the start of the model to the appropriate state
        for i, prob in enumerate( starts ):
            if prob != 0:
                model.add_transition( model.start, states[i], prob )

        # Connect all states to each other if they have a non-zero probability
        for i in range( n ):
            for j, prob in enumerate( transition_probabilities[i] ):
                if prob != 0.:
                    model.add_transition( states[i], states[j], prob )

        if ends is not None:
            # Connect states to the end of the model if a non-zero probability
            for i, prob in enumerate( ends ):
                if prob != 0:
                    model.add_transition( states[j], model.end, prob )

        model.bake()
        return model

    def concatenate(self, model):
        pass
  
    def viterbi(self):
        pass

    #def from_matrix(self,mat, distributions, starts, ends, name=name, state_names=state_names, merge=None):
    #    pass

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
