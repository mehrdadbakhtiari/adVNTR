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

        self.start_index = -1     # will be set in bake
        self.end_index = -1       # will be set in bake

        self.dynamic_table = []

        self.subModels = [self]
        self.n_subModels = 1

        self.is_baked = False

    def add_model(self, model):
        pass

    def add_state(self, state):
        self.states.append(state)
        self.n_states += 1
        # initialize transition map
        self.transition_map[state] = defaultdict(lambda: 0)

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
            self.transition_map[from_state][to_state] = probability

    def add_transitions(self, transitions):
        pass

    def log_probability(self, seq):
        T = len(seq)
        N = self.state_count() - 2 # exclude the first two states (model-start and model-end)
        prob_mat = np.zeros((N,2))
        for n in range(N):
            state = self.states[n]
            print("self.start.name: ", self.start.name)
            print("state.name: ", state.name)
            prob_mat[n,0] = self.transition_map[self.start][state] * state.distribution[seq[0]]
        for t in range(1,T):
            prob_mat[:,t%2] = 0.0
            for n in range(N):
                state = self.states[n]
                for n_prev in range(N):
                    state_prev = self.states[n_prev]
                    prob_mat[n,t%2] += prob_mat[n_prev,(t-1)%2] * self.transition_map[state_prev][state]
                prob_mat[n,t%2] *= state.distribution[seq[t]]
        for n in range(N):
            state = self.states[n]
            prob_mat[n,(T-1)%2] *= self.transition_map[state][self.end]
        prob = sum(prob_mat[:,(T-1)%2])
        return np.log(prob)

    def bake(self, merge=None):
        # set the topology and sort by the state name
        """
        In a model, start state comes the first and end state comes the last.
        Other states are in the middle, and they are sorted by their name.
        e.g.)
        start - I0 - D1 - M1, I1 - D2 - M2 - I2, ... D10 - M10 - I10 - end

        setting start_index and end_index

        setting connections between subModels
        """
        np.random.seed(0)
        #random.seed(0)

        # bake all the subModels
        for subModel in self.subModels:
            # Ordering states
            states_without_start_and_end = [state for state in subModel.states if state is not subModel.start and state is not subModel.end ]
            sorted_states = list(sorted( states_without_start_and_end, key=attrgetter('name')))

            subModel.states = [subModel.start] + sorted_states + [subModel.end]
            subModel._sort_states()
            indices = { subModel.states[i]: i for i in range(subModel.n_states) }

            subModel.start_index = indices[subModel.start]
            subModel.end_index = indices[subModel.end]

        self.is_baked = True

    def _sort_states(self):
        """
        Sort states in pre-defined (topology of our hmm model) order.

        State naming rule:
        There should be three types of state except start and end
        1. Insert
        2. Match
        3. Delete
        Insertion states start with I
        Match states start with M
        Delete states start with D

        Format:
        I/M/D[index]_[repeating_unit_index]

        Example:

        total_hmm_start
        --------------------------
        suffix_matcher_hmm_start
        ...
        suffix_matcher_hmm_end
        --------------------------
        unit_start (repeating unit matcher hmm)
        --------------------------
        repeat_start_1
        I0_1
        D1_1
        M1_1
        I1_1
        ...
        repeat_end_1
        --------------------------
        repeat_start_2
        I0_2
        D1_2
        M1_2
        I1_2
        ...
        repeat_end_2
        --------------------------
        unit_end (repeating unit matcher hmm)
        --------------------------
        prefix_matcher_hmm_start
        ...
        prefix_matcher_hmm_end
        --------------------------
        total_hmm_end

        :return: None
        """

        sorted_states = []

        insert_states = []
        match_states = []
        delete_states = []

        for state in self.states:
            if state.name.startswith("I"):
                insert_states.append(state)
            elif state.name.startswith("M"):
                match_states.append(state)
            elif state.name.startswith("D"):
                delete_states.append(state)
            else:
                assert (state == self.start or state == self.end), "State type should be in {I, M, D, start, end}"

        insert_states.sort(key=lambda x: x.name)
        match_states.sort(key=lambda x: x.name)
        delete_states.sort(key=lambda x: x.name)

        # 1. Start state
        sorted_states.append(self.start)
        # 2. Insert 0 state (number of repeating units)
        sorted_states.append(insert_states.pop(0))

        # 3. Delete, Match, Insert states
        assert (len(match_states) == len(delete_states))

        for i in range(len(match_states)):
            sorted_states.append(delete_states[i])
            sorted_states.append(match_states[i])
            sorted_states.append(insert_states[i])

        # 4. End state
        sorted_states.append(self.end)

        self.states = sorted_states

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
                if (self.transition_map[state1][state2] > 0.0):
                  transition_probabilities[i, n] = self.transition_map[state1][state2]

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

    def concatenate( self, other, suffix='', prefix='', transition_probability=1.0 ):
        """Concatenate this model to another model.

        Concatenate this model to another model in such a way that a single
        probability 1 edge is added between self.end and other.start. Rename
        all other states appropriately by adding a suffix or prefix if needed.

        Parameters
        ----------
        other : HiddenMarkovModel
            The other model to concatenate

        suffix : str, optional
            Add the suffix to the end of all state names in the other model.
            Default is ''.

        prefix : str, optional
            Add the prefix to the beginning of all state names in the other
            model. Default is ''.

        Returns
        -------
        None
        """

        ## "self" must be baked
        #if (not self.is_baked):
        #    print("ERROR: to call concatenate,the my model must be baked.")
        #    raise ValueError

        ## "other" must be baked
        #if (not other.is_baked):
        #    print("ERROR: to call concatenate,the other model must be baked.")
        #    raise ValueError

        other.name = "{}{}{}".format( prefix, other.name, suffix )
        for state in other.states:
            state.name = "{}{}{}".format( prefix, state.name, suffix )

        self.subModels.append(other)
        self.n_subModels += 1

        # setting connections between subModels if they exist
        subModel = self.subModels[-1]
        subModel_prev = self.subModels[-2]
        subModel_prev.transition_map[subModel_prev.end][subModel.start] = transition_probability

        # when concatenate with another model we need another bake
        self.is_baked = False 

        ## add states of other 
        #self.add_states(other.states)
    
        ## set the n_states
        #self.n_states += other.n_states
 
        ## set transition map
        #for state in other.states:
        #    for key, prob in other.transition_map[state].items():
        #        if (prob != 0):
        #            self.transition_map[state][key] = prob

        ## set transitional link
        #self.add_transition( self.end, other.start, 1.00 )
        #self.end = other.end

    def viterbi(self, sequence):
        """
        
        :param sequence:
        :return: log probability and viterbi path
        """

        if not self.is_baked:
            print ("ERROR: To call viterbi, the model must have been baked")
            raise ValueError

        # Return variables
        logp = 0  # Log probability of the best path
        vpath = []  # Viterbi path

        # build aggregate states and transition map from subModels
        state_to_index = {}
        n_states = 0
        states = []
        transition_map = dict()
        for subModel in self.subModels:
            state_to_index.update( dict(zip(subModel.states, range(n_states,n_states+subModel.n_states))) )
            n_states += subModel.n_states
            for state in subModel.states:
                states.append(state)
            transition_map.update(subModel.transition_map)

        # Initialize dynamic table
        self.dynamic_table = np.zeros((n_states, len(sequence) + 1))
        self.dynamic_table[state_to_index[self.start]][0] = 1

        # Storing previous states (row) and column (Easy version)
        vpath_table = np.zeros((n_states, len(sequence) + 1), dtype=object)

        for col in range(len(sequence)):
            for row, state in enumerate(states):
                #print ("state: {}, row: {}, col: {}".format(state.name, row, col))
                # for all neighboring states, update the probability
                for neighbor_state in transition_map[state]:
                    #print ("neighbor_state: {}".format(neighbor_state.name))
                    if not state.is_silent():
                        # Emit a character and move to the next column
                        prob = self.dynamic_table[row][col] * \
                               transition_map[state][neighbor_state] * \
                               state.distribution[sequence[col]]
                        neighbor_state_index = state_to_index[neighbor_state]

                        if self.dynamic_table[neighbor_state_index][col + 1] < prob:
                            self.dynamic_table[neighbor_state_index][col + 1] = prob
                            # update path table
                            vpath_table[neighbor_state_index][col + 1] = [row, col]

                    else:  # Silent state: Stay in the same column
                        prob = self.dynamic_table[row][col] * \
                               transition_map[state][neighbor_state]
                        neighbor_state_index = state_to_index[neighbor_state]

                        if self.dynamic_table[neighbor_state_index][col] < prob:
                            self.dynamic_table[neighbor_state_index][col] = prob
                            # update path table
                            vpath_table[neighbor_state_index][col] = [row, col]

        # Back tracking viterbi path
        # Start from end_state
        # NOTE This should be changed when there are many end state
        vpath.append((state_to_index[self.subModels[-1].end], self.subModels[-1].end))
        end_index = state_to_index[self.subModels[-1].end]
        previous_pointer = vpath_table[end_index][-1]
        row = previous_pointer[0]
        col = previous_pointer[1]

        while True:
            vpath.append((state_to_index[states[row]], states[row]))
            previous_pointer = vpath_table[row][col]
            if previous_pointer == 0:
                break
            row = previous_pointer[0]
            col = previous_pointer[1]

        logp = np.log(self.dynamic_table[state_to_index[self.subModels[-1].end]][-1])
        vpath = vpath[::-1]

        return logp, vpath


if __name__ == "__main__":
    pass
