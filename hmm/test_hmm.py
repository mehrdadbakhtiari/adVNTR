from collections import defaultdict

from hmm import State
from hmm import Model
from hmm import DiscreteDistribution

def test_hmm_add_transition():
    pass

def test_hmm_distribution():
    pass

def test_hmm_state():
    transitions = defaultdict(lambda: defaultdict(float))
    emissions = defaultdict(lambda: defaultdict(float))

    # State(insert_distribution, name='I%s_%s' % (i, repeat))
    distribution = {'A':0.2, 'C':0.3, 'G':0.3, 'T':0.2}
    state1 = State(distribution, name="intron")
    state2 = State(distribution)

    emissions['start']['intron'] = 0.5
    emissions['start']['exon'] = 0.5
    emissions['intron']['intron'] = 0.6
    emissions['intron']['exon'] = 0.4
    emissions['exon']['intron'] = 0.6
    emissions['exon']['exon'] = 0.4

    print("state1 name: ", state1.name)
    print("state2 name: ", state2.name)
    print("emissions: ", emissions)

def test_hmm_model():
    hmm = Model(name="HiddenMarkovModel")

    distribution = {'A':0.2, 'C':0.3, 'G':0.3, 'T':0.2}
    state1 = State(distribution, name="s1")
    state2 = State(distribution, name="s2")

    hmm.add_states(state1, state2)
    hmm.add_edge(hmm.start, state1)
    hmm.add_edge(state1, state2)
    hmm.add_edge(state2, hmm.start)

    print("states: ", hmm.states)
    for state in hmm.states:
        print("state name: ", state.name)
    print("n states: ", hmm.n_states)
    print("n edges: ", hmm.n_edges)
    pass

if __name__ == "__main__":
    test_hmm_state()
    test_hmm_model()
