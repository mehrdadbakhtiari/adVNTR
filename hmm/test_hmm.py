from collections import defaultdict

from hmm import State
from hmm import Model

def test_hmm_add_transition():
    pass

def test_hmm_state():
    transitions = defaultdict(lambda: defaultdict(float))
    emissions = defaultdict(lambda: defaultdict(float))

    hmm = Model(name="HiddenMarkovModel")

    # State(insert_distribution, name='I%s_%s' % (i, repeat))
    state1 = State(name="intron")
    state2 = State(name="exon")

    emissions['start']['intron'] = 0.5
    emissions['start']['exon'] = 0.5
    emissions['intron']['intron'] = 0.6
    emissions['intron']['exon'] = 0.4
    emissions['exon']['intron'] = 0.6
    emissions['exon']['exon'] = 0.6


if __name__ == "__main__":
    test_hmm_add_transition()