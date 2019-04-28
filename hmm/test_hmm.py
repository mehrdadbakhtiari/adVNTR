from collections import defaultdict

from hmm import State
from hmm import Model
from hmm import DiscreteDistribution

import unittest


class TestMethods(unittest.TestCase):

    def test_example_pomegranate(self):
        """
        This example is taken from https://pomegranate.readthedocs.io/en/latest/HiddenMarkovModel.html
        """

        from pomegranate import *
        d1 = DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        d2 = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        d3 = DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})

        s1 = State(d1, name="s1")
        s2 = State(d2, name="s2")
        s3 = State(d3, name="s3")

        model = Model(name='example')
        model.add_states([s1, s2, s3])
        model.add_transition(model.start, s1, 0.90)
        model.add_transition(model.start, s2, 0.10)
        model.add_transition(s1, s1, 0.80)
        model.add_transition(s1, s2, 0.20)
        model.add_transition(s2, s2, 0.90)
        model.add_transition(s2, s3, 0.10)
        model.add_transition(s3, s3, 0.70)
        model.add_transition(s3, model.end, 0.30)
        model.bake()

        print(model.log_probability(list('ACGACTATTCGAT')))
        # should be -22.73896159971087

        answer = model.log_probability(list('ACGACTATTCGAT'))
        expected = -22.73896159971087
        self.assertAlmostEqual(expected, answer)

        print(", ".join(state.name for i, state in model.viterbi(list('ACGACTATTCGAT'))[1]))
        # should be example - start, s1, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s3, example - end

        answer = ", ".join(state.name for i, state in model.viterbi(list('ACGACTATTCGAT'))[1])
        expected = "example - start, s1, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s3, example - end"
        self.assertEqual(expected, answer)


    def test_hmm_add_transition(self):
        pass

    def test_hmm_distribution(self):
        pass

    def test_hmm_state(self):
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

    def test_hmm_model(self):
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
    unittest.main()
