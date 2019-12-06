from collections import defaultdict
import sys
sys.path.append('../')

from hmm import State
from hmm import Model
from hmm import DiscreteDistribution

from advntr import profile_hmm

import unittest

import numpy as np

class TestMethods(unittest.TestCase):

    def setUp(self):
        pass


    def test_example_pomegranate(self):
        """
        This example is taken from https://pomegranate.readthedocs.io/en/latest/HiddenMarkovModel.html
        """

        from pomegranate import DiscreteDistribution, State, HiddenMarkovModel
        d1 = DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        d2 = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        d3 = DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})

        s1 = State(d1, name="s1")
        s2 = State(d2, name="s2")
        s3 = State(d3, name="s3")

        model = HiddenMarkovModel(name='example')
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

        answer = model.log_probability(list('ACGACTATTCGAT'))
        expected = -22.73896159971087

        print(" > log probability of pomegranate model for 'ACGACTATTCGAT': ", answer)
        # should be -22.73896159971087
        self.assertAlmostEqual(expected, answer)

        print(" > pomegranate viterbi states: ",", ".join(state.name for i, state in model.viterbi(list('ACGACTATTCGAT'))[1]))
        # should be example - start, s1, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s3, example - end

        answer = ", ".join(state.name for i, state in model.viterbi(list('ACGACTATTCGAT'))[1])
        expected = "example-start, s1, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s3, example-end"
        self.assertEqual(expected, answer)

    def test_log_probability(self):
        """
        This test is for computing probability of a sequence with the model
        """

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

        answer = model.log_probability(list('ACGACTATTCGAT'))
        expected = -22.73896159971087

        print(" > log probability of hmm model for 'ACGACTATTCGAT': ", answer)
        # should be -22.73896159971087
        self.assertAlmostEqual(expected, answer)

    def test_dense_transition_matrix(self):
        """
        This test checks the computation of the dense_transition_matrix method
        """
        from pomegranate import DiscreteDistribution as pome_DiscreteDistribution 
        from pomegranate import State as pome_State 
        from pomegranate import HiddenMarkovModel as pome_HiddenMarkovModel
        d1 = pome_DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        d2 = pome_DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        d3 = pome_DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})

        s1 = pome_State(d1, name="s1")
        s2 = pome_State(d2, name="s2")
        s3 = pome_State(d3, name="s3")

        pome_model = pome_HiddenMarkovModel(name='example')
        pome_model.add_states([s1, s2, s3])
        pome_model.add_transition(pome_model.start, s1, 0.90)
        pome_model.add_transition(pome_model.start, s2, 0.10)
        pome_model.add_transition(s1, s1, 0.80)
        pome_model.add_transition(s1, s2, 0.20)
        pome_model.add_transition(s2, s2, 0.90)
        pome_model.add_transition(s2, s3, 0.10)
        pome_model.add_transition(s3, s3, 0.70)
        pome_model.add_transition(s3, pome_model.end, 0.30)
        pome_model.bake()

        pome_mat = pome_model.dense_transition_matrix()
        #print(pome_mat)

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

        mat = model.dense_transition_matrix()
        #print(mat)

        self.assertAlmostEqual(np.sum(np.abs(mat-pome_mat)), 0.0)

    def test_from_matrix(self):

        d1 = DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        d2 = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        d3 = DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})
        distributions = [d1, d2, d3]
    
        transition_probabilities = np.array([[0.8, 0.2, 0.0],
                                             [0.0, 0.9, 0.1],
                                             [0.0, 0.0, 0.7]])

        starts = np.array([0.9, 0.1, 0.0])
        ends = np.array([0.0, 0.0, 0.3])
        state_names = ['s1', 's2', 's3']

        model = Model.from_matrix(transition_probabilities, distributions, starts, ends,
        state_names, name="from_matrix_model", verbose=False, merge=None )

        answer = model.log_probability(list('ACGACTATTCGAT'))
        expected = -22.73896159971087

        print(" > log probability of model from matrix for 'ACGACTATTCGAT': ", answer)
        # should be -22.73896159971087
        self.assertAlmostEqual(expected, answer)
        pass

    def test_concatenate_without_delete_state(self):

        #from pomegranate import DiscreteDistribution as pome_DiscreteDistribution 
        #from pomegranate import State as pome_State 
        #from pomegranate import HiddenMarkovModel as pome_HiddenMarkovModel
        #d1 = pome_DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        #d2 = pome_DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        #d3 = pome_DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})

        #s1_1 = pome_State(d1, name="s1_1")
        #s2_1 = pome_State(d2, name="s2_1")
        #s3_1 = pome_State(d3, name="s3_1")

        #s1_2 = pome_State(d3, name="s1_2")
        #s2_2 = pome_State(d1, name="s2_2")
        #s3_2 = pome_State(d2, name="s3_2")

        #pome_model_1 = pome_HiddenMarkovModel(name='model_1')
        #pome_model_1.add_states([s1_1, s2_1, s3_1])
        #pome_model_1.add_transition(pome_model_1.start, s1_1, 0.90)
        #pome_model_1.add_transition(pome_model_1.start, s2_1, 0.10)
        #pome_model_1.add_transition(s1_1, s1_1, 0.80)
        #pome_model_1.add_transition(s1_1, s2_1, 0.20)
        #pome_model_1.add_transition(s2_1, s2_1, 0.90)
        #pome_model_1.add_transition(s2_1, s3_1, 0.10)
        #pome_model_1.add_transition(s3_1, s3_1, 0.70)
        #pome_model_1.add_transition(s3_1, pome_model_1.end, 0.30)
        #pome_model_1.bake()

        #pome_model_2 = pome_HiddenMarkovModel(name='model_2')
        #pome_model_2.add_states([s1_2, s2_2, s3_2])
        #pome_model_2.add_transition(pome_model_2.start, s1_2, 0.90)
        #pome_model_2.add_transition(pome_model_2.start, s2_2, 0.10)
        #pome_model_2.add_transition(s1_2, s1_2, 0.80)
        #pome_model_2.add_transition(s1_2, s2_2, 0.20)
        #pome_model_2.add_transition(s2_2, s2_2, 0.90)
        #pome_model_2.add_transition(s2_2, s3_2, 0.10)
        #pome_model_2.add_transition(s3_2, s3_2, 0.70)
        #pome_model_2.add_transition(s3_2, pome_model_2.end, 0.30)
        #pome_model_2.bake()

        #pome_model_1.concatenate(pome_model_2)
        #pome_model_1.bake()

        #pome_log_prob, pome_vpath = pome_model_1.viterbi(list('ACGACTATTCGAT'))
        #print(" > Pomegranate model viterbi states: ", " ".join(state.name for i, state in pome_vpath))
        #print(" > Pomegranate model log probability: ", pome_log_prob)
        ##expected = "model_1-start s1_1 s2_1 s2_1 s2_1 s2_1 s2_1 s2_1 s2_1 s2_1 s3_1 model_2-start s1_2 s2_2 s3_2 model_2-end"

        d1 = DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        d2 = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        d3 = DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})

        s1_1 = State(d1, name="s1_1")
        s2_1 = State(d2, name="s2_1")
        s3_1 = State(d3, name="s3_1")

        s1_2 = State(d3, name="s1_2")
        s2_2 = State(d1, name="s2_2")
        s3_2 = State(d2, name="s3_2")

        model_1 = Model(name='model_1')
        model_1.add_states([s1_1, s2_1, s3_1])
        model_1.add_transition(model_1.start, s1_1, 0.90)
        model_1.add_transition(model_1.start, s2_1, 0.10)
        model_1.add_transition(s1_1, s1_1, 0.80)
        model_1.add_transition(s1_1, s2_1, 0.20)
        model_1.add_transition(s2_1, s2_1, 0.90)
        model_1.add_transition(s2_1, s3_1, 0.10)
        model_1.add_transition(s3_1, s3_1, 0.70)
        model_1.add_transition(s3_1, model_1.end, 0.30)
        model_1.bake()

        model_2 = Model(name='model_2')
        model_2.add_states([s1_2, s2_2, s3_2])
        model_2.add_transition(model_2.start, s1_2, 0.90)
        model_2.add_transition(model_2.start, s2_2, 0.10)
        model_2.add_transition(s1_2, s1_2, 0.80)
        model_2.add_transition(s1_2, s2_2, 0.20)
        model_2.add_transition(s2_2, s2_2, 0.90)
        model_2.add_transition(s2_2, s3_2, 0.10)
        model_2.add_transition(s3_2, s3_2, 0.70)
        model_2.add_transition(s3_2, model_2.end, 0.30)
        model_2.bake()

        model_1.concatenate(model_2)
        model_1.bake()

        log_prob, vpath = model_1.viterbi(list('ACGACTATTCGAT'))
        answer_vpath = " ".join(state.name for i, state in vpath)
        print(" > Our model viterbi states: ", answer_vpath)
        print(" > Our model log probability: ", log_prob)
        expected_vpath = "model_1-start s1_1 s2_1 s2_1 s2_1 s2_1 s2_1 s2_1 s2_1 s2_1 s3_1 model_1-end model_2-start s1_2 s2_2 s3_2 model_2-end"
        expected_logProb = -27.589111223253287

        self.assertEqual(expected_vpath, answer_vpath)
        self.assertAlmostEqual(expected_logProb, log_prob)
        
    #def test_hmm_add_transition_before_add_state(self):

    #    hmm = Model("add transtion test")
    #    emission = {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2}
    #    distribution = DiscreteDistribution(emission)
    #    s1 = State(distribution, name="s1")
    #    s2 = State(distribution, name="s2")
    #    transition = dict({})
    #    transition['s1'] = dict([('s2', 0.3)])

    #    with self.assertRaises(Exception) as context:
    #        hmm.add_transition(s1, s2, transition['s1']['s2'])  # Will raise an exception

    #    self.assertTrue("No such state" in context.exception)

    def test_hmm_add_transition(self):
        pass
        model = Model("transition test")
        distribution = {'A': 0.2, 'C': 0.3, 'G': 0.3, 'T': 0.2}
        s1 = State(distribution, name="s1")
        s2 = State(distribution, name="s2")

        # The model should have states first
        model.add_states(*[s1, s2])

        # Define transitions
        transition = dict()
        transition['s1'] = defaultdict(lambda: 0, type=float)
        transition['s2'] = defaultdict(lambda: 0, type=float)

        transition['s1']['s2'] = 0.6
        transition['s1']['s1'] = 0.4
        transition['s2']['s1'] = 0.6
        transition['s2']['s2'] = 0.3

        # Add transitions
        model.add_transition(model.start, s1, 1)
        model.add_transition(s1, s2, transition['s1']['s2'])
        model.add_transition(s1, s1, transition['s1']['s1'])
        model.add_transition(s2, s1, transition['s2']['s1'])
        model.add_transition(s2, s2, transition['s2']['s2'])
        model.add_transition(s2, s2, transition['s2']['s2'])
        model.add_transition(s2, model.end, 0.1)

    def test_hmm_distribution(self):
        pass

    def test_hmm_state(self):
        pass
        transitions = defaultdict(lambda: defaultdict(float))
        emissions = defaultdict(lambda: defaultdict(float))

        # State(insert_distribution, name='I%s_%s' % (i, repeat))
        distribution = {'A':0.2, 'C':0.3, 'G':0.3, 'T':0.2}
        state1 = State(distribution, name="intron")
        state2 = State(distribution, name="exon")

        emissions['start']['intron'] = 0.5
        emissions['start']['exon'] = 0.5
        emissions['intron']['intron'] = 0.6
        emissions['intron']['exon'] = 0.4
        emissions['exon']['intron'] = 0.6
        emissions['exon']['exon'] = 0.4

        self.assertEqual("intron", state1.name)
        self.assertEqual("exon", state2.name)

    def test_ordering_states(self):

        d1 = DiscreteDistribution({'A': 0.35, 'C': 0.20, 'G': 0.05, 'T': 0.40})
        d2 = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
        d3 = DiscreteDistribution({'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10})

        s0 = State(d1, name="I0_1")
        s1 = State(d1, name="I1_1")
        s2 = State(d2, name="M1_1")
        s3 = State(d3, name="D1_1")
        s4 = State(d1, name="I2_1")
        s5 = State(d2, name="M2_1")
        s6 = State(d3, name="D2_1")

        model = Model(name='ordering')
        model.add_states([s0, s1, s2, s3, s4, s5, s6])
        model.bake()

        names = [state.name for state in model.states]
        answer = " ".join(names)
        expected = "ordering-start I0_1 D1_1 M1_1 I1_1 D2_1 M2_1 I2_1 ordering-end"  # Only in our hmm model
        self.assertEqual(expected, answer)

    def test_viterbi_path_without_delete_state(self):

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

        print(" > Our model viterbi states: ", ", ".join(state.name for i, state in model.viterbi(list('ACGACTATTCGAT'))[1]))
        # should be example-start, s1, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s3, example-end
        expected = "example-start, s1, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s2, s3, example-end"
        answer = (", ".join(state.name for i, state in model.viterbi(list('ACGACTATTCGAT'))[1]))
        self.assertEqual(expected, answer)


if __name__ == "__main__":
    unittest.main()
