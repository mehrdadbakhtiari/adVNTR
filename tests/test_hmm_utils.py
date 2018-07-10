import json
import unittest

from advntr import hmm_utils


class TestHMMUtils(unittest.TestCase):

    def setUp(self):
        with open("tests/data/hmm_utils.json", "r") as read_file:
            self.data = json.load(read_file)

    def test_extract_repeating_segments_from_read(self):
        visited_states = self.data['visited_states'].split(',')
        repeats, states = hmm_utils.extract_repeating_segments_from_read(self.data['sequence'], visited_states)
        self.assertEqual(self.data['correct_repeats'], repeats)

    def test_multiple_alignment_for_real_data(self):
        visited_states = self.data['visited_states'].split(',')
        repeats, states = hmm_utils.extract_repeating_segments_from_read(self.data['sequence'], visited_states)
        alignment = hmm_utils.get_multiple_alignment_of_viterbi_paths(repeats, states)
        self.assertEqual(self.data['alignment'], alignment)

    def test_multiple_alignment_for_two_sequences(self):
        repeats = ['ACTTA', 'ATTGA']
        states = [['M1', 'M2', 'M3', 'M4', 'M5'],
                  ['M1', 'D2', 'M3', 'M4', 'I4', 'M5']]
        correct_alignment = ['ACTT-A', 'A-TTGA']
        alignment = hmm_utils.get_multiple_alignment_of_viterbi_paths(repeats, states)
        self.assertEqual(alignment, correct_alignment)
