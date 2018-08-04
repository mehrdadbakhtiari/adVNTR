import unittest

from advntr.reference_vntr import ReferenceVNTR
from advntr.vntr_finder import VNTRFinder


class TestFrameshiftIdentification(unittest.TestCase):

    def get_reference_vntr(self, ru_count=10):
        pattern = 'ACGTACGT'
        ref_vntr = ReferenceVNTR(1, pattern, 1000, 'chr1', None, None)
        ref_vntr.repeat_segments = [pattern] * ru_count
        return ref_vntr

    def get_vntr_finder(self):
        vntr_finder = VNTRFinder(self.get_reference_vntr())
        return vntr_finder

    def test_frameshift_in_uniform_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 14
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, True)

    def test_frameshift_with_high_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 18
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, True)

    def test_frameshift_with_low_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 7
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, True)

    def test_frameshift_with_extremely_low_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 3
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, True)

    def test_normal_vntr_with_high_error_in_uniform_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 2
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, False)

    def test_normal_vntr_with_low_error_in_uniform_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 1
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, False)

    def test_normal_vntr_without_error_in_uniform_coverage(self):
        vntr_finder = self.get_vntr_finder()

        avg_bp_coverage = 14.0
        observed_indels = 0
        expected_indel_transitions = 1 / avg_bp_coverage
        frameshift = vntr_finder.identify_frameshift(avg_bp_coverage, observed_indels, expected_indel_transitions)
        self.assertEqual(frameshift, False)
