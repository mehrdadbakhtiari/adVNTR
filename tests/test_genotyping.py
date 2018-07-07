import unittest

from advntr.reference_vntr import ReferenceVNTR
from advntr.vntr_finder import VNTRFinder


class TestGenotyping(unittest.TestCase):

    def get_reference_vntr(self):
        ref_vntr = ReferenceVNTR(1, 'CACA', 1000, 'chr1', None, None)
        return ref_vntr

    def test_statistical_model_for_haploid_case(self):
        vntr_finder = VNTRFinder(self.get_reference_vntr())
        genotype = vntr_finder.find_genotype_based_on_observed_repeats([3, 3, 3, 3, 3])
        self.assertEqual(genotype, (3, 3))

    def test_statistical_model_for_diploid_case(self):
        vntr_finder = VNTRFinder(self.get_reference_vntr())
        genotype = vntr_finder.find_genotype_based_on_observed_repeats([2, 2, 3, 3, 3])
        if genotype[0] > genotype[1]:
            genotype = (genotype[1], genotype[0])
        self.assertEqual(genotype, (2, 3))

    def test_statistical_model_for_erroneous_diploid_case(self):
        vntr_finder = VNTRFinder(self.get_reference_vntr())
        genotype = vntr_finder.find_genotype_based_on_observed_repeats([4, 5, 5, 5, 7, 8, 8, 8, 9])
        if genotype[0] > genotype[1]:
            genotype = (genotype[1], genotype[0])
        self.assertEqual(genotype, (5, 8))

    def test_recruit_read_for_positive_read(self):
        vntr_finder = VNTRFinder(self.get_reference_vntr())
        logp = -20
        vpath = []
        min_score_to_count_read = -50
        read_length = 100
        results = vntr_finder.recruit_read(logp, vpath, min_score_to_count_read, read_length)
        self.assertEqual(results, True)
