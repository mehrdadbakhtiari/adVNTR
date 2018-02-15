import unittest
from sam_utils import get_related_reads_and_read_count_in_samfile


class TestBlastFilteringMethod(unittest.TestCase):

    def test_full_sensitivity(self):
        self.assertEqual('foo'.upper(), 'FOO')
        blast_ids = self.filter_reads_with_keyword_matching()

        reference_end_pos = self.reference_vntr.start_point + self.reference_vntr.get_length()
        samfile = 'original_reads/paired_dat.sam'
        related_reads, read_count = get_related_reads_and_read_count_in_samfile(self.reference_vntr.pattern,
                                                                                self.reference_vntr.start_point,
                                                                                read_file=samfile,
                                                                                pattern_end=reference_end_pos)
        for re_read in related_reads:
            if re_read not in blast_ids:
                print('FN in filtering')


if __name__ == '__main__':
    unittest.main()
