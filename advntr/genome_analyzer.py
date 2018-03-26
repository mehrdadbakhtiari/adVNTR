from Bio import SeqIO

from advntr.profiler import time_usage
from advntr.sam_utils import extract_unmapped_reads_to_fasta_file
from advntr.vntr_finder import VNTRFinder


class GenomeAnalyzer:
    def __init__(self, reference_vntrs, target_vntr_ids, working_directory='./'):
        self.reference_vntrs = reference_vntrs
        self.target_vntr_ids = target_vntr_ids
        self.working_dir = working_directory

        self.vntr_finder = {}
        for ref_vntr in self.reference_vntrs:
            if ref_vntr.id in target_vntr_ids:
                self.vntr_finder[ref_vntr.id] = VNTRFinder(ref_vntr)

    @staticmethod
    def print_genotype(vntr_id, copy_numbers):
        print(vntr_id)
        if copy_numbers is not None:
            print('/'.join([str(cn) for cn in sorted(copy_numbers)]))
        else:
            print('None')

    @time_usage
    def get_vntr_filtered_reads_map(self, read_file, illumina=True):
        vntr_reads = {}
        vntr_read_ids = {}
        empty_set = True
        for vid in self.target_vntr_ids:
            vntr_reads[vid] = []
            read_ids = self.vntr_finder[vid].filter_reads_with_keyword_matching(self.working_dir, read_file, illumina)
            vntr_read_ids[vid] = read_ids
            if len(read_ids) > 0:
                empty_set = False

        if not empty_set:
            unmapped_reads = SeqIO.parse(read_file, 'fasta')
            for read in unmapped_reads:
                for vntr_id in vntr_read_ids.keys():
                    if read.id in vntr_read_ids[vntr_id]:
                        vntr_reads[vntr_id].append(read)
        return vntr_reads

    def find_repeat_counts_from_pacbio_alignment_file(self, alignment_file):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        vntr_reads = self.get_vntr_filtered_reads_map(unmapped_reads_file, False)

        for vid in self.target_vntr_ids:
            reads = vntr_reads[vid]
            copy_numbers = self.vntr_finder[vid].find_repeat_count_from_pacbio_alignment_file(alignment_file, reads)
            self.print_genotype(vid, copy_numbers)

    def find_repeat_counts_from_pacbio_reads(self, read_file, naive=False):
        vntr_reads = self.get_vntr_filtered_reads_map(read_file, False)
        for vid in self.target_vntr_ids:
            copy_numbers = self.vntr_finder[vid].find_repeat_count_from_pacbio_reads(vntr_reads[vid], naive)
            self.print_genotype(vid, copy_numbers)

    def find_frameshift_from_alignment_file(self, alignment_file):
        for vid in self.target_vntr_ids:
            result = self.vntr_finder[vid].find_frameshift_from_alignment_file(alignment_file, [])
            print(vid)
            print(result)

    def find_repeat_counts_from_alignment_file(self, alignment_file, average_coverage):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        vntr_reads = self.get_vntr_filtered_reads_map(unmapped_reads_file)
        for vid in self.target_vntr_ids:
            unmapped_reads = vntr_reads[vid]
            copy_number = self.vntr_finder[vid].find_repeat_count_from_alignment_file(alignment_file, unmapped_reads,
                                                                                      average_coverage)
            self.print_genotype(vid, copy_number)

    def find_repeat_counts_from_short_reads(self, read_file):
        for vid in self.target_vntr_ids:
            copy_number = self.vntr_finder[vid].find_repeat_count_from_short_reads(read_file)
            self.print_genotype(vid, copy_number)
