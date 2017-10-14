from Bio import SeqIO

from profiler import time_usage
from sam_utils import extract_unmapped_reads_to_fasta_file
from vntr_finder import VNTRFinder


class GenomeAnalyzer:
    def __init__(self, reference_vntrs, target_vntr_ids, working_directory='./'):
        self.reference_vntrs = reference_vntrs
        self.target_vntr_ids = target_vntr_ids
        self.working_dir = working_directory

        self.vntr_finder = {}
        for vntr_id in self.target_vntr_ids:
            self.vntr_finder[vntr_id] = VNTRFinder(self.reference_vntrs[vntr_id])

    @time_usage
    def get_vntr_filtered_reads_map(self, read_file, illumina=True):
        vntr_reads = {}
        vntr_read_ids = {}
        for id in self.target_vntr_ids:
            vntr_reads[id] = []
            read_ids = self.vntr_finder[id].filter_reads_with_keyword_matching(self.working_dir, read_file, illumina)
            vntr_read_ids[id] = read_ids

        unmapped_reads = SeqIO.parse(read_file, 'fasta')
        for read in unmapped_reads:
            for vntr_id in vntr_read_ids.keys():
                if read.id in vntr_read_ids[vntr_id]:
                    vntr_reads[vntr_id].append(read)
        return vntr_reads

    def find_repeat_counts_from_pacbio_alignment_file(self, alignment_file):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        vntr_reads = self.get_vntr_filtered_reads_map(unmapped_reads_file, False)

        for id in self.target_vntr_ids:
            reads = vntr_reads[id]
            copy_numbers = self.vntr_finder[id].find_repeat_count_from_pacbio_alignment_file(alignment_file, reads)
            print(id)
            print(copy_numbers)

    def find_repeat_counts_from_pacbio_reads(self, read_file):
        vntr_reads = self.get_vntr_filtered_reads_map(read_file, False)
        for id in self.target_vntr_ids:
            copy_numbers = self.vntr_finder[id].find_repeat_count_from_pacbio_reads(vntr_reads[id])
            print(id)
            print(copy_numbers)

    def find_frameshift_from_alignment_file(self, alignment_file):
        for id in self.target_vntr_ids:
            result = self.vntr_finder[id].find_frameshift_from_alignment_file(alignment_file, [])
            print(id)
            print(result)

    def find_repeat_counts_from_alignment_file(self, alignment_file):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        vntr_reads = self.get_vntr_filtered_reads_map(unmapped_reads_file)
        for id in self.target_vntr_ids:
            unmapped_reads = vntr_reads[id]
            copy_number = self.vntr_finder[id].find_repeat_count_from_alignment_file(alignment_file, unmapped_reads)
            print(id)
            print(copy_number)

    def find_repeat_counts_from_short_reads(self, read_file):
        for id in self.target_vntr_ids:
            copy_number = self.vntr_finder[id].find_repeat_count_from_short_reads(read_file)
            print(id)
            print(copy_number)
