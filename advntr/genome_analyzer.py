import logging
import os

from Bio.SeqRecord import SeqRecord

from advntr.profiler import time_usage
from advntr.sam_utils import extract_unmapped_reads_to_fasta_file
from advntr.vntr_finder import VNTRFinder


class GenomeAnalyzer:
    def __init__(self, reference_vntrs, target_vntr_ids, working_directory='./', outfmt='text', is_haploid=False):
        self.reference_vntrs = reference_vntrs
        self.target_vntr_ids = target_vntr_ids
        self.working_dir = working_directory
        self.outfmt = outfmt
        self.is_haploid = is_haploid

        self.vntr_finder = {}
        for ref_vntr in self.reference_vntrs:
            if ref_vntr.id in target_vntr_ids:
                self.vntr_finder[ref_vntr.id] = VNTRFinder(ref_vntr, is_haploid=is_haploid)

    def print_genotype(self, vntr_id, copy_numbers):
        if self.outfmt == 'bed':
            self.print_genotype_in_bed_format(vntr_id, copy_numbers)
        else:
            self.print_genotype_in_text_format(vntr_id, copy_numbers)

    def print_bed_header(self):
        repeats = 'R' if self.is_haploid else 'R1\tR2'
        print('#CHROM\tStart\tEnd\tVNTR_ID\tGene\tMotif\tRefCopy\t%s' % repeats)

    def print_genotype_in_bed_format(self, vntr_id, copy_numbers):
        chr = self.vntr_finder[vntr_id].reference_vntr.chromosome
        start = self.vntr_finder[vntr_id].reference_vntr.start_point
        end = start + self.vntr_finder[vntr_id].reference_vntr.get_length()
        gene = self.vntr_finder[vntr_id].reference_vntr.gene_name
        motif = self.vntr_finder[vntr_id].reference_vntr.pattern
        ref_copy = len(self.vntr_finder[vntr_id].reference_vntr.get_repeat_segments())
        if copy_numbers is None:
            repeats = 'None' if self.is_haploid else 'None\tNone'
        else:
            repeats = str(copy_numbers[0]) if self.is_haploid else '\t'.join([str(cn) for cn in sorted(copy_numbers)])
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chr, start, end, vntr_id, gene, motif, ref_copy, repeats))

    def print_genotype_in_text_format(self, vntr_id, copy_numbers):
        print(vntr_id)
        if copy_numbers is not None:
            if self.is_haploid:
                print(copy_numbers[0])
            else:
                print('/'.join([str(cn) for cn in sorted(copy_numbers)]))
        else:
            print('None')

    @time_usage
    def get_filtered_read_ids(self, read_file, illumina=True):
        vntr_read_ids = {}
        reads = []

        keywords_filename = self.working_dir + '/keywords.txt'
        with open(keywords_filename, 'w') as keywords_file:
            for vid in self.target_vntr_ids:
                keywords = self.vntr_finder[vid].get_keywords_for_filtering(illumina, keyword_size=15)
                keywords_file.write('%s %s\n' % (vid, ' '.join(keywords)))
                vntr_read_ids[vid] = []

        filtering_result = self.working_dir + 'filtering_out.txt'
        os.system('adVNTR-Filtering %s < %s > %s' % (read_file, keywords_filename, filtering_result))
        with open(filtering_result) as infile:
            lines = infile.readlines()
            for line in lines:
                line = line.split()
                if line[0].isdigit() and line[1].isdigit():
                    logging.info('filtering selected %s reads for %s' % (line[1], line[0]))
                    vntr_read_ids[int(line[0])] = set(line[2:])
                else:
                    read = SeqRecord(line[1], line[0])
                    reads.append(read)

        return reads, vntr_read_ids

    @time_usage
    def get_vntr_filtered_reads_map(self, read_file, illumina=True):
        """
        :param read_file: Fasta file containing unmapped reads of sample
        :param illumina: If this parameter is set, it shows that sample contains short reads
        :return: All filtered reads, and a map from VNTR id to read id
        """
        reads, vntr_read_ids = self.get_filtered_read_ids(read_file, illumina)

        return reads, vntr_read_ids

    def find_repeat_counts_from_pacbio_alignment_file(self, alignment_file):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(unmapped_reads_file, False)

        if self.outfmt == 'bed':
            self.print_bed_header()
        for vid in self.target_vntr_ids:
            reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            copy_numbers = self.vntr_finder[vid].find_repeat_count_from_pacbio_alignment_file(alignment_file, reads)
            self.print_genotype(vid, copy_numbers)

    def find_repeat_counts_from_pacbio_reads(self, read_file, naive=False):
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(read_file, False)
        if self.outfmt == 'bed':
            self.print_bed_header()
        for vid in self.target_vntr_ids:
            unmapped_reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            copy_numbers = self.vntr_finder[vid].find_repeat_count_from_pacbio_reads(unmapped_reads, naive)
            self.print_genotype(vid, copy_numbers)

    def find_frameshift_from_alignment_file(self, alignment_file):
        for vid in self.target_vntr_ids:
            result = self.vntr_finder[vid].find_frameshift_from_alignment_file(alignment_file, [])
            print(vid)
            print(result)

    def find_repeat_counts_from_alignment_file(self, alignment_file, average_coverage, update=False):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir)
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(unmapped_reads_file)
        if self.outfmt == 'bed':
            self.print_bed_header()
        for vid in self.target_vntr_ids:
            unmapped_reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            copy_number = self.vntr_finder[vid].find_repeat_count_from_alignment_file(alignment_file, unmapped_reads,
                                                                                      average_coverage, update)
            self.print_genotype(vid, copy_number)

    def find_repeat_counts_from_short_reads(self, read_file):
        if self.outfmt == 'bed':
            self.print_bed_header()
        for vid in self.target_vntr_ids:
            copy_number = self.vntr_finder[vid].find_repeat_count_from_short_reads(read_file)
            self.print_genotype(vid, copy_number)
