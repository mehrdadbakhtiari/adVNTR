import logging
import os

from Bio.SeqRecord import SeqRecord

from advntr.profiler import time_usage
from advntr.sam_utils import extract_unmapped_reads_to_fasta_file
from advntr.vntr_finder import VNTRFinder


class GenomeAnalyzer:
    def __init__(self, ref_vntrs, target_vntr_ids, working_dir='./', outfmt='text', is_haploid=False, ref_filename=None,
                 input_file=None):
        self.reference_vntrs = ref_vntrs
        self.target_vntr_ids = target_vntr_ids
        self.working_dir = working_dir
        self.outfmt = outfmt
        self.is_haploid = is_haploid
        self.ref_filename = ref_filename
        self.input_file = input_file

        self.vntr_finder = {}
        for ref_vntr in self.reference_vntrs:
            if ref_vntr.id in target_vntr_ids:
                self.vntr_finder[ref_vntr.id] = VNTRFinder(ref_vntr, is_haploid, ref_filename)

    def print_genotype(self, vntr_id, genotype_result):
        if self.outfmt == 'bed':
            self.print_genotype_in_bed_format(vntr_id, genotype_result.copy_numbers)
        elif self.outfmt == 'vcf':
            self.print_genotype_in_vcf(vntr_id, genotype_result)
        else:
            self.print_genotype_in_text_format(vntr_id, genotype_result.copy_numbers)

    def print_bed_header(self):
        repeats = 'R' if self.is_haploid else 'R1\tR2'
        print('#CHROM\tStart\tEnd\tVNTR_ID\tGene\tMotif\tRefCopy\t%s' % repeats)

    def print_genotype_in_bed_format(self, vntr_id, copy_numbers):
        chromosome = self.vntr_finder[vntr_id].reference_vntr.chromosome
        start = self.vntr_finder[vntr_id].reference_vntr.start_point
        end = start + self.vntr_finder[vntr_id].reference_vntr.get_length()
        gene = self.vntr_finder[vntr_id].reference_vntr.gene_name
        motif = self.vntr_finder[vntr_id].reference_vntr.pattern
        ref_copy = len(self.vntr_finder[vntr_id].reference_vntr.get_repeat_segments())
        if copy_numbers is None:
            repeats = 'None' if self.is_haploid else 'None\tNone'
        else:
            repeats = str(copy_numbers[0]) if self.is_haploid else '\t'.join([str(cn) for cn in sorted(copy_numbers)])
        print('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chromosome, start, end, vntr_id, gene, motif, ref_copy, repeats))

    def print_vcf_header(self):
        vcf_version = "VCFv4.2"
        from advntr import __version__
        print("##fileformat={}".format(vcf_version))
        print("##source=adVNTR ver. {}".format(__version__))

        print('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">')
        print('##INFO=<ID=VID,Number=1,Type=Integer,Description="VNTR ID">')
        print('##INFO=<ID=RU,Number=1,Type=String,Description="Repeat motif">')
        print('##INFO=<ID=RC,Number=1,Type=Integer,Description="Reference repeat unit count">')

        print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
        print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
        print('##FORMAT=<ID=SR,Number=1,Type=Integer,Description="Spanning read count">')
        print('##FORMAT=<ID=FR,Number=1,Type=Integer,Description="Flanking read count">')
        print('##FORMAT=<ID=ML,Number=1,Type=Float,Description="Maximum likelihood">')

        contigs = set()
        for vntr_id in self.target_vntr_ids:
            ref_vntr = self.vntr_finder[vntr_id].reference_vntr
            chromosome = ref_vntr.chromosome[3:]
            contigs.add(chromosome)
        for contig in sorted(list(contigs)):
            print('##contig=<ID={}>'.format(contig))

        sample = self.input_file.strip().split("/")[-1].split(".")[0]
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample)

    def print_genotype_in_vcf(self, vntr_id, genotype_result):
        vntr = self.vntr_finder[vntr_id].reference_vntr

        # POS
        start = vntr.start_point
        end = start + vntr.get_length()
        # ID
        id = '.'
        # REF
        ref = ''.join(vntr.get_repeat_segments())
        # ALT
        alt = ''
        consensus_seq = vntr.pattern  # TODO use consensus and confidence
        GT = []
        diff_count = 0
        diff_index = -1

        if genotype_result.copy_numbers is None:
            GT.append('.')
            GT.append('.')
        else:
            for index, copy_number in enumerate(genotype_result.copy_numbers):
                if copy_number != vntr.estimated_repeats:
                    diff_index = index
                    diff_count += 1
                    GT.append(diff_count)
                    if len(set(genotype_result.copy_numbers)) == 1:
                        GT.append(diff_count)
                        break
                else:
                    GT.append(0)

        if diff_count == 2:
            alt = consensus_seq * genotype_result.copy_numbers[0] + "," + consensus_seq * genotype_result.copy_numbers[1]
        elif diff_count == 1:
            alt = consensus_seq * genotype_result.copy_numbers[diff_index]
        else:
            alt = '.'

        # QUAL
        qual = '.'
        # FILTER
        filter = '.'

        # INFO
        advntr_vntr_ID = vntr_id  # ID
        reference_repeat_pattern = vntr.pattern  # RU
        reference_repeat_count = vntr.estimated_repeats  # RC

        info_cols = "END=" + str(end) + ";VID=" + str(advntr_vntr_ID) + ";RU=" + reference_repeat_pattern +\
                    ";RC=" + str(reference_repeat_count)

        # FORMAT
        format_cols = "GT:DP:SR:FR:ML"

        format_string = ""
        # GT
        format_string += str(GT[0]) + "/" + str(GT[1]) + ":"  # '/' for unphased, '|' for phased.
        # DP
        format_string += str(genotype_result.recruited_reads_count) + ":"
        # SR
        format_string += str(genotype_result.spanning_reads_count) + ":"
        # FR
        format_string += str(genotype_result.flanking_reads_count) + ":"
        # ML
        format_string += "{0:.4f}".format(genotype_result.maximum_likelihood)

        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(vntr.chromosome, vntr.start_point, id,\
              ref, alt, qual, filter, info_cols, format_cols, format_string))

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

        keywords_filename = self.working_dir + '/keywords_%s.txt' % os.path.basename(read_file)
        with open(keywords_filename, 'w') as keywords_file:
            for vid in self.target_vntr_ids:
                keywords = self.vntr_finder[vid].get_keywords_for_filtering(illumina, keyword_size=15)
                keywords_file.write('%s %s\n' % (vid, ' '.join(keywords)))
                vntr_read_ids[vid] = []

        filtering_result = self.working_dir + 'filtering_out_%s.txt' % os.path.basename(read_file)
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

    def find_repeat_counts_from_pacbio_alignment_file(self, alignment_file, log_pacbio_reads):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir, self.ref_filename)
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(unmapped_reads_file, False)

        if self.outfmt == 'bed':
            self.print_bed_header()
        if self.outfmt == 'vcf':
            self.print_vcf_header()
        for vid in self.target_vntr_ids:
            reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            try:
                genotype_result = self.vntr_finder[vid].find_repeat_count_from_pacbio_alignment_file(alignment_file, reads, log_pacbio_reads)
                self.print_genotype(vid, genotype_result)
            except UnboundLocalError as unbound_local_error:
                error_message = "UnboundLocalError when finding repeat count from pacbio alignment file for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, unbound_local_error)
                logging.warning(error_message)
            except Exception as error:
                error_message = "Error when finding repeat count from pacbio alignment file for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, error)
                logging.warning(error_message)

    def find_repeat_counts_from_pacbio_reads(self, read_file, log_pacbio_reads, naive=False):
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(read_file, False)
        if self.outfmt == 'bed':
            self.print_bed_header()
        if self.outfmt == 'vcf':
            self.print_vcf_header()
        for vid in self.target_vntr_ids:
            unmapped_reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            try:
                copy_numbers = self.vntr_finder[vid].find_repeat_count_from_pacbio_reads(unmapped_reads, log_pacbio_reads, naive)
                self.print_genotype(vid, copy_numbers)
            except UnboundLocalError as unbound_local_error:
                error_message = "UnboundLocalError when finding repeat count from pacbio reads for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, unbound_local_error)
                logging.warning(error_message)
            except Exception as error:
                error_message = "Error when finding repeat count from pacbio reads for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, error)
                logging.warning(error_message)

    def find_frameshift_from_alignment_file(self, alignment_file):
        for vid in self.target_vntr_ids:
            try:
                result = self.vntr_finder[vid].find_frameshift_from_alignment_file(alignment_file, [])
                print(vid)
                print(result)
            except UnboundLocalError as unbound_local_error:
                error_message = "UnboundLocalError when finding frameshift from alignment file for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, unbound_local_error)
                logging.warning(error_message)
            except Exception as error:
                error_message = "Error when finding frameshift from alignment file for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, error)
                logging.warning(error_message)

    def find_repeat_counts_from_alignment_file(self, alignment_file, average_coverage, update=False):
        unmapped_reads_file = extract_unmapped_reads_to_fasta_file(alignment_file, self.working_dir, self.ref_filename)
        filtered_reads, vntr_reads_ids = self.get_vntr_filtered_reads_map(unmapped_reads_file)
        if self.outfmt == 'bed':
            self.print_bed_header()
        if self.outfmt == 'vcf':
            self.print_vcf_header()
        for vid in self.target_vntr_ids:
            unmapped_reads = [read for read in filtered_reads if read.id in vntr_reads_ids[vid]]
            try:
                genotype_result = self.vntr_finder[vid].find_repeat_count_from_alignment_file(alignment_file, unmapped_reads,
                                                                                      average_coverage, update)
                self.print_genotype(vid, genotype_result)
            except UnboundLocalError as unbound_local_error:
                error_message = "UnboundLocalError when finding repeat count from alignment file (illumina) for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, unbound_local_error)
                logging.warning(error_message)
            except Exception as error:
                error_message = "Error when finding repeat count from alignment file (illumina) for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, error)
                logging.warning(error_message)

    def find_repeat_counts_from_short_reads(self, read_file):
        if self.outfmt == 'bed':
            self.print_bed_header()
        if self.outfmt == 'vcf':
            self.print_vcf_header()
        for vid in self.target_vntr_ids:
            try:
                copy_number = self.vntr_finder[vid].find_repeat_count_from_short_reads(read_file)
                self.print_genotype(vid, copy_number)
            except UnboundLocalError as unbound_local_error:
                error_message = "UnboundLocalError when finding repeat count for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, unbound_local_error)
                logging.warning(error_message)
            except Exception as error:
                error_message = "Error when finding repeat count for vntr id {}: {}. Skipping genotyping for this VNTR.".format(vid, error)
                logging.warning(error_message)
