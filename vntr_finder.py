import os
from multiprocessing import Process, Manager, Value, Semaphore
from random import random
from uuid import uuid4

import numpy
import pysam
from Bio import SeqIO, pairwise2

import settings
from blast_wrapper import get_blast_matched_ids, make_blast_database
from coverage_bias import CoverageBiasDetector, CoverageCorrector
from hmm_utils import *
from sam_utils import get_related_reads_and_read_count_in_samfile, extract_unmapped_reads_to_fasta_file


class VNTRFinder:
    """Find the VNTR structure of a reference VNTR in NGS data of the donor."""

    def __init__(self, reference_vntr):
        self.reference_vntr = reference_vntr
        self.min_repeat_bp_to_add_read = 2
        if len(self.reference_vntr.pattern) < 30:
            self.min_repeat_bp_to_add_read = 2
        self.min_repeat_bp_to_count_repeats = 2

    def build_vntr_matcher_hmm(self, copies, read_length=150):
        patterns = self.reference_vntr.get_repeat_segments() * 100
        flanking_region_size = read_length - 10
        left_flanking_region = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking_region = self.reference_vntr.right_flanking_region[:flanking_region_size]

        vntr_matcher = get_suffix_matcher_hmm(left_flanking_region)
        right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
        repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
        vntr_matcher.concatenate(repeats_matcher)
        vntr_matcher.concatenate(right_flanking_matcher)
        vntr_matcher.bake(merge=None)
        return vntr_matcher

    def get_vntr_matcher_hmm(self, read_length):
        """Try to load trained HMM for this VNTR
        If there was no trained HMM, it will build one and store it for later usage
        """
        copies = int(round(float(read_length) / len(self.reference_vntr.pattern) + 0.5))

        base_name = str(self.reference_vntr.id) + '_' + str(read_length) + '.json'
        stored_hmm_file = settings.TRAINED_HMMS_DIR + base_name
        if settings.USE_TRAINED_HMMS and os.path.isfile(stored_hmm_file):
            model = Model()
            model = model.from_json(stored_hmm_file)
            return model

        vntr_matcher = self.build_vntr_matcher_hmm(copies, read_length)

        json_str = vntr_matcher.to_json()
        with open(stored_hmm_file, 'w') as outfile:
            outfile.write(json_str)
        return vntr_matcher

    def filter_reads_with_keyword_matching(self, working_directory, read_file, short_reads=True):
        db_name = working_directory[:-1]
        blast_db_name = working_directory + db_name
        empty_db = False
        if not os.path.exists(blast_db_name + '.nsq') and not os.path.exists(blast_db_name + '.nal'):
            empty_db = make_blast_database(read_file, blast_db_name)

        word_size = int(len(self.reference_vntr.pattern)/3)
        if word_size > 11:
            word_size = 11
        word_size = str(word_size)

        blast_ids = set([])
        search_id = str(uuid4()) + str(self.reference_vntr.id)
        queries = self.reference_vntr.get_repeat_segments()
        identity_cutoff = '40'
        if not short_reads:
            queries = [self.reference_vntr.left_flanking_region[-80:], self.reference_vntr.right_flanking_region[:80]]
            word_size = str('10')
            identity_cutoff = '70'
        if not empty_db:
            for query in queries:
                search_result = get_blast_matched_ids(query, blast_db_name, max_seq='50000', word_size=word_size,
                                                      evalue=10, search_id=search_id, identity_cutoff=identity_cutoff)
                blast_ids |= search_result

        print('blast selected ', len(blast_ids), ' reads')
        if len(blast_ids) == len(self.reference_vntr.get_repeat_segments()) * 50 * 1000:
            with open('errors.txt', 'a') as out:
                out.write('maximum number of read selected in filtering for pattern %s\n' % self.reference_vntr.id)
        return blast_ids

    @staticmethod
    def add_hmm_score_to_list(sema, hmm, read, result_scores):
        logp, vpath = hmm.viterbi(str(read.seq))
        result_scores.append(logp)
        sema.release()

    def add_false_read_scores_of_chromosome(self, samfile, reference, hmm, false_scores):
        process_list = []
        sema = Semaphore(5)
        for read in samfile.fetch(reference, multiple_iterators=True):
            if read.is_unmapped:
                continue
            if random() > settings.SCORE_FINDING_READS_FRACTION:
                continue
            read_start = read.reference_start
            read_end = read.reference_end if read.reference_end else read_start + len(read.seq)
            vntr_start = self.reference_vntr.start_point
            vntr_end = vntr_start + self.reference_vntr.get_length()
            reference_name = read.reference_name
            if not reference_name.startswith('chr'):
                reference_name = 'chr' + reference_name
            if reference_name == self.reference_vntr.chromosome and (
                        vntr_start <= read_start < vntr_end or vntr_start < read_end <= vntr_end):
                continue
            if read.seq.count('N') > 0:
                continue
            sema.acquire()
            p = Process(target=VNTRFinder.add_hmm_score_to_list, args=(sema, hmm, read, false_scores))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

    def calculate_min_score_to_select_a_read(self, hmm, alignment_file):
        """Calculate the score distribution of false positive reads
        and return score to select the 0.0001 percentile of the distribution
        """
        process_list = []
        manager = Manager()
        false_scores = manager.list()
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        refs = [ref for ref in samfile.references if ref in settings.CHROMOSOMES or 'chr'+ref in settings.CHROMOSOMES]
        for ref in refs:
            p = Process(target=self.add_false_read_scores_of_chromosome, args=(samfile, ref, hmm, false_scores))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        score = numpy.percentile(false_scores, 100 - 0.0001)
        return score

    def get_min_score_to_select_a_read(self, hmm, alignment_file, read_length):
        """Try to load the minimum score for this VNTR

        If the score is not stored, it will compute the score and write it for this VNTR in precomputed data.
        """
        base_name = str(self.reference_vntr.id) + '_' + str(read_length) + '.scores'
        stored_scores_file = settings.TRAINED_HMMS_DIR + base_name
        if settings.USE_TRAINED_HMMS and os.path.isfile(stored_scores_file):
            with open(stored_scores_file, 'r') as infile:
                frac_score = [(line.split()[0], line.split()[1]) for line in infile.readlines() if line.strip() != '']
                fraction_score_map = {float(reads_fraction): float(score) for reads_fraction, score in frac_score}
            if settings.SCORE_FINDING_READS_FRACTION in fraction_score_map.keys():
                return fraction_score_map[settings.SCORE_FINDING_READS_FRACTION]

        print('Minimum score is not precomputed for vntr id: %s' % self.reference_vntr.id)
        score = self.calculate_min_score_to_select_a_read(hmm, alignment_file)
        print('computed score: ', score)
        with open(stored_scores_file, 'a') as outfile:
            outfile.write('%s %s\n' % (settings.SCORE_FINDING_READS_FRACTION, score))

        return score

    def process_unmapped_read(self, sema, read_segment, hmm, min_score_to_count_read,
                              vntr_bp_in_unmapped_reads, selected_reads):
        if read_segment.seq.count('N') <= 0:
            logp, vpath = hmm.viterbi(str(read_segment.seq))
            rev_logp, rev_vpath = hmm.viterbi(str(read_segment.seq.reverse_complement()))
            if logp < rev_logp:
                logp = rev_logp
                vpath = rev_vpath
            repeat_bps = get_number_of_repeat_bp_matches_in_vpath(vpath)
            if logp > min_score_to_count_read:
                if repeat_bps > self.min_repeat_bp_to_count_repeats:
                    vntr_bp_in_unmapped_reads.value += repeat_bps
                if repeat_bps > self.min_repeat_bp_to_add_read:
                    selected_reads.append(read_segment.id)
        sema.release()

    def check_if_read_spans_vntr(self, read, length_distribution):
        flanking_region_size = 100
        left_flanking = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking = self.reference_vntr.right_flanking_region[:flanking_region_size]
        left_align = pairwise2.align.localms(str(read.seq), left_flanking, 1, -1, -1, -1)[0]
        if left_align[2] < len(left_flanking) * 0.8:
            return
        right_align = pairwise2.align.localms(str(read.seq), right_flanking, 1, -1, -1, -1)[0]
        if right_align[2] < len(right_flanking) * 0.8:
            return
        if right_align[3] < left_align[3]:
            return
        print(str(read.seq)[left_align[3]:right_align[3]+flanking_region_size])
        length_distribution.append(right_align[3] - (left_align[3] + flanking_region_size))

    def find_repeat_count_from_pacbio_alignment_file(self, alignment_file, working_directory='./'):
        length_distribution = []

        unmapped_read_file = extract_unmapped_reads_to_fasta_file(alignment_file, working_directory)
        print('unmapped reads extracted')

        filtered_read_ids = self.filter_reads_with_keyword_matching(working_directory, unmapped_read_file, False)
        print('unmapped reads filtered')

        unmapped_reads = SeqIO.parse(unmapped_read_file, 'fasta')
        for read in unmapped_reads:
            if read.id in filtered_read_ids:
                self.check_if_read_spans_vntr(read, length_distribution)

        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        region_start = vntr_start - 30000
        region_end = vntr_end + 30000
        chromosome = self.reference_vntr.chromosome[3:]
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        for read in samfile.fetch(chromosome, region_start, region_end):
            self.check_if_read_spans_vntr(read, length_distribution)

        print('length_distribution: ', length_distribution)
        # max_copies = max(length_distribution) / float(len(self.reference_vntr.pattern))
        # vntr_matcher = self.build_vntr_matcher_hmm(max_copies)

        average_length = sum(length_distribution) / float(len(length_distribution)) if len(length_distribution) else 0
        copy_count = average_length / float(len(self.reference_vntr.pattern))
        print('copy_count: ', copy_count)
        return copy_count

    def find_repeat_count_from_alignment_file(self, alignment_file, working_directory='./'):
        unmapped_read_file = extract_unmapped_reads_to_fasta_file(alignment_file, working_directory)
        print('unmapped reads extracted')

        filtered_read_ids = self.filter_reads_with_keyword_matching(working_directory, unmapped_read_file)
        print('unmapped reads filtered')

        hmm = None
        min_score_to_count_read = None
        sema = Semaphore(settings.CORES)
        manager = Manager()
        selected_reads = manager.list()
        vntr_bp_in_unmapped_reads = Value('d', 0.0)

        number_of_reads = 0
        read_length = 150

        process_list = []

        unmapped_reads = SeqIO.parse(unmapped_read_file, 'fasta')
        for read_segment in unmapped_reads:
            if number_of_reads == 0:
                read_length = len(str(read_segment.seq))
            number_of_reads += 1
            if not hmm:
                hmm = self.get_vntr_matcher_hmm(read_length=read_length)
                min_score_to_count_read = self.get_min_score_to_select_a_read(hmm, alignment_file, read_length)

            if read_segment.id in filtered_read_ids:
                sema.acquire()
                p = Process(target=self.process_unmapped_read, args=(sema, read_segment, hmm, min_score_to_count_read,
                                                                     vntr_bp_in_unmapped_reads, selected_reads))
                process_list.append(p)
                p.start()
        for p in process_list:
            p.join()

        print('vntr base pairs in unmapped reads:', vntr_bp_in_unmapped_reads.value)

        vntr_bp_in_mapped_reads = 0
        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        region_start = vntr_start - settings.MAX_INSERT_SIZE
        region_end = vntr_end + settings.MAX_INSERT_SIZE
        chromosome = self.reference_vntr.chromosome[3:]
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        for read in samfile.fetch(chromosome, region_start, region_end):
            if read.is_unmapped:
                continue
            read_end = read.reference_end if read.reference_end else read.reference_start + len(read.seq)
            if vntr_start <= read.reference_start < vntr_end or vntr_start < read_end <= vntr_end:
                end = min(read_end, vntr_end)
                start = max(read.reference_start, vntr_start)
                vntr_bp_in_mapped_reads += end - start
        print('vntr base pairs in mapped reads:', vntr_bp_in_mapped_reads)

        total_counted_vntr_bp = vntr_bp_in_unmapped_reads.value + vntr_bp_in_mapped_reads
        pattern_occurrences = total_counted_vntr_bp / float(len(self.reference_vntr.pattern))
        bias_detector = CoverageBiasDetector(alignment_file, self.reference_vntr.chromosome, 'GRCh37')
        coverage_corrector = CoverageCorrector(bias_detector.get_gc_content_coverage_map())

        observed_copy_number = pattern_occurrences / coverage_corrector.get_sequencing_mean_coverage()
        scaled_copy_number = coverage_corrector.get_scaled_coverage(self.reference_vntr, observed_copy_number)
        print('scaled copy number and observed copy number: ', scaled_copy_number, observed_copy_number)
        print('unmapped reads influence: ', scaled_copy_number * vntr_bp_in_unmapped_reads.value /
              (vntr_bp_in_mapped_reads + vntr_bp_in_unmapped_reads.value))
        return scaled_copy_number

    def find_repeat_count_from_short_reads(self, short_read_files, working_directory='./'):
        """
        Map short read sequencing data to human reference genome (hg19) and call find_repeat_count_from_alignment_file
        :param short_read_files: short read sequencing data
        :param working_directory: directory for generating the outputs
        """
        alignment_file = '' + short_read_files
        # TODO: use bowtie2 to map short reads to hg19
        return self.find_repeat_count_from_alignment_file(alignment_file, working_directory)

    def find_accuracy(self, samfile='original_reads/paired_dat.sam'):
        """Find sensitivity and false positive reads for a set of simulated data
        """
        reference_end_pos = self.reference_vntr.start_point + self.reference_vntr.get_length()
        related_reads, read_count = get_related_reads_and_read_count_in_samfile(self.reference_vntr.pattern,
                                                                                self.reference_vntr.start_point,
                                                                                read_file=samfile,
                                                                                pattern_end=reference_end_pos)
        # TODO
        selected_reads = []
        occurrences = 0
        avg_coverage = 1
        true_positives = [read for read in selected_reads if read in related_reads]
        false_positives = [read for read in selected_reads if read not in true_positives]
        false_negatives = [read for read in related_reads if read not in selected_reads]
        # print('TP:', len(true_positives), 'FP:', len(false_positives), 'selected:', len(selected_reads))
        # print('FN:', len(false_negatives))
        sensitivity = float(len(true_positives)) / len(related_reads) if len(related_reads) > 0 else 0
        if sensitivity > 0.9:
            print(sensitivity, len(false_positives))
        if 1 > sensitivity > 0.9 and len(false_negatives) > 0 and len(false_positives) > 0:
            print('sensitivity ', sensitivity, ' FN:', false_negatives[0], ' FP:', false_positives[0])
        with open('FP_and_sensitivity_HMM_read_scoring_method.txt', 'a') as outfile:
            outfile.write('%s\t%s\t%s\t%s\t%s\n' % (
                len(false_positives), sensitivity, self.reference_vntr.id, len(self.reference_vntr.pattern),
                len(true_positives)))
        error = abs(len(self.reference_vntr.get_repeat_segments()) - occurrences / avg_coverage)
        print(error)
