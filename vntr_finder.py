import logging
import numpy
import os
from multiprocessing import Process, Manager, Value, Semaphore
from random import random
from uuid import uuid4

import pysam
from Bio import pairwise2
from Bio.Seq import Seq

from blast_wrapper import get_blast_matched_ids, make_blast_database
from coverage_bias import CoverageBiasDetector, CoverageCorrector
from hmm_utils import *
from pacbio_haplotyper import PacBioHaplotyper
from pomegranate import HiddenMarkovModel as Model
from profiler import time_usage
from sam_utils import get_reference_genome_of_alignment_file
from sam_utils import get_related_reads_and_read_count_in_samfile
import settings
from utils import is_low_quality_read


class SelectedRead:
    def __init__(self, sequence, logp, vpath, mapq=None, reference_start=None):
        self.sequence = sequence
        self.logp = logp
        self.vpath = vpath
        self.mapq = mapq
        self.is_mapped = reference_start is not None

    def is_mapped(self):
        return self.is_mapped


class VNTRFinder:
    """Find the VNTR structure of a reference VNTR in NGS data of the donor."""

    def __init__(self, reference_vntr):
        self.reference_vntr = reference_vntr
        self.min_repeat_bp_to_add_read = 2
        if len(self.reference_vntr.pattern) < 30:
            self.min_repeat_bp_to_add_read = 2
        self.min_repeat_bp_to_count_repeats = 2

        self.minimum_left_flanking_size = {}
        self.minimum_right_flanking_size = {119: 19, 1214: 12}

        self.vntr_start = self.reference_vntr.start_point
        self.vntr_end = self.vntr_start + self.reference_vntr.get_length()

    @time_usage
    def build_vntr_matcher_hmm(self, copies, flanking_region_size=100):
        patterns = self.reference_vntr.get_repeat_segments()
        left_flanking_region = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking_region = self.reference_vntr.right_flanking_region[:flanking_region_size]

        vntr_matcher = get_read_matcher_model(left_flanking_region, right_flanking_region, patterns, copies)
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

        flanking_region_size = read_length - 10
        vntr_matcher = self.build_vntr_matcher_hmm(copies, flanking_region_size)

        json_str = vntr_matcher.to_json()
        with open(stored_hmm_file, 'w') as outfile:
            outfile.write(json_str)
        return vntr_matcher

    @time_usage
    def filter_reads_with_keyword_matching(self, working_directory, read_file, short_reads=True):
        db_name = 'blast_db__' + os.path.basename(read_file)
        blast_db_name = working_directory + db_name
        empty_db = False
        if not os.path.exists(blast_db_name + '.nsq') and not os.path.exists(blast_db_name + '.nal'):
            empty_db = make_blast_database(read_file, blast_db_name)

        word_size = int(len(self.reference_vntr.pattern)/3)
        if word_size > 11:
            word_size = 11
        if word_size < 5:
            word_size = 5
        word_size = str(word_size)

        search_results = []
        blast_ids = set([])
        search_id = str(uuid4()) + str(self.reference_vntr.id)
        queries = self.reference_vntr.get_repeat_segments()
        if len(self.reference_vntr.pattern) < 10:
            min_copies = int(10 / len(self.reference_vntr.pattern))
            queries = [self.reference_vntr.pattern * min_copies]
        identity_cutoff = '40'
        if not short_reads:
            queries = [self.reference_vntr.left_flanking_region[-80:], self.reference_vntr.right_flanking_region[:80]]
            word_size = str('10')
            identity_cutoff = '70'
        if not empty_db:
            for query in queries:
                search_result = get_blast_matched_ids(query, blast_db_name, max_seq='50000', word_size=word_size,
                                                      evalue=10, search_id=search_id, identity_cutoff=identity_cutoff)
                search_results.append(search_result)

            if short_reads:
                for search_result in search_results:
                    blast_ids |= search_result
            else:
                blast_ids = search_results[0] & search_results[1]

        logging.info('blast selected %s reads for %s' % (len(blast_ids), self.reference_vntr.id))
        if len(blast_ids) == len(self.reference_vntr.get_repeat_segments()) * 50 * 1000:
            logging.error('maximum number of read selected in filtering for pattern %s' % self.reference_vntr.id)
        return blast_ids

    @staticmethod
    def add_hmm_score_to_list(sema, hmm, read, result_scores):
        logp, vpath = hmm.viterbi(str(read.seq))
        rev_logp, rev_vpath = hmm.viterbi(str(Seq(str(read.seq)).reverse_complement()))
        if logp < rev_logp:
            logp = rev_logp
        result_scores.append(logp)
        sema.release()

    def is_true_read(self, read):
        read_start = read.reference_start
        read_end = read.reference_end if read.reference_end else read_start + len(read.seq)
        reference_name = read.reference_name
        if not reference_name.startswith('chr'):
            reference_name = 'chr' + reference_name
        if reference_name == self.reference_vntr.chromosome and self.vntr_start - len(read.seq) < read_start < self.vntr_end:
            return True
        return False

    def find_score_distribution_of_ref(self, samfile, reference, hmm, false_scores, true_scores):
        process_list = []
        sema = Semaphore(settings.CORES)
        for read in samfile.fetch(reference, multiple_iterators=True):
            if read.is_unmapped:
                continue
            if read.seq.count('N') > 0:
                continue

            if self.is_true_read(read):
                sema.acquire()
                p = Process(target=VNTRFinder.add_hmm_score_to_list, args=(sema, hmm, read, true_scores))
            else:
                if random() > settings.SCORE_FINDING_READS_FRACTION:
                    continue
                sema.acquire()
                p = Process(target=VNTRFinder.add_hmm_score_to_list, args=(sema, hmm, read, false_scores))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

    def save_scores(self, true_scores, false_scores, alignment_file):
        with open('true_scores_dist_%s_%s' % (self.reference_vntr.id, os.path.basename(alignment_file)), 'w') as out:
            for score in true_scores:
                out.write('%.4f\n' % score)
        with open('false_scores_dist_%s_%s' % (self.reference_vntr.id, os.path.basename(alignment_file)), 'w') as out:
            for score in false_scores:
                out.write('%.4f\n' % score)

    @time_usage
    def calculate_min_score_to_select_a_read(self, hmm, alignment_file):
        """Calculate the score distribution of false positive reads
        and return score to select the 1e-8 percentile of the distribution
        """
        process_list = []
        manager = Manager()
        false_scores = manager.list()
        true_scores = manager.list()
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        refs = [ref for ref in samfile.references if ref in settings.CHROMOSOMES or 'chr' + ref in settings.CHROMOSOMES]
        for ref in refs:
            p = Process(target=self.find_score_distribution_of_ref, args=(samfile, ref, hmm, false_scores, true_scores))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        if settings.SAVE_SCORE_DISTRIBUTION:
            self.save_scores(true_scores, false_scores, alignment_file)

        score = numpy.percentile(false_scores, 100 - settings.SCORE_SELECTION_PERCENTILE)
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

        logging.debug('Minimum score is not precomputed for vntr id: %s' % self.reference_vntr.id)
        score = self.calculate_min_score_to_select_a_read(hmm, alignment_file)
        logging.debug('computed score: %s' % score)
        with open(stored_scores_file, 'a') as outfile:
            outfile.write('%s %s\n' % (settings.SCORE_FINDING_READS_FRACTION, score))

        return score

    def process_unmapped_read(self, sema, read_segment, hmm, min_score_to_count_read,
                              vntr_bp_in_unmapped_reads, selected_reads, best_seq):
        if read_segment.seq.count('N') <= 0:
            sequence = str(read_segment.seq)
            logp, vpath = hmm.viterbi(sequence)
            rev_logp, rev_vpath = hmm.viterbi(str(read_segment.seq.reverse_complement()))
            if logp < rev_logp:
                sequence = str(read_segment.seq.reverse_complement())
                logp = rev_logp
                vpath = rev_vpath
            if logp > best_seq['logp']:
                best_seq['logp'] = logp
                best_seq['seq'] = sequence
                best_seq['vpath'] = vpath
            repeat_bps = get_number_of_repeat_bp_matches_in_vpath(vpath)
            if logp > min_score_to_count_read:
                if repeat_bps > self.min_repeat_bp_to_count_repeats:
                    vntr_bp_in_unmapped_reads.value += repeat_bps
                if repeat_bps > self.min_repeat_bp_to_add_read:
                    selected_reads.append(SelectedRead(sequence, logp, vpath))
        sema.release()

    def find_frameshift_from_selected_reads(self, selected_reads):
        mutations = {}
        repeating_bps_in_data = 0
        repeats_lengths_distribution = []
        for read in selected_reads:
            visited_states = [state.name for idx, state in read.vpath[1:-1]]
            repeats_lengths = get_repeating_pattern_lengths(visited_states)
            repeats_lengths_distribution += repeats_lengths
            current_repeat = None
            repeating_bps_in_data += get_number_of_repeat_bp_matches_in_vpath(read.vpath)
            for i in range(len(visited_states)):
                if visited_states[i].endswith('fix') or visited_states[i].startswith('M'):
                    continue
                if visited_states[i].startswith('unit_start'):
                    if current_repeat is None:
                        current_repeat = 0
                    else:
                        current_repeat += 1
                if current_repeat is None or current_repeat >= len(repeats_lengths):
                    continue
                if not visited_states[i].startswith('I') and not visited_states[i].startswith('D'):
                    continue
                if repeats_lengths[current_repeat] == len(self.reference_vntr.pattern):
                    continue
                state = visited_states[i].split('_')[0]
                if state.startswith('I'):
                    state += get_emitted_basepair_from_visited_states(visited_states[i], visited_states, read.sequence)
                if abs(repeats_lengths[current_repeat] - len(self.reference_vntr.pattern)) <= 2:
                    if state not in mutations.keys():
                        mutations[state] = 0
                    mutations[state] += 1
        sorted_mutations = sorted(mutations.items(), key=lambda x: x[1])
        logging.debug('sorted mutations: %s ' % sorted_mutations)
        frameshift_candidate = sorted_mutations[-1] if len(sorted_mutations) else (None, 0)
        logging.info(sorted(repeats_lengths_distribution))
        logging.info('Frameshift Candidate and Occurrence %s: %s' % frameshift_candidate)
        logging.info('Observed repeating base pairs in data: %s' % repeating_bps_in_data)
        avg_bp_coverage = float(repeating_bps_in_data) / self.reference_vntr.get_length()
        logging.info('Average coverage for each base pair: %s' % avg_bp_coverage)
        if frameshift_candidate[1] > avg_bp_coverage / 4:
            logging.info('There is a frameshift at %s' % frameshift_candidate[0])
            return frameshift_candidate[0]
        return None

    def read_flanks_repeats_with_confidence(self, vpath):
        minimum_left_flanking = 5
        minimum_right_flanking = 5
        if self.reference_vntr.id in self.minimum_left_flanking_size:
            minimum_left_flanking = self.minimum_left_flanking_size[self.reference_vntr.id]
        if self.reference_vntr.id in self.minimum_right_flanking_size:
            minimum_right_flanking = self.minimum_right_flanking_size[self.reference_vntr.id]

        if get_left_flanking_region_size_in_vpath(vpath) > minimum_left_flanking:
            if get_right_flanking_region_size_in_vpath(vpath) > minimum_right_flanking:
                return True
        return False

    def check_if_flanking_regions_align_to_str(self, read_str, length_distribution, spanning_reads):
        flanking_region_size = 100
        left_flanking = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking = self.reference_vntr.right_flanking_region[:flanking_region_size]
        left_alignments = pairwise2.align.localms(read_str, left_flanking, 1, -1, -1, -1)
        if len(left_alignments) < 1:
            return
        left_align = left_alignments[0]
        if left_align[2] < len(left_flanking) * (1 - settings.MAX_ERROR_RATE):
            return
        right_alignments = pairwise2.align.localms(read_str, right_flanking, 1, -1, -1, -1)
        if len(right_alignments) < 1:
            return
        right_align = right_alignments[0]
        if right_align[2] < len(right_flanking) * (1 - settings.MAX_ERROR_RATE):
            return
        if right_align[3] < left_align[3]:
            return
        spanning_reads.append(read_str[left_align[3]:right_align[3]+flanking_region_size])
        length_distribution.append(right_align[3] - (left_align[3] + flanking_region_size))

    def check_if_pacbio_read_spans_vntr(self, sema, read, length_distribution, spanning_reads):
        self.check_if_flanking_regions_align_to_str(str(read.seq), length_distribution, spanning_reads)
        reverse_complement_str = str(Seq(str(read.seq)).reverse_complement())
        self.check_if_flanking_regions_align_to_str(reverse_complement_str, length_distribution, spanning_reads)
        sema.release()

    @time_usage
    def get_spanning_reads_of_unaligned_pacbio_reads(self, unmapped_filtered_reads):
        sema = Semaphore(settings.CORES)
        manager = Manager()
        shared_length_distribution = manager.list()
        shared_spanning_reads = manager.list()

        process_list = []
        for read in unmapped_filtered_reads:
            sema.acquire()
            p = Process(target=self.check_if_pacbio_read_spans_vntr, args=(sema, read, shared_length_distribution,
                                                                           shared_spanning_reads))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()
        logging.info('length_distribution of unmapped spanning reads: %s' % list(shared_length_distribution))
        return list(shared_spanning_reads), list(shared_length_distribution)

    @time_usage
    def get_spanning_reads_of_aligned_pacbio_reads(self, alignment_file):
        sema = Semaphore(settings.CORES)
        manager = Manager()
        shared_length_distribution = manager.list()
        mapped_spanning_reads = manager.list()

        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        region_start = vntr_start
        region_end = vntr_end
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        reference = get_reference_genome_of_alignment_file(samfile)
        chromosome = self.reference_vntr.chromosome if reference == 'HG19' else self.reference_vntr.chromosome[3:]
        process_list = []
        for read in samfile.fetch(chromosome, region_start, region_end):
            sema.acquire()
            p = Process(target=self.check_if_pacbio_read_spans_vntr, args=(sema, read, shared_length_distribution,
                                                                           mapped_spanning_reads))
            process_list.append(p)
            p.start()

        for p in process_list:
            p.join()

        logging.info('length_distribution of mapped spanning reads: %s' % list(shared_length_distribution))
        return list(mapped_spanning_reads)

    @time_usage
    def get_haplotype_copy_numbers_from_spanning_reads(self, spanning_reads):
        if len(spanning_reads) < 1:
            logging.info('There is no spanning read')
            return None
        max_length = 0
        import numpy
        median = numpy.median([len(l) for l in spanning_reads])
        spanning_reads = [r for r in spanning_reads if len(r) < median * 2]
        for read in spanning_reads:
            if len(read) - 100 > max_length:
                max_length = len(read) - 100
        max_copies = int(round(max_length / float(len(self.reference_vntr.pattern))))
        # max_copies = min(max_copies, 2 * len(self.reference_vntr.get_repeat_segments()))
        vntr_matcher = self.build_vntr_matcher_hmm(max_copies)
        haplotyper = PacBioHaplotyper(spanning_reads)
        haplotypes = haplotyper.get_error_corrected_haplotypes()
        copy_numbers = []
        for haplotype in haplotypes:
            # print('haplotype: %s' % haplotype)
            logp, vpath = vntr_matcher.viterbi(haplotype)
            rev_logp, rev_vpath = vntr_matcher.viterbi(str(Seq(haplotype).reverse_complement()))
            if logp < rev_logp:
                vpath = rev_vpath
            copy_numbers.append(get_number_of_repeats_in_vpath(vpath))
        return copy_numbers

    @time_usage
    def find_repeat_count_from_pacbio_alignment_file(self, alignment_file, unmapped_filtered_reads):
        logging.debug('finding repeat count from pacbio alignment file for %s' % self.reference_vntr.id)

        unaligned_spanning_reads, length_dist = self.get_spanning_reads_of_unaligned_pacbio_reads(unmapped_filtered_reads)
        mapped_spanning_reads = self.get_spanning_reads_of_aligned_pacbio_reads(alignment_file)

        spanning_reads = mapped_spanning_reads + unaligned_spanning_reads
        copy_numbers = self.get_haplotype_copy_numbers_from_spanning_reads(spanning_reads)
        return copy_numbers

    @time_usage
    def find_repeat_count_from_pacbio_reads(self, unmapped_filtered_reads, naive=False):
        logging.debug('finding repeat count from pacbio reads file for %s' % self.reference_vntr.id)
        spanning_reads, length_dist = self.get_spanning_reads_of_unaligned_pacbio_reads(unmapped_filtered_reads)
        if naive:
            if len(length_dist):
                copy_numbers = [round(sum(length_dist) / float(len(length_dist)) / len(self.reference_vntr.pattern))] * 2
            else:
                copy_numbers = None
        else:
            copy_numbers = self.get_haplotype_copy_numbers_from_spanning_reads(spanning_reads)
        return copy_numbers

    @time_usage
    def select_illumina_reads(self, alignment_file, unmapped_filtered_reads):
        hmm = None
        min_score_to_count_read = None
        sema = Semaphore(settings.CORES)
        manager = Manager()
        selected_reads = manager.list()
        vntr_bp_in_unmapped_reads = Value('d', 0.0)

        number_of_reads = 0
        read_length = 150

        process_list = []

        best_seq = manager.dict()
        best_seq['logp'] = -10e8
        best_seq['vpath'] = ''
        best_seq['seq'] = ''

        for read_segment in unmapped_filtered_reads:
            if number_of_reads == 0:
                read_length = len(str(read_segment.seq))
            number_of_reads += 1
            if not hmm:
                hmm = self.get_vntr_matcher_hmm(read_length=read_length)
                min_score_to_count_read = self.get_min_score_to_select_a_read(hmm, alignment_file, read_length)

            if len(read_segment.seq) < read_length:
                continue

            sema.acquire()
            p = Process(target=self.process_unmapped_read, args=(sema, read_segment, hmm, min_score_to_count_read,
                                                                 vntr_bp_in_unmapped_reads, selected_reads, best_seq))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        logging.debug('vntr base pairs in unmapped reads: %s' % vntr_bp_in_unmapped_reads.value)
        logging.debug('highest logp in unmapped reads: %s' % best_seq['logp'])
        logging.debug('best sequence %s' % best_seq['seq'])
        logging.debug('best vpath: %s' % [state.name for idx, state in list(best_seq['vpath'])[1:-1]])

        vntr_bp_in_mapped_reads = 0
        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(alignment_file, read_mode)
        reference = get_reference_genome_of_alignment_file(samfile)
        chromosome = self.reference_vntr.chromosome if reference == 'HG19' else self.reference_vntr.chromosome[3:]
        for read in samfile.fetch(chromosome, vntr_start, vntr_end):
            if not hmm:
                read_length = len(read.seq)
                hmm = self.get_vntr_matcher_hmm(read_length=read_length)
                min_score_to_count_read = self.get_min_score_to_select_a_read(hmm, alignment_file, read_length)

            if read.is_unmapped:
                continue
            if len(read.seq) < int(read_length * 0.9):
                logging.debug('Rejecting read for short length: %s' % read.seq)
                continue
            read_end = read.reference_end if read.reference_end else read.reference_start + len(read.seq)
            if vntr_start - read_length < read.reference_start < vntr_end or vntr_start < read_end < vntr_end:
                if read.seq.count('N') <= 0:
                    sequence = str(read.seq)
                    logp, vpath = hmm.viterbi(sequence)
                    rev_logp, rev_vpath = hmm.viterbi(str(Seq(read.seq).reverse_complement()))
                    if logp < rev_logp:
                        sequence = str(Seq(read.seq).reverse_complement())
                        logp = rev_logp
                        vpath = rev_vpath
                    if is_low_quality_read(read) and logp < min_score_to_count_read:
                        logging.debug('Rejected Read: %s' % sequence)
                        continue
                    selected_reads.append(SelectedRead(sequence, logp, vpath, read.mapq, read.reference_start))
                end = min(read_end, vntr_end)
                start = max(read.reference_start, vntr_start)
                vntr_bp_in_mapped_reads += end - start
        logging.debug('vntr base pairs in mapped reads: %s' % vntr_bp_in_mapped_reads)

        return selected_reads

    @time_usage
    def find_frameshift_from_alignment_file(self, alignment_file, unmapped_filtered_reads):
        logging.debug('finding frameshift from alignment file for %s' % self.reference_vntr.id)

        selected_reads = self.select_illumina_reads(alignment_file, unmapped_filtered_reads)
        return self.find_frameshift_from_selected_reads(selected_reads)

    @time_usage
    def find_repeat_count_from_alignment_file(self, alignment_file, unmapped_filtered_reads):
        logging.debug('finding repeat count from alignment file for %s' % self.reference_vntr.id)

        selected_reads = self.select_illumina_reads(alignment_file, unmapped_filtered_reads)

        flanked_repeats = []
        observed_repeats = []
        for selected_read in selected_reads:
            repeats = get_number_of_repeats_in_vpath(selected_read.vpath)
            logging.debug('logp of read: %s' % str(selected_read.logp))
            logging.debug('left flankign size: %s' % get_left_flanking_region_size_in_vpath(selected_read.vpath))
            logging.debug('right flanking size: %s' % get_right_flanking_region_size_in_vpath(selected_read.vpath))
            logging.debug('repeating bp: %s' % get_number_of_repeat_bp_matches_in_vpath(selected_read.vpath))
            logging.debug(selected_read.sequence)
            visited_states = [state.name for idx, state in selected_read.vpath[1:-1]]
            # logging.debug('%s' % visited_states)
            if self.read_flanks_repeats_with_confidence(selected_read.vpath):
                logging.debug('spanning read: %s ' % selected_read.sequence)
                logging.debug('visited states :%s' % [state.name for idx, state in selected_read.vpath[1:-1]])
                logging.debug('repeats: %s' % repeats)
                flanked_repeats.append(repeats)
            observed_repeats.append(repeats)
        print('flanked repeats:', flanked_repeats)
        print('observed repeats:', sorted(observed_repeats))
        observed_repeats = reversed(sorted(observed_repeats))
        result = set([])
        for repeat in flanked_repeats:
            result.add(repeat)
        result = list(result)
        for repeat in observed_repeats:
            if len(result) >= 2:
                break
            if len(result) > 0 and repeat < result[0]:
                result.add(result[0])
                break
            result.add(repeat)

        return list(result)
        # TODO: separate methods

        total_counted_vntr_bp = vntr_bp_in_unmapped_reads.value + vntr_bp_in_mapped_reads
        pattern_occurrences = total_counted_vntr_bp / float(len(self.reference_vntr.pattern))
        bias_detector = CoverageBiasDetector(alignment_file, self.reference_vntr.chromosome, reference)
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
