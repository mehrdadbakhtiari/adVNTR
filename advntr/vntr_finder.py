from collections import Counter
import logging
import numpy
import os
from multiprocessing import Process, Manager, Value, Semaphore
from random import random

import pysam
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO

from advntr.coverage_bias import CoverageBiasDetector, CoverageCorrector
from advntr.hmm_utils import *
from advntr.pacbio_haplotyper import PacBioHaplotyper
from advntr.profiler import time_usage
from advntr.sam_utils import get_reference_genome_of_alignment_file, get_related_reads_and_read_count_in_samfile
from advntr import settings
from advntr.utils import is_low_quality_read
from pomegranate import HiddenMarkovModel as Model


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

    def __init__(self, reference_vntr, is_haploid=False, reference_filename=None):
        self.reference_vntr = reference_vntr
        self.is_haploid = is_haploid
        self.reference_filename = reference_filename
        self.min_repeat_bp_to_add_read = 2
        if len(self.reference_vntr.pattern) < 30:
            self.min_repeat_bp_to_add_read = 2
        self.min_repeat_bp_to_count_repeats = 2

        self.minimum_left_flanking_size = {}
        self.minimum_right_flanking_size = {69212: 19, 532789: 12, 400825: 10, 468671: 10}

        self.vntr_start = self.reference_vntr.start_point
        self.vntr_end = self.vntr_start + self.reference_vntr.get_length()

    def get_copies_for_hmm(self, read_length):
        return int(round(float(read_length) / len(self.reference_vntr.pattern) + 0.5))

    @staticmethod
    def get_alignment_file_read_mode(alignment_file):
        read_mode = 'r' if alignment_file.endswith('sam') else 'rb'
        if alignment_file.endswith('cram'):
            read_mode = 'rc'
        return read_mode

    @time_usage
    def build_vntr_matcher_hmm(self, copies, flanking_region_size=100):
        patterns = self.reference_vntr.get_repeat_segments()
        left_flanking_region = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking_region = self.reference_vntr.right_flanking_region[:flanking_region_size]

        vntr_matcher = get_read_matcher_model(left_flanking_region, right_flanking_region, patterns, copies)
        return vntr_matcher

    def get_vntr_matcher_hmm(self, read_length):
        """Try to load trained HMM for this VNTR
        If there was no trained HMM, it will build one and store it for later usage
        """
        logging.info('Using read length %s' % read_length)
        copies = self.get_copies_for_hmm(read_length)

        base_name = str(self.reference_vntr.id) + '_' + str(read_length) + '.json'
        stored_hmm_file = settings.TRAINED_HMMS_DIR + base_name
        if settings.USE_TRAINED_HMMS and os.path.isfile(stored_hmm_file):
            model = Model()
            model = model.from_json(stored_hmm_file)
            return model

        flanking_region_size = read_length
        vntr_matcher = self.build_vntr_matcher_hmm(copies, flanking_region_size)

        if settings.USE_TRAINED_HMMS:
            json_str = vntr_matcher.to_json()
            with open(stored_hmm_file, 'w') as outfile:
                outfile.write(json_str)
        return vntr_matcher

    def get_keywords_for_filtering(self, short_reads=True, keyword_size=21):
        vntr = ''.join(self.reference_vntr.get_repeat_segments())
        if len(vntr) < keyword_size:
            min_copies = int(keyword_size / len(vntr)) + 1
            vntr = str(vntr) * min_copies
        locus = self.reference_vntr.left_flanking_region[-15:] + vntr + self.reference_vntr.right_flanking_region[:15]
        queries = []
        step_size = 5 if len(self.reference_vntr.pattern) != 5 else 6
        for i in range(0, len(locus) - keyword_size + 1, step_size):
            queries.append(locus[i:i+keyword_size])

        if not short_reads:
            queries = [self.reference_vntr.left_flanking_region[-80:], self.reference_vntr.right_flanking_region[:80]]
        queries = set(queries)
        return queries

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
        reference_name = read.reference_name
        if not reference_name.startswith('chr'):
            reference_name = 'chr' + reference_name
        if reference_name == self.reference_vntr.chromosome and self.vntr_start - len(read.seq) < read_start < self.vntr_end:
            return True
        return False

    def get_min_score_to_select_a_read(self, read_length):
        if self.reference_vntr.scaled_score is None or self.reference_vntr.scaled_score == 0:
            return None
        return self.reference_vntr.scaled_score * read_length

    @staticmethod
    def recruit_read(logp, vpath, min_score_to_count_read, read_length):
        if min_score_to_count_read is not None and logp > min_score_to_count_read:
            return True
        matches = get_number_of_matches_in_vpath(vpath)
        if min_score_to_count_read is None and matches >= 0.9 * read_length and logp > -read_length:
            return True
        return False

    def process_unmapped_read(self, sema, read_segment, hmm, recruitment_score, vntr_bp_in_unmapped_reads,
                              selected_reads, compute_reverse=True):
        if read_segment.count('N') <= 0:
            sequence = read_segment.upper()
            logp, vpath = hmm.viterbi(sequence)
            if compute_reverse:
                reverse_sequence = str(Seq(sequence).reverse_complement())
                rev_logp, rev_vpath = hmm.viterbi(reverse_sequence)
                if logp < rev_logp:
                    sequence = reverse_sequence
                    logp = rev_logp
                    vpath = rev_vpath
            repeat_bps = get_number_of_repeat_bp_matches_in_vpath(vpath)
            if self.recruit_read(logp, vpath, recruitment_score, len(sequence)):
                if repeat_bps > self.min_repeat_bp_to_count_repeats:
                    vntr_bp_in_unmapped_reads.value += repeat_bps
                if repeat_bps > self.min_repeat_bp_to_add_read:
                    selected_reads.append(SelectedRead(sequence, logp, vpath))
        if sema is not None:
            sema.release()

    def identify_frameshift(self, location_coverage, observed_indel_transitions, expected_indels, error_rate=0.01):
        if observed_indel_transitions >= location_coverage:
            return True
        from scipy.stats import binom
        sequencing_error_prob = binom.pmf(observed_indel_transitions, location_coverage, error_rate)
        frameshift_prob = binom.pmf(observed_indel_transitions, location_coverage, expected_indels)
        prob = sequencing_error_prob / frameshift_prob
        return prob < 0.01

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
        avg_bp_coverage = float(repeating_bps_in_data) / self.reference_vntr.get_length() / 2
        logging.info('Average coverage for each base pair: %s' % avg_bp_coverage)

        expected_indel_transitions = 1 / avg_bp_coverage
        if self.identify_frameshift(avg_bp_coverage, frameshift_candidate[1], expected_indel_transitions):
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
        min_left, max_left = 10e9, 0
        for aln in left_alignments:
            if aln[2] < len(left_flanking) * (1 - settings.MAX_ERROR_RATE):
                continue
            min_left = min(min_left, aln[3])
            max_left = max(max_left, aln[3])
        if max_left - min_left > 30:
            with open('vntr_complex.txt', 'a') as out:
                out.write('%s %s\n' % (self.reference_vntr.id, max_left - min_left))
        left_align = left_alignments[0]
        if left_align[2] < len(left_flanking) * (1 - settings.MAX_ERROR_RATE):
            return

        right_alignments = pairwise2.align.localms(read_str, right_flanking, 1, -1, -1, -1)
        if len(right_alignments) < 1:
            return
        min_right, max_right = 10e9, 0
        for aln in right_alignments:
            if aln[2] < len(right_flanking) * (1 - settings.MAX_ERROR_RATE):
                continue
            min_right = min(min_right, aln[3])
            max_right = max(max_right, aln[3])
        if max_right - min_right > 30:
            with open('vntr_complex.txt', 'a') as out:
                out.write('%s %s\n' % (self.reference_vntr.id, max_right - min_right))
        right_align = right_alignments[0]
        if right_align[2] < len(right_flanking) * (1 - settings.MAX_ERROR_RATE):
            return

        if right_align[3] < left_align[3]:
            return
        spanning_reads.append(read_str[left_align[3]:right_align[3]+flanking_region_size])
        length_distribution.append(right_align[3] - (left_align[3] + flanking_region_size))

    def check_if_pacbio_read_spans_vntr(self, sema, read, length_distribution, spanning_reads):
        self.check_if_flanking_regions_align_to_str(str(read.seq).upper(), length_distribution, spanning_reads)
        reverse_complement_str = str(Seq(str(read.seq)).reverse_complement())
        self.check_if_flanking_regions_align_to_str(reverse_complement_str.upper(), length_distribution, spanning_reads)
        sema.release()

    def check_if_pacbio_mapped_read_spans_vntr(self, sema, read, length_distribution, spanning_reads):
        flanking_region_size = 100
        region_start = self.reference_vntr.start_point - flanking_region_size
        region_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        if read.get_reference_positions()[0] < region_start and read.get_reference_positions()[-1] > region_end:
            read_region_start = None
            read_region_end = None
            for read_pos, ref_pos in enumerate(read.get_reference_positions()):
                if ref_pos >= region_start and read_region_start is None:
                    read_region_start = read_pos
                if ref_pos >= region_end and read_region_end is None:
                    read_region_end = read_pos
            if read_region_start is not None and read_region_end is not None:
                result = read.seq[read_region_start:read_region_end+flanking_region_size]
                if read.is_reverse:
                    result = str(Seq(result).reverse_complement())
                spanning_reads.append(result)
                length_distribution.append(len(result) - flanking_region_size * 2)
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
        length_distribution = manager.list()
        mapped_spanning_reads = manager.list()

        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        region_start = vntr_start
        region_end = vntr_end
        read_mode = self.get_alignment_file_read_mode(alignment_file)
        samfile = pysam.AlignmentFile(alignment_file, read_mode, reference_filename=self.reference_filename)
        reference = get_reference_genome_of_alignment_file(samfile)
        chromosome = self.reference_vntr.chromosome if reference == 'HG19' else self.reference_vntr.chromosome[3:]
        process_list = []
        for read in samfile.fetch(chromosome, region_start, region_end):
            sema.acquire()
            p = Process(target=self.check_if_pacbio_mapped_read_spans_vntr, args=(sema, read, length_distribution,
                                                                                  mapped_spanning_reads))
            process_list.append(p)
            p.start()

        for p in process_list:
            p.join()

        logging.info('length_distribution of mapped spanning reads: %s' % list(length_distribution))
        return list(mapped_spanning_reads)

    def get_conditional_likelihood(self, ck, ci, cj, ru_counts, r, r_e):
        if ck == ci == cj:
            return 1-r
        if cj == 0:  # CHECK LATER
            return 0.5 * (1-r)
        if ck == ci:
            return 0.5 * ((1-r) + r_e ** abs(ck-cj))
        if ck == cj:
            return 0.5 * ((1-r) + r_e ** abs(ck-ci))
        if ck != ci and ck != cj:
            return 0.5 * (r_e ** abs(ck-ci) + r_e ** abs(ck-cj))

    def find_genotype_based_on_observed_repeats(self, observed_copy_numbers):
        ru_counts = {}
        for cn in observed_copy_numbers:
            if cn not in ru_counts.keys():
                ru_counts[cn] = 0
            ru_counts[cn] += 1
        if len(ru_counts.keys()) < 2:
            priors = 0.5
            ru_counts[0] = 1
        else:
            priors = 1.0 / (len(ru_counts.keys()) * (len(ru_counts.keys())-1) / 2)
        import operator
        ru_counts = sorted(ru_counts.items(), key=operator.itemgetter(1), reverse=True)
        r = 0.03
        r_e = r / (2 + r)
        prs = {}
        for ck, occ in ru_counts:
            if ck == 0:
                continue
            for i in range(len(ru_counts)):
                ci = ru_counts[i][0]
                for j in range(len(ru_counts)):
                    if j < i:
                        continue
                    if self.is_haploid and i != j:
                        continue
                    cj = ru_counts[j][0]
                    if (ci, cj) not in prs.keys():
                        prs[(ci, cj)] = []
                    prs[(ci, cj)].append(self.get_conditional_likelihood(ck, ci, cj, ru_counts, r, r_e) ** occ)

        posteriors = {}
        import numpy
        for key in prs.keys():
            prs[key] = numpy.prod(numpy.array(prs[key]))
            posteriors[key] = prs[key] * priors

        sum_of_probs = sum(posteriors.values())

        max_prob = 1e-20
        result = None
        for key, value in posteriors.items():
            if value / sum_of_probs > max_prob:
                max_prob = value / sum_of_probs
                result = key

        logging.info('Maximum probability for genotyping: %s' % max_prob)
        return result

    def get_dominant_copy_numbers_from_spanning_reads(self, spanning_reads):
        if len(spanning_reads) < 1:
            logging.info('There is no spanning read')
            return None
        max_length = 0
        for read in spanning_reads:
            if len(read) - 100 > max_length:
                max_length = len(read) - 100
        max_copies = int(round(max_length / float(len(self.reference_vntr.pattern))))
        # max_copies = min(max_copies, 2 * len(self.reference_vntr.get_repeat_segments()))
        vntr_matcher = self.build_vntr_matcher_hmm(max_copies)
        observed_copy_numbers = []
        for haplotype in spanning_reads:
            logp, vpath = vntr_matcher.viterbi(haplotype)
            rev_logp, rev_vpath = vntr_matcher.viterbi(str(Seq(haplotype).reverse_complement()))
            if logp < rev_logp:
                vpath = rev_vpath
            observed_copy_numbers.append(get_number_of_repeats_in_vpath(vpath))

        logging.info('flanked repeats: %s' % observed_copy_numbers)
        return self.find_genotype_based_on_observed_repeats(observed_copy_numbers)

    @time_usage
    def get_haplotype_copy_numbers_from_spanning_reads(self, spanning_reads):
        if len(spanning_reads) < 1:
            logging.info('There is no spanning read')
            return None
        max_length = 0
        for read in spanning_reads:
            if len(read) - 100 > max_length:
                max_length = len(read) - 100
        max_copies = int(round(max_length / float(len(self.reference_vntr.pattern))))
        max_copies = min(max_copies, 2 * len(self.reference_vntr.get_repeat_segments()))
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

    def find_ru_counts_with_naive_approach(self, length_dist, spanning_reads):
        haplotyper = PacBioHaplotyper(spanning_reads)
        haplotypes = haplotyper.get_error_corrected_haplotypes(1)
        flanking_region_lengths = []
        new_spanning_reads = []
        if len(haplotypes) == 0:
            return None
        self.check_if_flanking_regions_align_to_str(haplotypes[0].upper(), flanking_region_lengths, new_spanning_reads)
        reverse_complement_str = str(Seq(haplotypes[0]).reverse_complement())
        self.check_if_flanking_regions_align_to_str(reverse_complement_str.upper(), flanking_region_lengths, new_spanning_reads)
        if len(flanking_region_lengths) > 0:
            return [round(flanking_region_lengths[0] / len(self.reference_vntr.pattern))] * 2
        else:
            return None

    def find_ru_counts_from_average_flanking_region_distance(self, length_dist):
        if len(length_dist):
            ru_counts_list = [round(length / len(self.reference_vntr.pattern)) for length in length_dist]
            ru_count_frequencies = Counter(ru_counts_list)
            copy_numbers = [ru_count_frequencies[0][0]]
            if len(ru_count_frequencies.keys()) > 1 and ru_count_frequencies[1][1] > ru_count_frequencies[0][1] / 5:
                copy_numbers.append(ru_count_frequencies[1][0])
            else:
                copy_numbers = copy_numbers * 2
        else:
            copy_numbers = None
        return copy_numbers

    @time_usage
    def find_repeat_count_from_pacbio_alignment_file(self, alignment_file, unmapped_filtered_reads):
        logging.debug('finding repeat count from pacbio alignment file for %s' % self.reference_vntr.id)

        unaligned_spanning_reads, length_dist = self.get_spanning_reads_of_unaligned_pacbio_reads(unmapped_filtered_reads)
        mapped_spanning_reads = self.get_spanning_reads_of_aligned_pacbio_reads(alignment_file)

        spanning_reads = mapped_spanning_reads + unaligned_spanning_reads
        copy_numbers = self.get_dominant_copy_numbers_from_spanning_reads(spanning_reads)
        return copy_numbers

    @time_usage
    def find_repeat_count_from_pacbio_reads(self, unmapped_filtered_reads, naive=False):
        logging.debug('finding repeat count from pacbio reads file for %s' % self.reference_vntr.id)
        spanning_reads, length_dist = self.get_spanning_reads_of_unaligned_pacbio_reads(unmapped_filtered_reads)
        if naive:
            copy_numbers = self.find_ru_counts_with_naive_approach(length_dist, spanning_reads)
        else:
            copy_numbers = self.get_dominant_copy_numbers_from_spanning_reads(spanning_reads)
        return copy_numbers

    @time_usage
    def iteratively_update_model(self, alignment_file, unmapped_filtered_reads, selected_reads, hmm):
        updated_selected_reads = selected_reads
        fitness = sum([read.logp for read in selected_reads])
        read_length = len(selected_reads[0].sequence)

        reference_repeats = []
        for reference_repeat in self.reference_vntr.get_repeat_segments():
            sequence = str(reference_repeat).upper()
            logp, vpath = hmm.viterbi(sequence)
            reference_repeats.append(SelectedRead(sequence, logp, vpath))

        logging.info('initial fitness: %s' % fitness)

        flanking_region_size = read_length
        left_flanking_region = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking_region = self.reference_vntr.right_flanking_region[:flanking_region_size]
        copies = self.get_copies_for_hmm(read_length)
        max_steps = 1000
        min_improvement = 1
        for i in range(max_steps):
            old_fitness = fitness
            current_vpaths = [(read.sequence, read.vpath) for read in updated_selected_reads + reference_repeats]
            hmm = get_read_matcher_model(left_flanking_region, right_flanking_region, None, copies, current_vpaths)
            updated_selected_reads = self.select_illumina_reads(alignment_file, unmapped_filtered_reads, False, hmm)
            fitness = sum([read.logp for read in selected_reads])

            if fitness - old_fitness < min_improvement:
                break

        logging.info('final fitness: %s' % fitness)
        return updated_selected_reads

    @time_usage
    def select_illumina_reads(self, alignment_file, unmapped_filtered_reads, update=False, hmm=None):
        recruitment_score = None
        sema = Semaphore(settings.CORES)
        manager = Manager()
        selected_reads = manager.list()
        vntr_bp_in_unmapped_reads = Value('d', 0.0)

        number_of_reads = 0
        read_length = 150

        process_list = []

        for read_segment in unmapped_filtered_reads:
            if number_of_reads == 0:
                read_length = len(str(read_segment.seq))
            number_of_reads += 1
            if not hmm:
                hmm = self.get_vntr_matcher_hmm(read_length=read_length)
            if not recruitment_score:
                recruitment_score = self.get_min_score_to_select_a_read(read_length)

            if len(read_segment.seq) < read_length:
                continue

            sema.acquire()
            p = Process(target=self.process_unmapped_read, args=(sema, str(read_segment.seq), hmm, recruitment_score,
                                                                 vntr_bp_in_unmapped_reads, selected_reads))
            process_list.append(p)
            p.start()
        for p in process_list:
            p.join()

        logging.debug('vntr base pairs in unmapped reads: %s' % vntr_bp_in_unmapped_reads.value)

        vntr_bp_in_mapped_reads = 0
        vntr_start = self.reference_vntr.start_point
        vntr_end = self.reference_vntr.start_point + self.reference_vntr.get_length()
        read_mode = self.get_alignment_file_read_mode(alignment_file)
        samfile = pysam.AlignmentFile(alignment_file, read_mode, reference_filename=self.reference_filename)
        reference = get_reference_genome_of_alignment_file(samfile)
        chromosome = self.reference_vntr.chromosome if reference == 'HG19' else self.reference_vntr.chromosome[3:]
        for read in samfile.fetch(chromosome, vntr_start, vntr_end):
            if not recruitment_score:
                read_length = len(read.seq)
                recruitment_score = self.get_min_score_to_select_a_read(read_length)
            if not hmm:
                hmm = self.get_vntr_matcher_hmm(read_length=read_length)

            if read.is_unmapped:
                continue
            if len(read.seq) < int(read_length * 0.9):
                logging.debug('Rejecting read for short length: %s' % read.seq)
                continue
            read_end = read.reference_end if read.reference_end else read.reference_start + len(read.seq)
            if vntr_start - read_length < read.reference_start < vntr_end or vntr_start < read_end < vntr_end:
                if read.seq.count('N') <= 0:
                    sequence = str(read.seq).upper()
                    logp, vpath = hmm.viterbi(sequence)
                    rev_logp, rev_vpath = hmm.viterbi(str(Seq(read.seq).reverse_complement()).upper())
                    if logp < rev_logp:
                        sequence = str(Seq(read.seq).reverse_complement()).upper()
                        logp = rev_logp
                        vpath = rev_vpath
                    length = len(sequence)
                    if is_low_quality_read(read) and not self.recruit_read(logp, vpath, recruitment_score, length):
                        logging.debug('Rejected Read: %s' % sequence)
                        continue
                    selected_reads.append(SelectedRead(sequence, logp, vpath, read.mapq, read.reference_start))
                end = min(read_end, vntr_end)
                start = max(read.reference_start, vntr_start)
                vntr_bp_in_mapped_reads += end - start
        logging.debug('vntr base pairs in mapped reads: %s' % vntr_bp_in_mapped_reads)

        if update:
            selected_reads = self.iteratively_update_model(alignment_file, unmapped_filtered_reads, selected_reads, hmm)

        return selected_reads

    @time_usage
    def find_frameshift_from_alignment_file(self, alignment_file, unmapped_filtered_reads):
        logging.debug('finding frameshift from alignment file for %s' % self.reference_vntr.id)

        selected_reads = self.select_illumina_reads(alignment_file, unmapped_filtered_reads)
        return self.find_frameshift_from_selected_reads(selected_reads)

    @time_usage
    def get_ru_count_with_coverage_method(self, pattern_occurrences, total_counted_vntr_bp, average_coverage):
        haplotypes = 1 if self.is_haploid else 2
        estimate = [int(pattern_occurrences / (float(average_coverage) * haplotypes))] * 2
        return estimate
        pattern_occurrences = total_counted_vntr_bp / float(len(self.reference_vntr.pattern))
        read_mode = self.get_alignment_file_read_mode(alignment_file)
        samfile = pysam.AlignmentFile(alignment_file, read_mode, reference_filename=self.reference_filename)
        reference = get_reference_genome_of_alignment_file(samfile)
        bias_detector = CoverageBiasDetector(alignment_file, self.reference_vntr.chromosome, reference)
        coverage_corrector = CoverageCorrector(bias_detector.get_gc_content_coverage_map())

        logging.info('Sequencing mean coverage: %s' % coverage_corrector.get_sequencing_mean_coverage())
        observed_copy_number = pattern_occurrences / coverage_corrector.get_sequencing_mean_coverage()
        scaled_copy_number = coverage_corrector.get_scaled_coverage(self.reference_vntr, observed_copy_number)
        logging.info('scaled copy number and observed copy number: %s, %s' % (scaled_copy_number, observed_copy_number))
        return [scaled_copy_number]

    @time_usage
    def find_repeat_count_from_alignment_file(self, alignment_file, unmapped_filtered_reads, average_coverage=None,
                                              update=False):
        logging.debug('finding repeat count from alignment file for %s' % self.reference_vntr.id)

        selected_reads = self.select_illumina_reads(alignment_file, unmapped_filtered_reads, update)

        covered_repeats = []
        flanking_repeats = []
        total_counted_vntr_bp = 0
        for selected_read in selected_reads:
            repeats = get_number_of_repeats_in_vpath(selected_read.vpath)
            total_counted_vntr_bp += get_number_of_repeat_bp_matches_in_vpath(selected_read.vpath)
            logging.debug('logp of read: %s' % str(selected_read.logp))
            logging.debug('left flankign size: %s' % get_left_flanking_region_size_in_vpath(selected_read.vpath))
            logging.debug('right flanking size: %s' % get_right_flanking_region_size_in_vpath(selected_read.vpath))
            logging.debug(selected_read.sequence)
            visited_states = [state.name for idx, state in selected_read.vpath[1:-1]]
            if self.read_flanks_repeats_with_confidence(selected_read.vpath):
                logging.debug('spanning read visited states :%s' % visited_states)
                logging.debug('repeats: %s' % repeats)
                covered_repeats.append(repeats)
            else:
                flanking_repeats.append(repeats)
        flanking_repeats = sorted(flanking_repeats)
        logging.info('covered repeats: %s' % covered_repeats)
        logging.info('flanking repeats: %s' % flanking_repeats)
        min_valid_flanked = max(covered_repeats) if len(covered_repeats) > 0 else 0
        max_flanking_repeat = [r for r in flanking_repeats if r == max(flanking_repeats) and r >= min_valid_flanked]
        if len(max_flanking_repeat) < 5:
            max_flanking_repeat = []

        exact_genotype = self.find_genotype_based_on_observed_repeats(covered_repeats + max_flanking_repeat)
        if exact_genotype is not None:
            exact_genotype_log = '/'.join([str(cn) for cn in sorted(exact_genotype)])
        else:
            exact_genotype_log = 'None'
        logging.info('RU count lower bounds: %s' % exact_genotype_log)
        if average_coverage is None:
            return exact_genotype

        pattern_occurrences = sum(flanking_repeats) + sum(covered_repeats)
        return self.get_ru_count_with_coverage_method(pattern_occurrences, total_counted_vntr_bp, average_coverage)

    def find_repeat_count_from_short_reads(self, short_read_files, working_directory='./'):
        """
        Map short read sequencing data to human reference genome (hg19) and call find_repeat_count_from_alignment_file
        :param short_read_files: short read sequencing data
        :param working_directory: directory for generating the outputs
        """
        alignment_file = '' + short_read_files
        # TODO: use bowtie2 to map short reads to hg19
        return self.find_repeat_count_from_alignment_file(alignment_file, working_directory)

    @time_usage
    def train_classifier_threshold(self, reference_file, read_length=150):
        hmm = self.get_vntr_matcher_hmm(read_length=read_length)
        simulated_true_reads = self.simulate_true_reads(read_length)
        simulated_false_filtered_reads = self.simulate_false_filtered_reads(reference_file)

        processed_true_reads = self.find_hmm_score_of_simulated_reads(hmm, simulated_true_reads)
        processed_false_reads = self.find_hmm_score_of_simulated_reads(hmm, simulated_false_filtered_reads)

        recruitment_score = self.find_recruitment_score_threshold(processed_true_reads, processed_false_reads)
        return recruitment_score / float(read_length)

    @time_usage
    def find_hmm_score_of_simulated_reads(self, hmm, reads):
        initial_recruitment_score = -10000
        manager = Manager()
        processed_reads = manager.list([])
        vntr_bp_in_reads = Value('d', 0.0)
        for read_segment in reads:
            self.process_unmapped_read(None, read_segment, hmm, initial_recruitment_score, vntr_bp_in_reads, processed_reads, False)
        return processed_reads

    @time_usage
    def simulate_false_filtered_reads(self, reference_file, min_match=3):
        alphabet = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        m = 4194301

        def get_hash(string):
            result = 0
            for k in range(len(string)):
                result = (result + alphabet[string[k].upper()] * (4 ** (keyword_size - k - 1))) % m
            return result

        false_filtered_reads = []
        MAX_FALSE_READS = 10000
        read_size = 150
        keyword_size = 11
        keywords = self.get_keywords_for_filtering(True, keyword_size)
        hashed_keywords = set([get_hash(keyword) for keyword in keywords])
        match_positions = []
        vntr_start = self.reference_vntr.start_point
        vntr_end = vntr_start + self.reference_vntr.get_length()
        fasta_sequences = SeqIO.parse(open(reference_file), 'fasta')
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if name != self.reference_vntr.chromosome:
                continue
            window_hash = None
            for i in range(0, len(sequence) - keyword_size):
                if sequence[i].upper() not in 'ACTG' or sequence[i - 1 + keyword_size].upper() not in 'ACTG':
                    continue
                if window_hash is None or sequence[i - 1].upper() not in 'ACTG':
                    if 'N' in sequence[i:i + keyword_size].upper():
                        window_hash = None
                        continue
                    window_hash = get_hash(sequence[i:i + keyword_size])
                    continue
                window_hash -= alphabet[sequence[i - 1].upper()] * (4 ** (keyword_size - 1))
                window_hash = (window_hash * 4 + alphabet[sequence[i - 1 + keyword_size].upper()]) % m
                if window_hash in hashed_keywords:
                    if name == self.reference_vntr.chromosome and vntr_start - read_size < i < vntr_end:
                        continue
                    if sequence[i:i + keyword_size].upper() in keywords:
                        match_positions.append(i)
                        if len(match_positions) >= min_match and match_positions[-1] - match_positions[-min_match] < read_size:
                            for j in range(match_positions[-1] - read_size, match_positions[-min_match], 5):
                                if 'N' not in sequence[j:j + read_size].upper():
                                    false_filtered_reads.append(sequence[j:j + read_size])
                if len(false_filtered_reads) > MAX_FALSE_READS:
                    break
        return false_filtered_reads

    def simulate_true_reads(self, read_length):
        vntr = ''.join(self.reference_vntr.get_repeat_segments())
        right_flank = self.reference_vntr.right_flanking_region
        left_flank = self.reference_vntr.left_flanking_region
        locus = left_flank[-read_length:] + vntr + right_flank[:read_length]
        step_size = 1
        alphabet = ['A', 'C', 'G', 'T']
        sim_reads = []
        for i in range(0, len(locus) - read_length + 1, step_size):
            sim_reads.append(locus[i:i+read_length].upper())
        # add 4 special reads to sim_read
        for copies in range(1, len(self.reference_vntr.get_repeat_segments()) - 1):
            vntr_section = ''.join(self.reference_vntr.get_repeat_segments()[:copies])
            for i in range(1, 11):
                sim_reads.append((left_flank[-i:] + vntr_section + right_flank)[:read_length])
                sim_reads.append((left_flank + vntr_section + right_flank[:i])[-read_length:])
        min_copies = int(read_length / len(vntr)) + 1
        for i in range(1, 21):
            # print(len((vntr * min_copies)[i:read_length+i]))
            sim_reads.append((vntr * min_copies)[i:read_length+i])
            # print(len((vntr * min_copies)[-read_length-i:-i]))
            sim_reads.append((vntr * min_copies)[-read_length-i:-i])
        simulated_true_reads = []
        for sim_read in sim_reads:
            from random import randint
            for i in range(randint(1, 2)):
                temp_read = list(sim_read)
                temp_read[randint(0, len(sim_read)-1)] = alphabet[randint(0, 3)]
                sim_read = ''.join(temp_read)
            simulated_true_reads.append(sim_read)
        return simulated_true_reads

    @time_usage
    def find_recruitment_score_threshold(self, processed_true_reads, processed_false_reads):
        from sklearn.linear_model import LogisticRegression
        true_scores = [read.logp for read in processed_true_reads]
        false_scores = [read.logp for read in processed_false_reads]
        if len(false_scores) == 0:
            false_scores = [min(true_scores) - 2]
        clf = LogisticRegression()
        x = [[score] for score in true_scores + false_scores]
        y = [1] * len(true_scores) + [0] * len(false_scores)
        clf.fit(x, y)
        recruitment_score = max(true_scores)
        for i in range(-1, -300, -1):
            if int(clf.predict([[i]])) == 0:
                recruitment_score = i
                break
        return recruitment_score
