from math import log
from Bio import SeqIO
from repeat_finder import get_blast_matched_ids
from hmm_utils import *
from sam_utils import get_related_reads_and_read_count_in_samfile, get_VNTR_coverage_over_total_coverage


class VNTRFinder:
    def __init__(self, id, pattern, ref_start_pos, ref_repeat_count, ref_visited_states, ref_file_name='chr15.fa'):
        self.id = id
        self.pattern = pattern
        self.ref_start_pos = ref_start_pos
        self.ref_repeat_count = ref_repeat_count
        self.ref_visited_states = ref_visited_states
        self.reference_file_name = ref_file_name

    def get_VNTR_matcher_hmm(self, patterns, copies, left_flanking_region, right_flanking_region):
        left_flanking_matcher = get_suffix_matcher_hmm(left_flanking_region)
        right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
        repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
        left_flanking_matcher.concatenate(repeats_matcher)
        left_flanking_matcher.concatenate(right_flanking_matcher)
        left_flanking_matcher.bake(merge=None)
        return left_flanking_matcher

    def extract_repeat_segments_from_visited_states(self, ref_file_name='chr15.fa'):
        fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
        ref_sequence = ''
        for fasta in fasta_sequences:
            name, ref_sequence = fasta.id, str(fasta.seq)
        pattern_start = self.ref_start_pos
        corresponding_region_in_ref = ref_sequence[pattern_start:pattern_start + (len(self.pattern) + 5) * self.ref_repeat_count].upper()

        lengths = []
        prev_start = None
        for i in range(len(self.ref_visited_states)):
            if self.ref_visited_states[i].startswith('unit_end') and prev_start is not None:
                current_len = 0
                for j in range(prev_start, i):
                    if is_matching_state(self.ref_visited_states[j]):
                        current_len += 1
                lengths.append(current_len)
            if self.ref_visited_states[i].startswith('unit_start'):
                prev_start = i

        repeat_segments = []
        added = 0
        for l in lengths:
            repeat_segments.append(corresponding_region_in_ref[added:added+l])
            added += l
        return repeat_segments

    def get_flanking_regions(self, start_point, end_point, flanking_region_size):
        fasta_sequences = SeqIO.parse(open(self.ref_file_name), 'fasta')
        ref_sequence = ''
        for fasta in fasta_sequences:
            name, ref_sequence = fasta.id, str(fasta.seq)
        left_flanking = ref_sequence[start_point - flanking_region_size:start_point].upper()
        right_flanking = ref_sequence[end_point:end_point + flanking_region_size].upper()
        return left_flanking, right_flanking

    def find_repeat_count(self, read_files):
        repeat_segments = self.extract_repeat_segments_from_visited_states()
        copies = int(round(150.0 / len(self.pattern) + 0.5))
        flanking_region_size = 150 - 1
        end_point = self.ref_start_pos + sum([len(e) for e in repeat_segments])
        left_flanking_region, right_flanking_region = self.get_flanking_regions(self.ref_start_pos, end_point, flanking_region_size)
        hmm = self.get_VNTR_matcher_hmm(repeat_segments * 100, copies, left_flanking_region, right_flanking_region)
        total_occurrences = 0

        word_size = int(len(self.pattern)/3)
        if word_size > 11:
            word_size = 11
        word_size = str(word_size)
        blast_ids = set([])
        for repeat_segment in repeat_segments:
            blast_ids |= get_blast_matched_ids(repeat_segment, 'original_reads/original_reads', max_seq='50000',
                                               evalue=1, word_size=word_size, search_id=str(self.id))

        print('blast selected ', len(blast_ids), ' reads')
        if len(blast_ids) == 50 * 1000:
            with open('errors.txt', 'a') as out:
                out.write('maximum number of read selected in filtering for pattern %s\n' % self.id)

        samfile = 'original_reads/paired_dat.sam'
        related_reads, read_count = get_related_reads_and_read_count_in_samfile(self.pattern, self.ref_start_pos, read_file=samfile, pattern_end=end_point)
        for re_read in related_reads:
            if re_read not in blast_ids:
                print('FN in filtering')

        min_score = 0
        for seg in repeat_segments:
            min_score = min(min_score, hmm.viterbi((seg * copies)[:150])[0])
            print(min_score)
        min_score_of_single_unit = min_score / copies - 5
        different_read_score_reads = {}
        different_flanking_score_reads = {}
        different_read_score_occurrences = {}
        for i in range(-13, 13):
            different_read_score_occurrences[int(min_score) + i * 8] = 0
        different_flanking_score_occurrences = {log(0.8): 0, log(0.7): 0, log(0.6): 0, log(0.5): 0, log(0.4): 0}
        print('different_read_score_occurrences: ', different_read_score_occurrences)
        print('different_flanking_score_occurrences: ', different_flanking_score_occurrences)

        number_of_reads = 0
        read_length = 0
        total_length = 102531392
        for read_file in read_files:
            print('opening read file')
            reads = SeqIO.parse(read_file, 'fasta')
            for read_segment in reads:
                if number_of_reads == 0:
                    read_length = len(str(read_segment.seq))
                number_of_reads += 1
                if read_segment.id not in blast_ids and read_segment.id not in related_reads:
                    continue
                logp, vpath = hmm.viterbi(str(read_segment.seq))
                rev_logp, rev_vpath = hmm.viterbi(str(read_segment.seq.reverse_complement()))
                if logp < rev_logp:
                    logp = rev_logp
                    vpath = rev_vpath
                occurrence = get_number_of_matched_repeats_in_a_vpath(vpath)
                min_occ_to_add_read = 1
                if len(self.pattern) < 50:
                    min_occ_to_add_read = 2
                if occurrence >= min_occ_to_add_read:
                    min_occ_to_count = min_occ_to_add_read
                    if len(self.pattern) < 24:
                        min_occ_to_count = 3
                    for s_threshold in different_read_score_occurrences.keys():
                        if logp > s_threshold:
                            different_read_score_occurrences[s_threshold] += occurrence if occurrence > min_occ_to_count else 0
                            if s_threshold not in different_read_score_reads.keys():
                                different_read_score_reads[s_threshold] = []
                            different_read_score_reads[s_threshold].append(read_segment.id)
                    for s_threshold in different_flanking_score_occurrences.keys():
                        if logp > min_score_of_single_unit * occurrence + s_threshold * (150-occurrence*len(self.pattern)):
                            different_flanking_score_occurrences[s_threshold] += occurrence if occurrence > min_occ_to_count else 0
                            if s_threshold not in different_flanking_score_reads.keys():
                                different_flanking_score_reads[s_threshold] = []
                            different_flanking_score_reads[s_threshold].append(read_segment.id)
                total_occurrences += occurrence

                number_of_reads += 1

        avg_coverage = float(number_of_reads * read_length) / total_length

        cn = 10000
        min_error = 1000
        for s_threshold in different_read_score_reads.keys():
            selected_reads = different_read_score_reads[s_threshold]
            TP = [read for read in selected_reads if read in related_reads]
            FP = [read for read in selected_reads if read not in TP]
            FN = [read for read in related_reads if read not in selected_reads]
            # print('TP:', len(TP), 'FP:', len(FP), 'selected:', len(selected_reads))
            # print('FN:', len(FN))
            sensitivity = float(len(TP)) / len(related_reads) if len(related_reads) > 0 else 0
            if sensitivity > 0.85 and sensitivity < 1 and len(FN) > 0 and len(FP) > 0:
                print('sensitivity ', sensitivity, ' FN:', FN[0], ' FP:', FP[0])
            with open('FP_and_sensitivity_HMM_read_scoring_method.txt', 'a') as outfile:
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (len(FP), sensitivity, s_threshold, self.id, len(self.pattern), len(TP)))
            occurrences = different_read_score_occurrences[s_threshold]
            if sensitivity > 0.6 and abs(self.ref_repeat_count - occurrences / avg_coverage) < min_error:
                min_error = abs(self.ref_repeat_count - occurrences / avg_coverage)
                cn = occurrences / avg_coverage
        for s_threshold in different_flanking_score_reads.keys():
            selected_reads = different_flanking_score_reads[s_threshold]
            TP = [read for read in selected_reads if read in related_reads]
            FP = [read for read in selected_reads if read not in TP]
            sensitivity = float(len(TP)) / len(related_reads) if len(related_reads) > 0 else 0
            with open('FP_and_sensitivity_HMM_joint_scoring_method.txt', 'a') as outfile:
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (len(FP), sensitivity, s_threshold, self.id, len(self.pattern), len(TP)))
            occurrences = different_flanking_score_occurrences[s_threshold]
            if sensitivity > 0.6 and abs(self.ref_repeat_count - occurrences / avg_coverage) < min_error:
                min_error = abs(self.ref_repeat_count - occurrences / avg_coverage)
                cn = occurrences / avg_coverage

        return cn


with open('patterns.txt') as input:
    patterns = input.readlines()
    patterns = [pattern.strip() for pattern in patterns]
with open('start_points.txt') as input:
    lines = input.readlines()
    start_points = [int(num.strip())-1 for num in lines]
with open('pattern_repeat_counts.txt') as input:
    lines = input.readlines()
    repeat_counts = [int(num.strip()) for num in lines]
with open('visited_states.txt') as input:
    lines = input.readlines()
    visited_states_list = [states.strip().split() for states in lines]

read_files = ['original_reads/paired_dat1.fasta', 'original_reads/paired_dat2.fasta']
for i in range(len(patterns)):
    print(i)
    if repeat_counts[i] == 0:
        continue
    vntr_fidner = VNTRFinder(i+1, patterns[i], start_points[i], repeat_counts[i], visited_states_list[i])
    cn = vntr_fidner.find_repeat_count(read_files)
    with open('hmm_repeat_count.txt', 'a') as output:
        output.write('%s %s\n' % (i, cn / repeat_counts[i]))
    # repeat_segments = extract_repeat_segments_from_visited_states(patterns[i], start_points[i], repeat_counts[i], visited_states_list[i])
    # end_point = start_points[i] + sum([len(e) for e in repeat_segments])
    # VNTR_coverage_ratio = get_VNTR_coverage_over_total_coverage(start_points[i], end_point)
    # with open('vntr_coverage_ratio.txt', 'a') as output:
    #     output.write('%s %s\n' % (i, VNTR_coverage_ratio))
