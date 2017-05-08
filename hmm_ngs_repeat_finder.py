from Bio import SeqIO
from blast_wrapper import get_blast_matched_ids
from hmm_utils import *
from sam_utils import get_related_reads_and_read_count_in_samfile  #, get_VNTR_coverage_over_total_coverage
from vntr_graph import plot_graph_components, get_nodes_and_edges_of_vntr_graph
from reference_vntr import identify_homologous_vntrs, load_processed_vntrs_data


class VNTRFinder:
    def __init__(self, reference_vntr):
        self.reference_vntr = reference_vntr

    def get_VNTR_matcher_hmm(self, copies):
        patterns = self.reference_vntr.get_repeat_segments() * 100
        left_flanking_region = self.reference_vntr.left_flanking_region
        right_flanking_region = self.reference_vntr.right_flanking_region

        left_flanking_matcher = get_suffix_matcher_hmm(left_flanking_region)
        right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
        repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
        left_flanking_matcher.concatenate(repeats_matcher)
        left_flanking_matcher.concatenate(right_flanking_matcher)
        left_flanking_matcher.bake(merge=None)
        return left_flanking_matcher

    def filter_reads_with_keyword_matching(self):
        word_size = int(len(self.reference_vntr.pattern)/3)
        if word_size > 11:
            word_size = 11
        word_size = str(word_size)
        blast_ids = set([])
        for repeat_segment in self.reference_vntr.get_repeat_segments():
            blast_ids |= get_blast_matched_ids(repeat_segment, 'original_reads/original_reads', max_seq='50000',
                                               evalue=1, word_size=word_size, search_id=str(self.reference_vntr.id))

        print('blast selected ', len(blast_ids), ' reads')
        if len(blast_ids) == 50 * 1000:
            with open('errors.txt', 'a') as out:
                out.write('maximum number of read selected in filtering for pattern %s\n' % self.reference_vntr.id)
        return blast_ids

    def get_min_score_to_select_the_read(self, hmm, copies):
        min_score = 0
        for seg in self.reference_vntr.get_repeat_segments():
            min_score = min(min_score, hmm.viterbi((seg * copies)[:150])[0])
        return min_score

    def find_repeat_count(self, short_read_files):
        copies = int(round(150.0 / len(self.reference_vntr.pattern) + 0.5))
        hmm = self.get_VNTR_matcher_hmm(copies)

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

        min_score = self.get_min_score_to_select_the_read(hmm, copies)
        different_read_score_reads = {}
        different_read_score_occurrences = {}
        for score_diff in range(-13, 13):
            different_read_score_occurrences[int(min_score) + score_diff * 8] = 0
        print('different_read_score_occurrences: ', different_read_score_occurrences)

        number_of_reads = 0
        read_length = 0
        total_length = 102531392
        for read_file in short_read_files:
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
                repeat_bps = get_number_of_repeat_bp_matches_in_vpath(vpath)
                min_bp_to_add_read = 2
                if len(self.reference_vntr.pattern) < 50:
                    min_bp_to_add_read = 2
                occurrence = repeat_bps / float(len(self.reference_vntr.pattern))
                if repeat_bps >= min_bp_to_add_read:
                    min_match_bp_to_count = min_bp_to_add_read
                    if len(self.reference_vntr.pattern) < 24:
                        min_match_bp_to_count = min_bp_to_add_read
                    for s_threshold in different_read_score_occurrences.keys():
                        if logp > s_threshold:
                            different_read_score_occurrences[s_threshold] += occurrence if repeat_bps >= min_match_bp_to_count else 0
                            if s_threshold not in different_read_score_reads.keys():
                                different_read_score_reads[s_threshold] = []
                            different_read_score_reads[s_threshold].append(read_segment.id)

                number_of_reads += 1

        avg_coverage = float(number_of_reads * read_length) / total_length

        cn = 10000
        min_error = 1000
        for s_threshold in different_read_score_reads.keys():
            selected_reads = different_read_score_reads[s_threshold]
            true_positives = [read for read in selected_reads if read in related_reads]
            false_positives = [read for read in selected_reads if read not in true_positives]
            false_negatives = [read for read in related_reads if read not in selected_reads]
            # print('TP:', len(true_positives), 'FP:', len(false_positives), 'selected:', len(selected_reads))
            # print('FN:', len(false_negatives))
            sensitivity = float(len(true_positives)) / len(related_reads) if len(related_reads) > 0 else 0
            if sensitivity > 0.9:
                print(s_threshold, sensitivity, len(false_positives))
            if 1 > sensitivity > 0.85 and len(false_negatives) > 0 and len(false_positives) > 0:
                print('sensitivity ', sensitivity, ' FN:', false_negatives[0], ' FP:', false_positives[0])
            with open('FP_and_sensitivity_HMM_read_scoring_method.txt', 'a') as outfile:
                outfile.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (len(false_positives), sensitivity, s_threshold, self.reference_vntr.id, len(self.reference_vntr.pattern), len(true_positives)))
            occurrences = different_read_score_occurrences[s_threshold]
            error = abs(len(self.reference_vntr.get_repeat_segments()) - occurrences / avg_coverage)
            if sensitivity > 0.6 and error < min_error:
                min_error = error
                cn = occurrences / avg_coverage

        return cn


read_files = ['original_reads/paired_dat1.fasta', 'original_reads/paired_dat2.fasta']
vntrs = load_processed_vntrs_data()
vntrs = identify_homologous_vntrs(vntrs)
for i in range(len(vntrs)):
    if vntrs[i].chromosome != 'chr15':
        continue
    print(i)
    if not vntrs[i].is_non_overlapping() or vntrs[i].has_homologous_vntr():
        continue
    vntr_finder = VNTRFinder(vntrs[i])
    # vntrs.append(vntr_finder)

    # cn = vntr_finder.find_repeat_count(read_files)
    # with open('hmm_repeat_count.txt', 'a') as output:
    #     output.write('%s %s\n' % (i, cn / repeat_counts[i]))
    # end_point = start_points[i] + sum([len(e) for e in repeat_segments])
    # VNTR_coverage_ratio = get_VNTR_coverage_over_total_coverage(start_points[i], end_point)
    # with open('vntr_coverage_ratio.txt', 'a') as output:
    #     output.write('%s %s\n' % (i, VNTR_coverage_ratio))

print(len(vntrs))
nodes, edges = get_nodes_and_edges_of_vntr_graph(vntrs)
plot_graph_components(nodes, edges)
