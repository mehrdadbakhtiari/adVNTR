from Bio import SeqIO
from blast_wrapper import get_blast_matched_ids, make_blast_database
from hmm_utils import *
from sam_utils import get_related_reads_and_read_count_in_samfile  #, get_VNTR_coverage_over_total_coverage
import settings
# from vntr_graph import plot_graph_components, get_nodes_and_edges_of_vntr_graph
from reference_vntr import identify_homologous_vntrs, load_processed_vntrs_data

import os


class VNTRFinder:
    """Find the VNTR structure of a reference VNTR in NGS data of the donor."""

    def __init__(self, reference_vntr):
        self.reference_vntr = reference_vntr

    def build_vntr_matcher_hmm(self, copies):
        patterns = self.reference_vntr.get_repeat_segments() * 100
        flanking_region_size = 140
        left_flanking_region = self.reference_vntr.left_flanking_region[-flanking_region_size:]
        right_flanking_region = self.reference_vntr.right_flanking_region[:flanking_region_size]

        vntr_matcher = get_suffix_matcher_hmm(left_flanking_region)
        right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
        repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
        vntr_matcher.concatenate(repeats_matcher)
        vntr_matcher.concatenate(right_flanking_matcher)
        vntr_matcher.bake(merge=None)
        return vntr_matcher

    def get_vntr_matcher_hmm(self, copies):
        """Try to load trained HMM for this VNTR
        If there was no trained HMM, it will build one and store it for later usage
        """
        stored_hmm_file = settings.TRAINED_HMMS_DIR + str(self.reference_vntr.id) + '.json'
        if settings.USE_TRAINED_HMMS and os.path.isfile(stored_hmm_file):
            model = Model()
            model = model.from_json(stored_hmm_file)
            return model

        vntr_matcher = self.build_vntr_matcher_hmm(copies)

        json_str = vntr_matcher.to_json()
        with open(stored_hmm_file, 'w') as outfile:
            outfile.write(json_str)
        return vntr_matcher

    def filter_reads_with_keyword_matching(self, short_read_files):
        word_size = int(len(self.reference_vntr.pattern)/3)
        if word_size > 11:
            word_size = 11
        word_size = str(word_size)
        blast_ids = set([])

        dir = '/'.join(short_read_files[0].split('/')[:-1]) + '/'
        db_name = dir[:-1]
        blast_db_name = dir + db_name
        if not os.path.exists(blast_db_name + '.nal'):
            make_blast_database(short_read_files, blast_db_name)
        for repeat_segment in self.reference_vntr.get_repeat_segments():
            blast_ids |= get_blast_matched_ids(repeat_segment, blast_db_name, max_seq='50000',
                                               evalue=10, word_size=word_size, search_id=str(self.reference_vntr.id))

        print('blast selected ', len(blast_ids), ' reads')
        if len(blast_ids) == len(self.reference_vntr.get_repeat_segments()) * 50 * 1000:
            with open('errors.txt', 'a') as out:
                out.write('maximum number of read selected in filtering for pattern %s\n' % self.reference_vntr.id)
        return blast_ids

    def calculate_min_score_to_select_a_read(self, hmm, true_reads, update=True):
        """Calculate the score distribution of false positive reads
        and return score to select the 0.0001 percentile of the distribution

        It will update the score for this VNTR in precomputed data if update variable is True, otherwise, it will
        only writes the result if there is no score for this VNTR in precomputed data.
        """
        with open('id_score_to_select.txt', 'r') as infile:
            id_score_map = {int(vntr_id): float(score) for vntr_id, score in infile.readlines()}
        # TODO
        with open('id_score_to_select.txt', 'w') as outfile:
            for vntr_id, score in id_score_map.items():
                outfile.write('%s %s\n' % (vntr_id, score))
        return 0

    def get_min_score_to_select_a_read(self):
        """Try to load the minimum score for this VNTR
        If the score was not precomputed, it writes an error and returns 0
        """
        with open('id_score_to_select.txt', 'r') as infile:
            id_score_pairs = [(line.split()[0], line.split()[1]) for line in infile.readlines() if line.strip() != '']
            id_score_map = {int(vntr_id): float(score) for vntr_id, score in id_score_pairs}

        if self.reference_vntr.id not in id_score_map:
            print('Minimum score is not precomputed for vntr id: %s' % self.reference_vntr.id)
            id_score_map[reference_vntrs.id] = 0

        return id_score_map[self.reference_vntr.id]

    def find_repeat_count(self, short_read_files):
        read_length = 150
        copies = int(round(float(read_length) / len(self.reference_vntr.pattern) + 0.5))
        hmm = self.get_vntr_matcher_hmm(copies)

        filtered_read_ids = self.filter_reads_with_keyword_matching(short_read_files)

        min_score_to_count_read = self.get_min_score_to_select_a_read()
        selected_reads = []
        counted_vntr_base_pairs = 0

        number_of_reads = 0
        read_length = 0
        min_repeat_bp_to_add_read = 2
        if len(self.reference_vntr.pattern) < 30:
            min_repeat_bp_to_add_read = 2
        min_repeat_bp_to_count_repeats = 2
        for read_file in short_read_files:
            print('opening a read file')
            reads = SeqIO.parse(read_file, 'fasta')
            for read_segment in reads:
                if number_of_reads == 0:
                    read_length = len(str(read_segment.seq))
                number_of_reads += 1
                if read_segment.id not in filtered_read_ids:
                    continue
                logp, vpath = hmm.viterbi(str(read_segment.seq))
                rev_logp, rev_vpath = hmm.viterbi(str(read_segment.seq.reverse_complement()))
                if logp < rev_logp:
                    logp = rev_logp
                    vpath = rev_vpath
                repeat_bps = get_number_of_repeat_bp_matches_in_vpath(vpath)
                if logp > min_score_to_count_read:
                    if repeat_bps > min_repeat_bp_to_count_repeats:
                        counted_vntr_base_pairs += repeat_bps
                    if repeat_bps > min_repeat_bp_to_add_read:
                        selected_reads.append(read_segment.id)

                number_of_reads += 1

        avg_coverage = float(number_of_reads * read_length) / settings.GENOME_LENGTH
        pattern_occurrences = counted_vntr_base_pairs / float(len(self.reference_vntr.pattern))
        print('pattern_occurrences and avg_coverage: ', pattern_occurrences, avg_coverage)
        return pattern_occurrences / avg_coverage

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
            outfile.write('%s\t%s\t%s\t%s\t%s\n' % (len(false_positives), sensitivity, self.reference_vntr.id, len(self.reference_vntr.pattern), len(true_positives)))
        error = abs(len(self.reference_vntr.get_repeat_segments()) - occurrences / avg_coverage)
        print(error)


read_files = ['original_reads/paired_dat1.fasta', 'original_reads/paired_dat2.fasta']
reference_vntrs = load_processed_vntrs_data()
reference_vntrs = identify_homologous_vntrs(reference_vntrs, 'chr15')
accurate_vntr_list = [271, 281, 283, 287, 288, 325, 327, 328, 329]
for i in range(len(reference_vntrs)):
    if reference_vntrs[i].chromosome != 'chr15':
        continue
    print(i)
    if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
        continue
    if reference_vntrs[i].id not in accurate_vntr_list:
        continue
    vntr_finder = VNTRFinder(reference_vntrs[i])
    copy_number = vntr_finder.find_repeat_count(read_files)
    # vntr_finder.find_accuracy()

    with open('hmm_repeat_count.txt', 'a') as output:
        output.write('%s %s\n' % (i, copy_number / len(reference_vntrs[i].get_repeat_segments())))
    # end_point = start_points[i] + sum([len(e) for e in repeat_segments])
    # VNTR_coverage_ratio = get_VNTR_coverage_over_total_coverage(start_points[i], end_point)
    # with open('vntr_coverage_ratio.txt', 'a') as output:
    #     output.write('%s %s\n' % (i, VNTR_coverage_ratio))

# print(len(reference_vntrs))
# nodes, edges = get_nodes_and_edges_of_vntr_graph(reference_vntrs)
# plot_graph_components(nodes, edges)
