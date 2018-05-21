from Bio import pairwise2

from advntr.hmm_utils import build_reference_repeat_finder_hmm, get_repeat_segments_from_visited_states_and_region
from advntr.utils import get_chromosome_reference_sequence


class ReferenceVNTR:
    def __init__(self, vntr_id, pattern, start_point, chromosome, gene_name, annotation, estimated_repeats=None,
                 chromosome_sequence=None, scaled_score=0):
        self.non_overlapping = True
        self.has_homologous = False
        self.id = vntr_id
        self.pattern = pattern
        self.start_point = start_point
        self.chromosome = chromosome
        self.gene_name = gene_name
        self.annotation = annotation
        self.estimated_repeats = estimated_repeats
        self.repeat_segments = []
        self.left_flanking_region = None
        self.right_flanking_region = None
        self.chromosome_sequence = chromosome_sequence
        self.scaled_score = scaled_score

    def init_from_vntrseek_data(self):
        corresponding_region_in_ref = self.get_corresponding_region_in_ref()
        repeat_segments = self.find_repeat_segments(corresponding_region_in_ref)
        self.repeat_segments = repeat_segments
        flanking_region_size = 500
        self.left_flanking_region, self.right_flanking_region = self.get_flanking_regions(flanking_region_size)

    def init_from_xml(self, repeat_segments, left_flanking_region, right_flanking_region):
        self.repeat_segments = repeat_segments
        self.left_flanking_region = left_flanking_region
        self.right_flanking_region = right_flanking_region

    def is_non_overlapping(self):
        return self.non_overlapping

    def has_homologous_vntr(self):
        return self.has_homologous

    def get_length(self):
        return sum([len(e) for e in self.repeat_segments])

    def get_repeat_segments(self):
        return self.repeat_segments

    def is_homologous_vntr(self, another):
        structure1 = self.left_flanking_region[-20:] + self.pattern + self.right_flanking_region[:20]
        structure2 = another.left_flanking_region[-20:] + another.pattern + another.right_flanking_region[:20]
        alignment_score = pairwise2.align.localms(structure1, structure2, 1, -1, -1, -1, score_only=True)
        if float(alignment_score) / len(structure1) > 0.66 or float(alignment_score) / len(structure2) > 0.66:
            return True
        return False

    def find_repeat_segments(self, region_in_ref):
        patterns = [self.pattern]
        model = build_reference_repeat_finder_hmm(patterns, copies=self.estimated_repeats)
        logp, path = model.viterbi(region_in_ref)
        visited_states = [state.name for idx, state in path[1:-1]]
        repeat_segments = get_repeat_segments_from_visited_states_and_region(visited_states, region_in_ref)

        return repeat_segments

    def __get_chromosome_reference_sequence(self):
        if self.chromosome_sequence is not None:
            return self.chromosome_sequence
        return get_chromosome_reference_sequence(self.chromosome)

    def get_corresponding_region_in_ref(self):
        ref_sequence = self.__get_chromosome_reference_sequence()
        estimated_length = len(self.pattern) * self.estimated_repeats
        corresponding_region_in_ref = ref_sequence[self.start_point:self.start_point + estimated_length].upper()
        while corresponding_region_in_ref.find('N') != -1:
            n_index = corresponding_region_in_ref.find('N')
            corresponding_region_in_ref = corresponding_region_in_ref[:n_index]
        return corresponding_region_in_ref

    def get_flanking_regions(self, flanking_region_size=140):
        ref_sequence = self.__get_chromosome_reference_sequence()
        left_flanking = ref_sequence[self.start_point - flanking_region_size:self.start_point].upper()
        end_of_repeats = self.start_point + self.get_length()
        right_flanking = ref_sequence[end_of_repeats:end_of_repeats + flanking_region_size].upper()
        return left_flanking, right_flanking
