from Bio import SeqIO, pairwise2
from hmm_utils import build_reference_repeat_finder_hmm, get_repeat_segments_from_visited_states_and_region
from settings import *


class ReferenceVNTR:
    def __init__(self, vntr_id, pattern, start_point, chromosome, estimated_repeats=None):
        self.non_overlapping = True
        self.has_homologous = False
        self.id = vntr_id
        self.pattern = pattern
        self.start_point = start_point
        self.chromosome = chromosome
        self.estimated_repeats = estimated_repeats
        self.repeat_segments = None
        self.left_flanking_region = None
        self.right_flanking_region = None

    def init_from_vntrseek_data(self):
        corresponding_region_in_ref = self.get_corresponding_region_in_ref()
        repeat_segments = self.find_repeat_segments(corresponding_region_in_ref)
        self.repeat_segments = repeat_segments
        flanking_region_size = 150 - 10
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
        structure1 = self.left_flanking_region[120:] + self.pattern + self.right_flanking_region[:20]
        structure2 = another.left_flanking_region[120:] + another.pattern + another.right_flanking_region[:20]
        alignment_score = pairwise2.align.localms(structure1, structure2, 1, -1, -1, -1, score_only=True)
        if float(alignment_score) / len(structure1) > 0.8 or float(alignment_score) / len(structure2) > 0.8:
            return True
        return False

    def find_repeat_segments(self, region_in_ref):
        patterns = [self.pattern]
        model = build_reference_repeat_finder_hmm(patterns, copies=self.estimated_repeats)
        logp, path = model.viterbi(region_in_ref)
        visited_states = [state.name for idx, state in path[1:-1]]
        repeat_segments = get_repeat_segments_from_visited_states_and_region(visited_states, region_in_ref)

        return repeat_segments

    def __get_reference_chromosome_sequence(self):
        ref_file_name = HG19_DIR + self.chromosome + '.fa'
        fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
        ref_sequence = ''
        for fasta in fasta_sequences:
            name, ref_sequence = fasta.id, str(fasta.seq)
        return ref_sequence

    def get_corresponding_region_in_ref(self):
        ref_sequence = self.__get_reference_chromosome_sequence()
        estimated_length = len(self.pattern) * self.estimated_repeats
        corresponding_region_in_ref = ref_sequence[self.start_point:self.start_point + estimated_length].upper()
        return corresponding_region_in_ref

    def get_flanking_regions(self, flanking_region_size=140):
        ref_sequence = self.__get_reference_chromosome_sequence()
        left_flanking = ref_sequence[self.start_point - flanking_region_size:self.start_point].upper()
        end_of_repeats = self.start_point + self.get_length()
        right_flanking = ref_sequence[end_of_repeats:end_of_repeats + flanking_region_size].upper()
        return left_flanking, right_flanking


def load_unprocessed_vntrseek_data(vntrseek_output, chr=None):
    vntrs = []
    with open(vntrseek_output) as input_file:
        input_lines = [line.strip() for line in input_file.readlines() if line.strip() != '']
        for vntr_id, line in enumerate(input_lines):
            vntrseek_repeat, _, pattern, chromosome, start = line.split()
            estimated_repeats = int(float(vntrseek_repeat) + 5)
            if chr and chromosome != chr:
                continue
            vntrs.append(ReferenceVNTR(vntr_id, pattern, int(start)-1, chromosome, estimated_repeats))
    return vntrs


def find_non_overlapping_vntrs(vntrseek_output='repeats_length_patterns_chromosomes_starts.txt'):
    vntrs = load_unprocessed_vntrseek_data(vntrseek_output)
    skipped_vntrs = []
    for i in range(len(vntrs)):
        print(i)
        estimated_end_point = len(vntrs[i].pattern) * vntrs[i].estimated_repeats + vntrs[i].start_point
        if i < len(vntrs) - 1 and vntrs[i].chromosome == vntrs[i+1].chromosome and estimated_end_point > vntrs[i+1].start_point:
            vntrs[i].estimated_repeats += vntrs[i+1].estimated_repeats
        vntrs[i].init_from_vntrseek_data()
        repeat_segments = vntrs[i].get_repeat_segments()
        if i in skipped_vntrs:
            vntrs[i].non_overlapping = False
        else:
            j = i + 1
            end_point = len(vntrs[i].pattern) * len(repeat_segments) + vntrs[i].start_point
            while j < len(vntrs) and vntrs[i].chromosome == vntrs[j].chromosome and end_point > vntrs[j].start_point:
                skipped_vntrs.append(j)
                j += 1
    return vntrs


def identify_homologous_vntrs(vntrs, chr=None):
    for i in range(len(vntrs)):
        for j in range(i + 1, len(vntrs)):
            if chr and (chr != vntrs[i].chromosome or chr != vntrs[j].chromosome):
                continue
            if vntrs[i].is_homologous_vntr(vntrs[j]):
                vntrs[i].has_homologous = True
                vntrs[j].has_homologous = True
    return vntrs


def process_vntrseek_data():
    vntrs = find_non_overlapping_vntrs()
    for vntr in vntrs:
        comma_separated_segments = ','.join(vntr.get_repeat_segments())
        with open('repeats_and_segments.txt', 'a') as out:
            out.write('%s %s %s %s %s\n' % (vntr.id, vntr.is_non_overlapping(), vntr.left_flanking_region,
                                            vntr.right_flanking_region, comma_separated_segments))


def load_processed_vntrs_data(vntrseek_output='repeats_length_patterns_chromosomes_starts.txt'):
    vntrs = []
    with open(vntrseek_output) as input_file:
        vntrseek_data = [line.strip() for line in input_file.readlines() if line.strip() != '']
    with open('repeats_and_segments.txt') as input_file:
        segments_lines = input_file.readlines()
    for vntr_id in range(len(vntrseek_data)):
        _id, non_overlapping, left_flanking_region, right_flanking_region, segments = segments_lines[vntr_id].split()
        repeat_segments = segments.split(',')
        vntrseek_repeat, _, pattern, chromosome, start = vntrseek_data[vntr_id].split()
        vntr = ReferenceVNTR(vntr_id, pattern, int(start)-1, chromosome, vntrseek_repeat)
        vntr.init_from_xml(repeat_segments, left_flanking_region, right_flanking_region)
        vntr.non_overlapping = True if non_overlapping == 'True' else False
        vntrs.append(vntr)
    return vntrs


# process_vntrseek_data()
