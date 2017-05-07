from Bio import SeqIO, pairwise2
from hmm_utils import build_reference_repeat_finder_hmm, get_repeat_segments_from_visited_states_and_region


class ReferenceVNTR:
    def __init__(self, id, pattern, start_point, estimated_repeats=None, chromosome='chr15'):
        self.non_overlapping = True
        self.has_homologous = False
        self.id = id
        self.pattern = pattern
        self.start_point = start_point
        self.estimated_repeats = estimated_repeats
        self.chromosome = chromosome
        self.repeat_segments = None

    def init_from_vntrseek_data(self):
        corresponding_region_in_ref = self.get_corresponding_region_in_ref()
        repeat_segments = self.find_repeat_segments(corresponding_region_in_ref)
        self.repeat_segments = repeat_segments

    def init_from_xml(self, repeat_segments):
        self.repeat_segments = repeat_segments

    def is_non_overlapping(self):
        return self.non_overlapping

    def has_homologous_vntr(self):
        return self.has_homologous

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

    def get_corresponding_region_in_ref(self):
        ref_file_name = self.chromosome + '.fa'
        fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
        ref_sequence = ''
        for fasta in fasta_sequences:
            name, ref_sequence = fasta.id, str(fasta.seq)
        estimated_length = len(self.pattern) * self.estimated_repeats
        corresponding_region_in_ref = ref_sequence[self.start_point:self.start_point + estimated_length].upper()
        return corresponding_region_in_ref


def find_non_overlapping_vntrs(vntrseek_output='repeats_length_patterns_chromosomes_starts.txt'):
    vntrs = []
    with open (vntrseek_output) as input:
        input_lines = [line.strip() for line in input.readlines() if line.strip() != '']
        for vntr_id, line in enumerate(input_lines):
            vntrseek_repeat, _, pattern, chr, start = line.split()
            estimated_repeats = int(vntrseek_repeat) + 5
            vntrs.append(ReferenceVNTR(vntr_id, pattern, int(start)-1, estimated_repeats, chr))

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


def identify_homologous_vntrs(vntrs):
    for i in range(len(vntrs)):
        for j in range(i + 1, len(vntrs)):
            if vntrs[i].is_homologous_vntr(vntrs[j]):
                vntrs[i].has_homologous = True
                vntrs[j].has_homologous = True
    return vntrs


def process_vntrseek_data():
    vntrs = find_non_overlapping_vntrs()
    for vntr in vntrs:
        comma_separated_segments = ','.join(vntr.get_repeat_segments())
        with open('repeats_and_segments.txt', 'a') as out:
            out.write('%s %s %s' % (vntr.id, vntr.is_non_overlapping(), comma_separated_segments))


def load_non_overlapping_vntrs_data(vntrseek_output='repeats_length_patterns_chromosomes_starts.txt'):
    vntrs = []
    with open (vntrseek_output) as input:
        vntrseek_data = [line.strip() for line in input.readlines() if line.strip() != '']
    with open('repeats_and_segments.txt') as input:
        segments_lines = input.readlines()
    for vntr_id in range(len(vntrseek_data)):
        id, repeats, segments = segments_lines[vntr_id].split()
        repeat_segments = segments.split(',')
        vntrseek_repeat, _, pattern, chr, start = vntrseek_data[vntr_id].split()
        vntr = ReferenceVNTR(vntr_id, pattern, int(start)-1, vntrseek_repeat, chr)
        vntr.init_from_xml(repeat_segments)
        vntrs.append(vntr)
    return vntrs


process_vntrseek_data()
