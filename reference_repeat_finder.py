from Bio import SeqIO
from hmm_utils import build_reference_repeat_finder_hmm, get_repeat_segments_from_visited_states_and_region


class ReferenceVNTR:
    def __init__(self, id, pattern, start_point, estimated_repeats=None, chromosome='chr15'):
        self.id = id
        self.pattern = pattern
        self.start_point = start_point
        self.estimated_repeats = estimated_repeats
        self.chromosome = chromosome.lower()

    def get_scores_and_segments(self, region_in_ref):
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

    def find_number_of_tandem_repeats_in_reference(self):
        corresponding_region_in_ref = self.get_corresponding_region_in_ref()
        repeat_segments = self.get_scores_and_segments(corresponding_region_in_ref)
        hmm = build_reference_repeat_finder_hmm(repeat_segments, copies=1)
        scores = []
        for seg in repeat_segments:
            score, temp_path = hmm.viterbi(seg)
            scores.append(score)

        return len(repeat_segments), repeat_segments


def find_true_repeat_counts():

    vntrs = []
    with open ('repeats_length_patterns_chromosomes_starts.txt') as input:
        input_lines = [line.strip() for line in input.readlines() if line.strip() != '']
        for vntr_id, line in enumerate(input_lines):
            vntrseek_repeat, _, pattern, chr, start = line.split()
            estimated_repeats = int(vntrseek_repeat) + 2
            vntrs.append(ReferenceVNTR(vntr_id, pattern, int(start), estimated_repeats, chr))

    skipped_vntrs = []
    for i in range(len(vntrs)):
        print(i)
        estimated_end_point = len(vntrs[i].pattern) * vntrs[i].estimated_repeats + vntrs[i].start_point
        if i < len(vntrs) - 1 and vntrs[i].chromosome == vntrs[i+1].chromosome and estimated_end_point > vntrs[i+1].start_point:
            vntrs[i].estimated_repeats += vntrs[i+1].estimated_repeats
        vntrs[i].estimated_repeats += 3
        repeat_count, repeat_segments = vntrs[i].find_number_of_tandem_repeats_in_reference()
        if i in skipped_vntrs:
            repeat_count = 0
        else:
            j = i + 1
            end_point = len(vntrs[i].pattern) * repeat_count + vntrs[i].start_point
            while j < len(vntrs) and vntrs[i].chromosome == vntrs[j].chromosome and end_point > vntrs[j].start_point:
                skipped_vntrs.append(j)
                j += 1
        comma_separated_segments = ','.join(repeat_segments)
        with open('repeats_and_segments.txt', 'a') as out:
            out.write('%s %s' % (repeat_count, comma_separated_segments))

find_true_repeat_counts()
