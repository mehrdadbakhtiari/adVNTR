from Bio import SeqIO
from hmm_utils import build_reference_repeat_finder_hmm, get_repeat_segments_from_visited_states_and_region


class ReferenceVNTR:
    def __init__(self, pattern, start_point, estimated_repeats, chromosome='chr15'):
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

        return repeat_segments, visited_states

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
        repeat_segments, states = self.get_scores_and_segments(corresponding_region_in_ref)
        hmm = build_reference_repeat_finder_hmm(repeat_segments, copies=1)
        scores = []
        for seg in repeat_segments:
            score, temp_path = hmm.viterbi(seg)
            scores.append(score)

        return len(repeat_segments), repeat_segments, states


def find_true_repeat_counts():
    with open('patterns.txt') as input:
        patterns = input.readlines()
        patterns = [pattern.strip() for pattern in patterns]
    with open('start_points.txt') as input:
        lines = input.readlines()
        start_points = [int(num.strip()) - 1 for num in lines]
    with open('vntrseek_repeats.txt') as input:
        lines = input.readlines()
        vntrseek_repeats = [int(float(num.strip())) - 1 for num in lines]

    repeats = []
    skipped_vntrs = []
    for i in range(len(patterns)):
        print(i)
        estimated_repeats = vntrseek_repeats[i] + 2
        if i < len(patterns) - 1 and len(patterns[i]) * estimated_repeats + start_points[i] > start_points[i+1]:
            estimated_repeats += vntrseek_repeats[i+1]
        estimated_repeats += 3
        reference_vntr = ReferenceVNTR(patterns[i], start_points[i], estimated_repeats)
        repeat_count, repeat_segments, states = reference_vntr.find_number_of_tandem_repeats_in_reference()
        repeats.append(repeat_count)
        if i in skipped_vntrs:
            repeat_count = 0
        else:
            j = i + 1
            while j < len(patterns) and len(patterns[i]) * repeat_count + start_points[i] > start_points[j]:
                skipped_vntrs.append(j)
                j += 1
        with open('pattern_repeat_counts.txt', 'a') as out:
            out.write('%s\n' % repeat_count)
        with open('visited_states.txt', 'a') as out:
            out.write('%s\n' % ' '.join(states))

find_true_repeat_counts()
