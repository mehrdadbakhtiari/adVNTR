from Bio import SeqIO
from hmm_utils import is_matching_state, build_reference_repeat_finder_hmm


def get_scores_and_segments(corresponding_region_in_ref, patterns, copies=40):
    model = build_reference_repeat_finder_hmm(patterns, copies=copies)
    single_model = build_reference_repeat_finder_hmm(patterns, copies=1)
    logp, path = model.viterbi(corresponding_region_in_ref)
    visited_states = [state.name for idx, state in path[1:-1]]

    lengths = []
    prev_start = None
    for i in range(len(visited_states)):
        if visited_states[i].startswith('unit_end') and prev_start is not None:
            current_len = 0
            for j in range(prev_start, i):
                if is_matching_state(visited_states[j]):
                    current_len += 1
            lengths.append(current_len)
        if visited_states[i].startswith('unit_start'):
            prev_start = i

    repeat_segments = []
    scores = []
    added = 0
    for l in lengths:
        repeat_segments.append(corresponding_region_in_ref[added:added+l])
        score, temp_path = single_model.viterbi(corresponding_region_in_ref[added:added+l])
        scores.append(score)
        added += l
    return scores, repeat_segments, visited_states


def find_number_of_tandem_repeats_in_reference(pattern, pattern_start, copies, ref_file_name='chr15.fa'):
    fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
    ref_sequence = ''
    for fasta in fasta_sequences:
        name, ref_sequence = fasta.id, str(fasta.seq)
    corresponding_region_in_ref = ref_sequence[pattern_start:pattern_start + len(pattern) * copies].upper()

    original_scores, repeat_segments, states = get_scores_and_segments(corresponding_region_in_ref, [pattern], copies)
    hmm = build_reference_repeat_finder_hmm(repeat_segments, copies=1)
    scores = []
    for seg in repeat_segments:
        score, temp_path = hmm.viterbi(seg)
        scores.append(score)

    return len(repeat_segments), repeat_segments, states


def plot_segment_scores(scores, pattern):
    import matplotlib.pyplot as plt
    X = [i for i, score in enumerate(scores)]
    Y = [score for i, score in enumerate(scores)]
    plt.plot(X, Y, color='blue', label=pattern)
    plt.xlabel('pattern number')
    plt.ylabel('logp')
    plt.legend(loc=0)
    plt.savefig('hmm_p_%s_fitted_by_first20.png' % pattern)  # save the figure to file
    plt.close()


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
        copies = vntrseek_repeats[i] + 2
        if i < len(patterns) - 1 and len(patterns[i]) * copies + start_points[i] > start_points[i+1]:
            copies += vntrseek_repeats[i+1]
        copies += 3
        repeat_count, repeat_segments, states = find_number_of_tandem_repeats_in_reference(patterns[i], start_points[i], copies)
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
