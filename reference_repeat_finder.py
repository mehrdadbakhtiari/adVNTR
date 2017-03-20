from Bio import Seq, SeqIO
from pomegranate import DiscreteDistribution, DirichletDistribution, State
from pomegranate import HiddenMarkovModel as Model


def path_to_alignment(x, y, path):
    for i, (index, state) in enumerate(path[1:-1]):
        name = state.name

        if name.startswith('D'):
            y = y[:i] + '-' + y[i:]
        elif name.startswith('I'):
            x = x[:i] + '-' + x[i:]

    return x, y


def build_hmm(patterns, copies=1):
    pattern = patterns[0]
    model = Model(name="HMM Model")
    insert_distribution = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})

    last_end = None
    start_random_matches = State(insert_distribution, name='start_random_matches')
    end_random_matches = State(insert_distribution, name='end_random_matches')
    model.add_states([start_random_matches, end_random_matches])
    for repeat in range(copies):
        insert_states = []
        match_states = []
        delete_states = []
        for i in range(len(pattern) + 1):
            insert_states.append(State(insert_distribution, name='I%s_%s' % (i, repeat)))

        for i in range(len(pattern)):
            distribution_map = {'A': 0.01, 'C': 0.01, 'G': 0.01, 'T': 0.01}
            distribution_map[pattern[i]] = 0.97
            match_states.append(State(DiscreteDistribution(distribution_map), name='M%s_%s' % (str(i + 1), repeat)))

        for i in range(len(pattern)):
            delete_states.append(State(None, name='D%s_%s' % (str(i + 1), repeat)))

        unit_start = State(None, name='unit_start_%s' % repeat)
        unit_end = State(None, name='unit_end_%s' % repeat)
        model.add_states(insert_states + match_states + delete_states + [unit_start, unit_end])
        last = len(delete_states)-1

        if repeat > 0:
            model.add_transition(last_end, unit_start, 0.5)
        else:
            model.add_transition(model.start, unit_start, 0.5)
            model.add_transition(model.start, start_random_matches, 0.5)
            model.add_transition(start_random_matches, unit_start, 0.5)
            model.add_transition(start_random_matches, start_random_matches, 0.5)

        model.add_transition(unit_end, end_random_matches, 0.5)
        if repeat == copies - 1:
            model.add_transition(unit_end, model.end, 0.5)
            model.add_transition(end_random_matches, end_random_matches, 0.5)
            model.add_transition(end_random_matches, model.end, 0.5)

        model.add_transition(unit_start, match_states[0], 0.98)
        model.add_transition(unit_start, delete_states[0], 0.01)
        model.add_transition(unit_start, insert_states[0], 0.01)

        model.add_transition(insert_states[0], insert_states[0], 0.01)
        model.add_transition(insert_states[0], delete_states[0], 0.01)
        model.add_transition(insert_states[0], match_states[0], 0.98)

        model.add_transition(delete_states[last], unit_end, 0.99)
        model.add_transition(delete_states[last], insert_states[last+1], 0.01)

        model.add_transition(match_states[last], unit_end, 0.99)
        model.add_transition(match_states[last], insert_states[last+1], 0.01)

        model.add_transition(insert_states[last+1], insert_states[last+1], 0.01)
        model.add_transition(insert_states[last+1], unit_end, 0.99)

        for i in range(0, len(pattern)):
            model.add_transition(match_states[i], insert_states[i+1], 0.01)
            model.add_transition(delete_states[i], insert_states[i+1], 0.01)
            model.add_transition(insert_states[i+1], insert_states[i+1], 0.01)
            if i < len(pattern) - 1:
                model.add_transition(insert_states[i+1], match_states[i+1], 0.98)
                model.add_transition(insert_states[i+1], delete_states[i+1], 0.01)

                model.add_transition(match_states[i], match_states[i+1], 0.98)
                model.add_transition(match_states[i], delete_states[i+1], 0.01)

                model.add_transition(delete_states[i], delete_states[i+1], 0.01)
                model.add_transition(delete_states[i], match_states[i+1], 0.98)

        last_end = unit_end

    model.bake()
    if len(patterns) > 1:
        # model.fit(patterns, algorithm='baum-welch', transition_pseudocount=1, use_pseudocount=True)
        fit_patterns = [pattern * copies for pattern in patterns]
        model.fit(fit_patterns, algorithm='viterbi', transition_pseudocount=1, use_pseudocount=True)

    return model


def is_matching_state(state_name):
    if state_name.startswith('M') or state_name.startswith('I') or state_name.startswith('start_random_matches')or state_name.startswith('end_random_matches'):
        return True
    return False


def get_scores_and_segments(corresponding_region_in_ref, patterns, copies=40):
    model = build_hmm(patterns, copies=copies)
    single_model = build_hmm(patterns, copies=1)
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
    hmm = build_hmm(repeat_segments, copies=1)
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
        j = i + 1
        while j < len(patterns) and len(patterns[i]) * repeat_count + start_points[i] > start_points[j]:
            skipped_vntrs.append(j)
            j += 1
        if i in skipped_vntrs:
            repeat_count = 0
        with open('pattern_repeat_counts.txt', 'a') as out:
            out.write('%s\n' % repeat_count)
        with open('visited_states.txt', 'a') as out:
            out.write('%s\n' % ' '.join(states))

find_true_repeat_counts()
