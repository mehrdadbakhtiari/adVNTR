from pomegranate import DiscreteDistribution, State
from pomegranate import HiddenMarkovModel as Model
import numpy as np

from profile_hmm import build_profile_hmm_for_repeats
import settings


def path_to_alignment(x, y, path):
    for i, (index, state) in enumerate(path[1:-1]):
        name = state.name

        if name.startswith('D'):
            y = y[:i] + '-' + y[i:]
        elif name.startswith('I'):
            x = x[:i] + '-' + x[i:]

    return x, y


def get_emitted_basepair_from_visited_states(state, visited_states, sequence):
    base_pair_idx = 0
    for visited_state in visited_states:
        if visited_state == state:
            return sequence[base_pair_idx]
        if is_matching_state(visited_state):
            base_pair_idx += 1
    return None


def is_matching_state(state_name):
    if state_name.startswith('M') or state_name.startswith('I') or state_name.startswith('start_random_matches') \
            or state_name.startswith('end_random_matches'):
        return True
    return False


def get_repeating_pattern_lengths(visited_states):
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
    return lengths


def get_repeat_segments_from_visited_states_and_region(visited_states, region):
    lengths = get_repeating_pattern_lengths(visited_states)

    repeat_segments = []
    added = 0
    for l in lengths:
        repeat_segments.append(region[added:added + l])
        added += l
    return repeat_segments


def get_number_of_repeats_in_vpath(vpath):
    starts = 0
    ends = 0
    visited_states = [state.name for idx, state in vpath[1:-1]]
    for i in range(len(visited_states)):
        if visited_states[i].startswith('unit_start'):
            starts += 1
        if visited_states[i].startswith('unit_end'):
            ends += 1
    return max(starts, ends)


def get_number_of_repeat_bp_matches_in_vpath(vpath):
    visited_states = [state.name for idx, state in vpath[1:-1]]
    result = 0
    for i in range(len(visited_states)):
        if is_matching_state(visited_states[i]) and not visited_states[i].endswith('fix'):
            result += 1
    return result


def get_left_flanking_region_size_in_vpath(vpath):
    visited_states = [state.name for idx, state in vpath[1:-1]]
    result = 0
    for i in range(len(visited_states)):
        if is_matching_state(visited_states[i]) and visited_states[i].endswith('suffix'):
            result += 1
    return result


def get_right_flanking_region_size_in_vpath(vpath):
    visited_states = [state.name for idx, state in vpath[1:-1]]
    result = 0
    for i in range(len(visited_states)):
        if is_matching_state(visited_states[i]) and visited_states[i].endswith('prefix'):
            result += 1
    return result


def get_prefix_matcher_hmm(pattern):
    model = Model(name="Prefix Matcher HMM Model")
    insert_distribution = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
    insert_states = []
    match_states = []
    delete_states = []
    hmm_name = 'prefix'
    for i in range(len(pattern) + 1):
        insert_states.append(State(insert_distribution, name='I%s_%s' % (i, hmm_name)))

    for i in range(len(pattern)):
        distribution_map = {'A': 0.01, 'C': 0.01, 'G': 0.01, 'T': 0.01}
        distribution_map[pattern[i]] = 0.97
        match_states.append(State(DiscreteDistribution(distribution_map), name='M%s_%s' % (str(i + 1), hmm_name)))

    for i in range(len(pattern)):
        delete_states.append(State(None, name='D%s_%s' % (str(i + 1), hmm_name)))

    unit_start = State(None, name='prefix_start_%s' % hmm_name)
    unit_end = State(None, name='prefix_end_%s' % hmm_name)
    model.add_states(insert_states + match_states + delete_states + [unit_start, unit_end])
    last = len(delete_states)-1

    model.add_transition(model.start, unit_start, 1)

    model.add_transition(unit_end, model.end, 1)

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

            model.add_transition(match_states[i], match_states[i+1], 0.97)
            model.add_transition(match_states[i], delete_states[i+1], 0.01)
            model.add_transition(match_states[i], unit_end, 0.01)

            model.add_transition(delete_states[i], delete_states[i+1], 0.01)
            model.add_transition(delete_states[i], match_states[i+1], 0.98)

    model.bake(merge=None)

    return model


def get_suffix_matcher_hmm(pattern):
    model = Model(name="Suffix Matcher HMM Model")
    insert_distribution = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
    insert_states = []
    match_states = []
    delete_states = []
    hmm_name = 'suffix'
    for i in range(len(pattern) + 1):
        insert_states.append(State(insert_distribution, name='I%s_%s' % (i, hmm_name)))

    for i in range(len(pattern)):
        distribution_map = {'A': 0.01, 'C': 0.01, 'G': 0.01, 'T': 0.01}
        distribution_map[pattern[i]] = 0.97
        match_states.append(State(DiscreteDistribution(distribution_map), name='M%s_%s' % (str(i + 1), hmm_name)))

    for i in range(len(pattern)):
        delete_states.append(State(None, name='D%s_%s' % (str(i + 1), hmm_name)))

    unit_start = State(None, name='suffix_start_%s' % hmm_name)
    unit_end = State(None, name='suffix_end_%s' % hmm_name)
    model.add_states(insert_states + match_states + delete_states + [unit_start, unit_end])
    last = len(delete_states)-1

    model.add_transition(model.start, unit_start, 1)

    model.add_transition(unit_end, model.end, 1)

    model.add_transition(unit_start, delete_states[0], 0.01)
    model.add_transition(unit_start, insert_states[0], 0.01)
    for i in range(len(pattern)):
        model.add_transition(unit_start, match_states[i], 0.98 / len(pattern))

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

    model.bake(merge=None)

    return model


def get_trained_model_for_one_mixed_repeat(patterns):
    pattern = patterns[0]
    model = Model(name="Transition Finder HMM")
    insert_distribution = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})
    insert_states = []
    match_states = []
    delete_states = []
    for i in range(len(pattern) + 1):
        insert_states.append(State(insert_distribution, name='I%s' % i))

    for i in range(len(pattern)):
        distribution_map = {'A': 0.01, 'C': 0.01, 'G': 0.01, 'T': 0.01}
        distribution_map[pattern[i]] = 0.97
        match_states.append(State(DiscreteDistribution(distribution_map), name='M%s' % str(i)))

    for i in range(len(pattern)):
        delete_states.append(State(None, name='D%s' % str(i)))

    unit_start = State(None, name='unit_start')
    unit_end = State(None, name='unit_end')
    model.add_states(insert_states + match_states + delete_states + [unit_start, unit_end])
    last = len(delete_states)-1

    model.add_transition(model.start, unit_start, 1)
    model.add_transition(unit_end, model.end, 1)

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

    model.bake(merge=None)
    if len(patterns) > 1:
        model.fit(patterns, algorithm='viterbi', transition_pseudocount=1, use_pseudocount=True, verbose=False,
                  n_jobs=settings.CORES)
    return model


def get_constant_number_of_repeats_matcher_hmm(patterns, copies):
    pattern = patterns[0]
    model = Model(name="Repeating Pattern Matcher HMM Model")

    transitions, emissions = build_profile_hmm_for_repeats(patterns, settings.MAX_ERROR_RATE)

    last_end = None
    for repeat in range(copies):
        insert_states = []
        match_states = []
        delete_states = []
        for i in range(len(pattern) + 1):
            insert_distribution = DiscreteDistribution(emissions['I%s' % i])
            insert_states.append(State(insert_distribution, name='I%s_%s' % (i, repeat)))

        for i in range(1, len(pattern) + 1):
            match_distribution = DiscreteDistribution(emissions['M%s' % i])
            match_states.append(State(match_distribution, name='M%s_%s' % (str(i), repeat)))

        for i in range(1, len(pattern) + 1):
            delete_states.append(State(None, name='D%s_%s' % (str(i), repeat)))

        unit_start = State(None, name='unit_start_%s' % repeat)
        unit_end = State(None, name='unit_end_%s' % repeat)
        model.add_states(insert_states + match_states + delete_states + [unit_start, unit_end])
        n = len(delete_states)-1

        if repeat > 0:
            model.add_transition(last_end, unit_start, 1)
        else:
            model.add_transition(model.start, unit_start, 1)

        if repeat == copies - 1:
            model.add_transition(unit_end, model.end, 1)

        model.add_transition(unit_start, match_states[0], transitions['unit_start']['M1'])
        model.add_transition(unit_start, delete_states[0], transitions['unit_start']['D1'])
        model.add_transition(unit_start, insert_states[0], transitions['unit_start']['I0'])

        model.add_transition(insert_states[0], insert_states[0], transitions['I0']['I0'])
        model.add_transition(insert_states[0], delete_states[0], transitions['I0']['D1'])
        model.add_transition(insert_states[0], match_states[0], transitions['I0']['M1'])

        model.add_transition(delete_states[n], unit_end, transitions['D%s' % (n+1)]['unit_end'])
        model.add_transition(delete_states[n], insert_states[n+1], transitions['D%s' % (n+1)]['I%s' % (n+1)])

        model.add_transition(match_states[n], unit_end, transitions['M%s' % (n+1)]['unit_end'])
        model.add_transition(match_states[n], insert_states[n+1], transitions['M%s' % (n+1)]['I%s' % (n+1)])

        model.add_transition(insert_states[n+1], insert_states[n+1], transitions['I%s' % (n+1)]['I%s' % (n+1)])
        model.add_transition(insert_states[n+1], unit_end, transitions['I%s' % (n+1)]['unit_end'])

        for i in range(1, len(pattern)+1):
            model.add_transition(match_states[i-1], insert_states[i], transitions['M%s' % i]['I%s' % i])
            model.add_transition(delete_states[i-1], insert_states[i], transitions['D%s' % i]['I%s' % i])
            model.add_transition(insert_states[i], insert_states[i], transitions['I%s' % i]['I%s' % i])
            if i < len(pattern):
                model.add_transition(insert_states[i], match_states[i], transitions['I%s' % i]['M%s' % (i+1)])
                model.add_transition(insert_states[i], delete_states[i], transitions['I%s' % i]['D%s' % (i+1)])

                model.add_transition(match_states[i-1], match_states[i], transitions['M%s' % i]['M%s' % (i+1)])
                model.add_transition(match_states[i-1], delete_states[i], transitions['M%s' % i]['D%s' % (i+1)])

                model.add_transition(delete_states[i-1], match_states[i], transitions['D%s' % i]['M%s' % (i+1)])
                model.add_transition(delete_states[i-1], delete_states[i], transitions['D%s' % i]['D%s' % (i+1)])

        last_end = unit_end

    model.bake(merge=None)
    return model


def get_variable_number_of_repeats_matcher_hmm(patterns, copies=1):
    model = get_constant_number_of_repeats_matcher_hmm(patterns, copies)

    start_repeats_matches = State(None, name='start_repeating_pattern_match')
    end_repeats_matches = State(None, name='end_repeating_pattern_match')
    mat = model.dense_transition_matrix()
    states = model.states
    states.append(start_repeats_matches)
    states.append(end_repeats_matches)
    states_count = len(mat)
    start_repeats_ind = states_count
    end_repeats_ind = states_count + 1
    mat = np.c_[mat, np.zeros(states_count), np.zeros(states_count)]
    mat = np.r_[mat, [np.zeros(states_count + 2)]]
    mat = np.r_[mat, [np.zeros(states_count + 2)]]

    unit_ends = []
    for i, state in enumerate(model.states):
        if state.name.startswith('unit_end'):
            unit_ends.append(i)

    for i in range(len(mat[model.start_index])):
        if mat[model.start_index][i] != 0:
            first_unit_start = i
    mat[model.start_index][first_unit_start] = 0.0
    mat[model.start_index][start_repeats_ind] = 1
    mat[start_repeats_ind][first_unit_start] = 1

    for unit_end in unit_ends:
        for j in range(len(mat[unit_end])):
            if mat[unit_end][j] != 0:
                next_state = j
        mat[unit_end][next_state] = 0.5
        mat[unit_end][end_repeats_ind] = 0.5

    mat[end_repeats_ind][model.end_index] = 1

    starts = np.zeros(states_count + 2)
    starts[model.start_index] = 1.0
    ends = np.zeros(states_count + 2)
    ends[model.end_index] = 1.0
    state_names = [state.name for state in states]
    distributions = [state.distribution for state in states]
    new_model = Model.from_matrix(mat, distributions, starts, ends, name='Repeat Matcher HMM Model', state_names=state_names, merge=None)
    new_model.bake(merge=None)
    return new_model


def get_read_matcher_model(left_flanking_region, right_flanking_region, patterns, copies=1):
    model = get_suffix_matcher_hmm(left_flanking_region)
    right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
    repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
    model.concatenate(repeats_matcher)
    model.concatenate(right_flanking_matcher)
    model.bake(merge=None)

    mat = model.dense_transition_matrix()

    first_repeat_matches = []
    repeat_match_states = []
    suffix_start = None
    for i, state in enumerate(model.states):
        if state.name[0] == 'M' and state.name.split('_')[-1] == '0':
            first_repeat_matches.append(i)
        if state.name[0] == 'M' and state.name.split('_')[-1] not in ['prefix', 'suffix']:
            repeat_match_states.append(i)
        if state.name == 'suffix_start_suffix':
            suffix_start = i

    mat[model.start_index][suffix_start] = 0.3
    for first_repeat_match in first_repeat_matches:
        mat[model.start_index][first_repeat_match] = 0.7 / len(first_repeat_matches)

    for match_state in repeat_match_states:
        to_end = 0.7 / len(repeat_match_states)
        total = 1 + to_end
        for next in range(len(mat[match_state])):
            if mat[match_state][next] != 0:
                mat[match_state][next] /= total
        mat[match_state][model.end_index] = to_end

    starts = np.zeros(len(model.states))
    starts[model.start_index] = 1.0
    ends = np.zeros(len(model.states))
    ends[model.end_index] = 1.0
    state_names = [state.name for state in model.states]
    distributions = [state.distribution for state in model.states]
    new_model = Model.from_matrix(mat, distributions, starts, ends, name='Read Matcher', state_names=state_names, merge=None)
    return new_model


def build_reference_repeat_finder_hmm(patterns, copies=1):
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
