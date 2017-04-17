from pomegranate import DiscreteDistribution, State
from pomegranate import HiddenMarkovModel as Model
import numpy as np


def is_matching_state(state_name):
    if state_name.startswith('M') or state_name.startswith('I') or state_name.startswith('start_random_matches')or state_name.startswith('end_random_matches'):
        return True
    return False


def get_number_of_repeat_bp_matches_in_vpath(vpath):
    visited_states = [state.name for idx, state in vpath[1:-1]]
    result = 0
    counting = False
    for i in range(len(visited_states)):
        if visited_states[i].startswith('start_repeating_pattern_match'):
            counting = True
        if visited_states[i] == 'end_repeating_pattern_match':
            counting = False
        if is_matching_state(visited_states[i]) and counting:
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


def get_constant_number_of_repeats_matcher_hmm(patterns, copies):
    pattern = patterns[0]
    model = Model(name="Repeating Pattern Matcher HMM Model")
    insert_distribution = DiscreteDistribution({'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25})

    last_end = None
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
            model.add_transition(last_end, unit_start, 1)
        else:
            model.add_transition(model.start, unit_start, 1)

        if repeat == copies - 1:
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

        last_end = unit_end

    model.bake(merge=None)
    if len(patterns) > 1:
        # model.fit(patterns, algorithm='baum-welch', transition_pseudocount=1, use_pseudocount=True)
        fit_patterns = [pattern * copies for pattern in patterns]
        model.fit(fit_patterns, algorithm='viterbi', transition_pseudocount=1, use_pseudocount=True)
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
    first_repeat_matches = []
    matches_states = []
    for i, state in enumerate(model.states):
        if state.name.startswith('unit_end'):
            unit_ends.append(i)
        if state.name[0] == 'M' and state.name.split('_')[-1] == '0':
            first_repeat_matches.append(i)
        if state.name[0] == 'M':
            matches_states.append(i)

    for i in range(len(mat[model.start_index])):
        if mat[model.start_index][i] != 0:
            first_unit_start = i
    mat[model.start_index][first_unit_start] = 0.0
    mat[model.start_index][start_repeats_ind] = 1
    mat[start_repeats_ind][first_unit_start] = 0.3
    for first_repeat_match in first_repeat_matches:
        mat[start_repeats_ind][first_repeat_match] = 0.7 / len(first_repeat_matches)

    for match_state in matches_states:
        to_end = 0.7 / len(matches_states)
        total = 1 + to_end
        for next in range(len(mat[match_state])):
            if mat[match_state][next] != 0:
                mat[match_state][next] /=  total
        mat[match_state][end_repeats_ind] = to_end

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
    return new_model
