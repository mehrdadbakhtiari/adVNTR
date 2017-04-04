from typing import re

from Bio import Seq, SeqIO
from pomegranate import DiscreteDistribution, State
from pomegranate import HiddenMarkovModel as Model
import numpy as np
from repeat_finder import get_blast_matched_ids


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

    for i in range(len(mat[model.start_index])):
        if mat[model.start_index][i] != 0:
            first_unit_start = i
    mat[model.start_index][first_unit_start] = 0.0
    mat[model.start_index][start_repeats_ind] = 1
    mat[start_repeats_ind][first_unit_start] = 1

    unit_ends = []
    for i, state in enumerate(model.states):
        if state.name.startswith('unit_end'):
            unit_ends.append(i)

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


def get_VNTR_matcher_hmm(patterns, copies, left_flanking_region, right_flanking_region):
    left_flanking_matcher = get_suffix_matcher_hmm(left_flanking_region)
    right_flanking_matcher = get_prefix_matcher_hmm(right_flanking_region)
    repeats_matcher = get_variable_number_of_repeats_matcher_hmm(patterns, copies)
    left_flanking_matcher.concatenate(repeats_matcher)
    left_flanking_matcher.concatenate(right_flanking_matcher)
    left_flanking_matcher.bake(merge=None)
    return left_flanking_matcher


def is_matching_state(state_name):
    if state_name.startswith('M') or state_name.startswith('I') or state_name.startswith('start_random_matches')or state_name.startswith('end_random_matches'):
        return True
    return False


def extract_repeat_segments_from_visited_states(pattern, pattern_start, copies, visited_states, ref_file_name='chr15.fa'):
    fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
    ref_sequence = ''
    for fasta in fasta_sequences:
        name, ref_sequence = fasta.id, str(fasta.seq)
    corresponding_region_in_ref = ref_sequence[pattern_start:pattern_start + (len(pattern) + 5) * copies].upper()

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
    added = 0
    for l in lengths:
        repeat_segments.append(corresponding_region_in_ref[added:added+l])
        added += l
    return repeat_segments


def get_flanking_regions(start_point, end_point, flanking_region_size, ref_file_name='chr15.fa'):
    fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
    ref_sequence = ''
    for fasta in fasta_sequences:
        name, ref_sequence = fasta.id, str(fasta.seq)
    left_flanking = ref_sequence[start_point - flanking_region_size:start_point].upper()
    right_flanking = ref_sequence[end_point:end_point + flanking_region_size].upper()
    return left_flanking, right_flanking


def get_number_of_matches_in_a_read(vpath):
    visited_states = [state.name for idx, state in vpath[1:-1]]
    result = 0
    for i in range(len(visited_states)):
        if visited_states[i].startswith('unit_end'):
            result += 1
    if result < 2:
        result = 0
    return result


def find_repeat_count(pattern_num, pattern, start_point, repeat_count, visited_states, read_files):
    repeat_segments = extract_repeat_segments_from_visited_states(pattern, start_point, repeat_count, visited_states)
    copies = int(round(150.0 / len(pattern)))
    flanking_region_size = 150 - len(pattern)
    end_point = start_point + sum([len(e) for e in repeat_segments])
    left_flanking_region, right_flanking_region = get_flanking_regions(start_point, end_point, flanking_region_size)
    hmm = get_VNTR_matcher_hmm(repeat_segments * 100, copies, left_flanking_region, right_flanking_region)
    total_occurrences = 0

    word_size = '6'
    if len(pattern) > 30:
        word_size = '11'
    blast_ids = get_blast_matched_ids(pattern, 'original_reads/original_reads', max_seq='10000', evalue=10,
                                      word_size=word_size, search_id=str(pattern_num))
    if len(blast_ids) == 10000:
        with open('errors.txt', 'a') as out:
            out.write('%s\n' % pattern_num)

    number_of_reads = 0
    read_length = 0
    total_length = 100 * 1000 * 1000
    for read_file in read_files:
        print('opening read file')
        reads = SeqIO.parse(read_file, 'fasta')
        for read_segment in reads:
            if number_of_reads == 0:
                read_length = len(str(read_segment.seq))
            number_of_reads += 1
            if read_segment.id not in blast_ids:
                continue
            occurrence = 0
            logp, vpath = hmm.viterbi(str(read_segment.seq))
            if vpath and logp > -170:
                occurrence = get_number_of_matches_in_a_read(vpath)
            logp, vpath = hmm.viterbi(str(read_segment.seq.reverse_complement()))
            if vpath and logp > -170:
                occurrence = max(occurrence, get_number_of_matches_in_a_read(vpath))
            total_occurrences += occurrence

            number_of_reads += 1
    avg_coverage = float(number_of_reads * read_length) / total_length
    print('occurrence: %s, avg_coverage: %s' % (total_occurrences, avg_coverage))
    return total_occurrences / avg_coverage


with open('patterns.txt') as input:
    patterns = input.readlines()
    patterns = [pattern.strip() for pattern in patterns]
with open('start_points.txt') as input:
    lines = input.readlines()
    start_points = [int(num.strip())-1 for num in lines]
with open('pattern_repeat_counts.txt') as input:
    lines = input.readlines()
    repeat_counts = [int(num.strip()) for num in lines]
with open('visited_states.txt') as input:
    lines = input.readlines()
    visited_states_list = [states.strip().split() for states in lines]

read_files = ['original_reads/paired_dat1.fasta', 'original_reads/paired_dat2.fasta']
for i in range(len(patterns)):
    print(i)
    if repeat_counts[i] == 0:
        continue
    cn = find_repeat_count(i, patterns[i], start_points[i], repeat_counts[i], visited_states_list[i], read_files)
    with open('hmm_repeat_count.txt', 'a') as output:
        output.write('%s %s\n' % (i, cn / repeat_counts[i]))
