try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

from advntr.profiler import time_usage


@time_usage
def build_profile_hmm_pseudocounts_for_alignment(thresh, pseu, alphabet, alignment):
    thresh = thresh * len(alignment)
    states = []
    insert_index = []
    for i in range(len(alignment[0])):
        total = 0
        for j in alignment:
            if j[i] == '-':
                total += 1
        if total >= thresh:
            insert_index.append(i)

    emission = dict({})
    emission['unit_start'] = dict((x, 0) for x in alphabet)
    emission['unit_end'] = dict((x, 0) for x in alphabet)
    emission['I0'] = dict((x, 0) for x in alphabet)
    for i in range(1, len(alignment[0]) - len(insert_index) + 1):
        emission['I' + str(i)] = dict((x, 0) for x in alphabet)
        emission['M' + str(i)] = dict((x, 0) for x in alphabet)
        emission['D' + str(i)] = dict((x, 0) for x in alphabet)

    for i in range(len(alignment)):
        per_state = []
        current_index = 1
        for j in range(len(alignment[i])):
            if j not in insert_index:
                if alignment[i][j] == '-':
                    per_state.append('D' + str(current_index))
                else:
                    per_state.append('M' + str(current_index))
                    emission['M' + str(current_index)][alignment[i][j]] += 1
                current_index += 1
            else:
                if alignment[i][j] != '-':
                    per_state.append('I' + str(current_index - 1))
                    emission['I' + str(current_index - 1)][alignment[i][j]] += 1
        states.append(per_state)

    for key in emission:
        if key != 'unit_start' and key != 'unit_end' and not key.startswith('D'):
            total = 0
            for sub_key in emission[key]:
                total += emission[key][sub_key]
            if total > 0:
                sub_total = 0
                for sub_key in emission[key]:
                    emission[key][sub_key] = (1.0 * emission[key][sub_key]) / total + pseu
                    # emission[key][sub_key] = (emission[key][sub_key]+pseu)/(1+pseu*len(emission[key]))
                    # emission[key][sub_key] = (emission[key][sub_key]+pseu)/(1+pseu*len(alphabet))
                    sub_total += 1.0 * emission[key][sub_key]
                for sub_key in emission[key]:
                    emission[key][sub_key] = emission[key][sub_key] / sub_total
            else:
                for sub_key in emission[key]:
                    emission[key][sub_key] = 1.0 / len(alphabet)
    transition = dict({})
    transition['unit_start'] = {}
    transition['unit_start']['I0'] = 0
    transition['unit_start']['D1'] = 0
    transition['unit_start']['M1'] = 0
    for i in states:
        transition['unit_start'][i[0]] += 1
    transition['I0'] = {}
    transition['I0']['I0'] = 0
    transition['I0']['D1'] = 0
    transition['I0']['M1'] = 0
    # for i in states:
    #     transition['I0'][i[0]] += 1

    for i in range(len(states)):
        for j in range(len(states[i]) - 1):
            if states[i][j] not in transition:
                transition[states[i][j]] = {}
            if states[i][j + 1] not in transition[states[i][j]]:
                transition[states[i][j]][states[i][j + 1]] = 0
            transition[states[i][j]][states[i][j + 1]] += 1
        if not states[i][len(states[i]) - 1] in transition:
            transition[states[i][len(states[i]) - 1]] = {}
        if 'unit_end' not in transition[states[i][len(states[i]) - 1]]:
            transition[states[i][len(states[i]) - 1]]['unit_end'] = 0
        transition[states[i][len(states[i]) - 1]]['unit_end'] += 1

    if 'I' + str(len(alignment[0]) - len(insert_index)) not in transition:
        transition['I' + str(len(alignment[0]) - len(insert_index))] = {}
        if 'unit_end' not in transition['I' + str(len(alignment[0]) - len(insert_index))]:
            transition['I' + str(len(alignment[0]) - len(insert_index))]['unit_end'] = 0

    if 'D' + str(len(alignment[0]) - len(insert_index)) not in transition:
        transition['D' + str(len(alignment[0]) - len(insert_index))] = {}
        if 'unit_end' not in transition['D' + str(len(alignment[0]) - len(insert_index))]:
            transition['D' + str(len(alignment[0]) - len(insert_index))]['unit_end'] = 0

    if 'M' + str(len(alignment[0]) - len(insert_index)) not in transition:
        transition['M' + str(len(alignment[0]) - len(insert_index))] = {}
        if 'unit_end' not in transition['M' + str(len(alignment[0]) - len(insert_index))]:
            transition['M' + str(len(alignment[0]) - len(insert_index))]['unit_end'] = 0

    for i in range(1, len(alignment[0]) - len(insert_index) + 1):
        if not 'I' + str(i) in transition:
            transition['I' + str(i)] = {}
        if not 'M' + str(i) in transition:
            transition['M' + str(i)] = {}
        if not 'D' + str(i) in transition:
            transition['D' + str(i)] = {}

    for key in transition:
        if key != 'unit_end':
            total = 0
            for sub_key in transition[key]:
                total += transition[key][sub_key]
            if key != 'unit_start' and key != 'I0':
                if key[1:] != str(len(alignment[0]) - len(insert_index)):
                    if not 'I' + key[1:] in transition[key]:
                        transition[key]['I' + key[1:]] = 0
                    if not 'D' + str(int(key[1:]) + 1) in transition[key]:
                        transition[key]['D' + str(int(key[1:]) + 1)] = 0
                    if not 'M' + str(int(key[1:]) + 1) in transition[key]:
                        transition[key]['M' + str(int(key[1:]) + 1)] = 0
                else:
                    if not 'I' + key[1:] in transition[key]:
                        transition[key]['I' + key[1:]] = 0
                    if 'unit_end' not in transition[key]:
                        transition[key]['unit_end'] = 0

            for sub_key in transition[key]:
                if total > 0:
                    transition[key][sub_key] = 1.0 * transition[key][sub_key] / total
                    transition[key][sub_key] = (transition[key][sub_key] + pseu) / (1 + pseu * len(transition[key]))
                else:
                    if len(transition[key]) == 3:
                        transition[key][sub_key] = 1.0 / 3
                    if len(transition[key]) == 2:
                        transition[key][sub_key] = 1.0 / 2

    index_list = ['unit_start', 'I0']
    for i in range(1, len(alignment[0])-len(insert_index)+1):
        index_list.extend(['M'+str(i), 'D'+str(i), 'I'+str(i)])
    index_list.append('unit_end')
    for key1 in index_list:
        if key1 not in transition.keys():
            transition[key1] = {}
        for key2 in index_list:
            if key2 not in transition[key1].keys():
                transition[key1][key2] = 0
    return transition, emission


@time_usage
def build_profile_hmm_for_repeats(repeats, error_rate):
    alphabet = 'ACGT'
    pseudocounts = (len(repeats) / 4.0) * (error_rate / 10)
    threshold = 0.5

    muscle_cline = MuscleCommandline('muscle', clwstrict=True)
    data = '\n'.join(['>%s\n' % str(i) + repeats[i] for i in range(len(repeats))])
    stdout, stderr = muscle_cline(stdin=data)
    alignment = AlignIO.read(StringIO(stdout), "clustal")
    aligned_repeats = [str(aligned.seq) for aligned in alignment]

    return build_profile_hmm_pseudocounts_for_alignment(threshold, pseudocounts, alphabet, aligned_repeats)
