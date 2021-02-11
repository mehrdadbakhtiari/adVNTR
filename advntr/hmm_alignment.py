from advntr.models import load_unique_vntrs_data

from collections import defaultdict
from glob import glob


def get_first_match_state_position(visited_states):
    for state in visited_states:
        if state.startswith("M"):
            return int(state.split("_")[0][1:])


def get_aligned_read(ru_index, ru_sequence, seq, visited_states, mutation):
    seq_index = 0
    aligned_seq = ""

    current_ru = None
    has_mutation_in_current_ru = False
    for state in visited_states:
        if "unit_start_{}".format(ru_index) in state:
            current_ru = state.split("_")[-1]
            has_mutation_in_current_ru = False
        if "unit_end_{}".format(ru_index) in state:
            if not has_mutation_in_current_ru:
                aligned_seq = ""
                continue

            if current_ru is None:  # partial start
                # Pad from the beginning
                start_position = get_first_match_state_position(visited_states)
                start_position = start_position - 1 if start_position > 0 else start_position
                aligned_seq = "-" * start_position + aligned_seq
                # padding_size = len(ru_sequence) - len(aligned_seq)
                # aligned_seq = " " * padding_size + aligned_seq
            return aligned_seq

        if state in mutation:
            has_mutation_in_current_ru = True

        if state.startswith("I") or state.startswith("M"):
            state_ru_index = state.split("_")[-1]
            if state_ru_index == str(ru_index):
                aligned_seq += seq[seq_index]
            seq_index += 1

        if state.startswith("D"):
            aligned_seq += "-"

    # A unit_start has found, but unit_end was not found
    # Padding to the last
    if current_ru == str(ru_index):
        aligned_seq += " " * (len(ru_sequence) - len(aligned_seq))

    return aligned_seq


def is_multiple_mutation(mutation):
    return True if "&" in mutation else False


def get_repeating_unit_state_count(visited_states):
    state_count_for_ru = {}
    complete_ru_index = 0
    prev_start_index = None
    last_end_index = 0
    for i in range(len(visited_states)):
        if visited_states[i].startswith('unit_end'):
            last_end_index = i
            if prev_start_index is None:
                # Haven't seen the start (partially mapped read) Ah only partial start? (error?)
                match_count = 0
                insert_count = 0
                delete_count = 0
                for j in range(0, i):
                    if visited_states[j].startswith("M"):
                        match_count += 1
                    if visited_states[j].startswith("I"):
                        insert_count += 1
                    if visited_states[j].startswith("D"):
                        delete_count += 1
                state_count_for_ru['partial_start'] = {'M': match_count, 'I': insert_count, 'D': delete_count}
            else:
                match_count = 0
                insert_count = 0
                delete_count = 0
                for j in range(prev_start_index, i):
                    if visited_states[j].startswith("M"):
                        match_count += 1
                    if visited_states[j].startswith("I"):
                        insert_count += 1
                    if visited_states[j].startswith("D"):
                        delete_count += 1
                state_count_for_ru[complete_ru_index] = {'M': match_count, 'I': insert_count, 'D': delete_count}
                complete_ru_index += 1
        if visited_states[i].startswith('unit_start'):
            prev_start_index = i

    if last_end_index == 0 and prev_start_index is None:  # When a read is completely within a repeating unit
        match_count = 0
        insert_count = 0
        delete_count = 0
        for j in range(len(visited_states)):
            if visited_states[j].startswith("M"):
                match_count += 1
            if visited_states[j].startswith("I"):
                insert_count += 1
            if visited_states[j].startswith("D"):
                delete_count += 1
        state_count_for_ru['partial_start'] = {'M': match_count, 'I': insert_count, 'D': delete_count}
        return state_count_for_ru

    if prev_start_index is not None:
        if prev_start_index > last_end_index:
            # if met unit start but not unit_end - update the length
            match_count = 0
            insert_count = 0
            delete_count = 0
            for j in range(prev_start_index, len(visited_states)):
                if visited_states[j].startswith("M"):
                    match_count += 1
                if visited_states[j].startswith("I"):
                    insert_count += 1
                if visited_states[j].startswith("D"):
                    delete_count += 1
            state_count_for_ru['partial_end'] = {'M': match_count, 'I': insert_count, 'D': delete_count}

    return state_count_for_ru


def is_matching_state(state):
    return True if state.startswith("M") or state.startswith("I") else False


def get_emitted_basepair_from_visited_states(state, visited_states, sequence):
    if state.startswith("D"):
        raise ValueError("Deletion state doesn't emit base")
    base_pair_idx = 0
    for visited_state in visited_states:
        if visited_state == state:
            break
        if is_matching_state(visited_state):
            base_pair_idx += 1
    return sequence[base_pair_idx]


def get_modified_base_count(detected_mutation):
    insertion_count = int(detected_mutation[-1])  # IX_X_LENX
    return insertion_count
    # deletion_count = detected_mutation.count("D")
    # if insertion_count > deletion_count:
    #     return insertion_count - deletion_count
    # return deletion_count - insertion_count


def generate_aln(advntr_logfile, output_mutations=None, out_folder="", reference_vntr_db=None, ref_vntr_dict=None):
    if reference_vntr_db is None and ref_vntr_dict is None:
        raise ValueError("Either reference DB or dictionary should be given")
    if reference_vntr_db is not None:
        ref_vntrs = load_unique_vntrs_data(reference_vntr_db)
        ref_vntr_dict = {ref_vntr.id: ref_vntr for ref_vntr in ref_vntrs}

    out = open(out_folder + advntr_logfile + ".aln", 'w')

    # Find mutations
    detected_mutations = set()
    with open(advntr_logfile, "r") as f:
        for line in f:
            if "There is a mutation at" in line:
                detected_mutation = line.strip().split(" ")[-1]
                detected_mutations.add(detected_mutation)

    target_mutations = detected_mutations
    if output_mutations:
        target_mutations = output_mutations.intersection(detected_mutations)
        if len(target_mutations) == 0:
            raise ValueError("There is no intersection with the target mutations")
    out.write("Target mutations {}\n".format(target_mutations))

    # Find related sequences
    vid_read_length = defaultdict(int)
    vid_to_aln_info = defaultdict(lambda: defaultdict(list))
    with open(advntr_logfile, "r") as f:
        for line in f:
            if "INFO:Using read length" in line:  # INFO:Using read length [read_length]
                read_length = int(line.split(" ")[-1])
                vid_read_length[vid] = read_length
            if "DEBUG:finding" in line:  # DEBUG:finding repeat count from alignment file for [vid]
                vid = int(line.split(" ")[-1])
                reference_vntr = ref_vntr_dict[vid]
                patterns = reference_vntr.get_repeat_segments()
                pattern_clusters = [[pattern] * patterns.count(pattern) for pattern in set(sorted(patterns))]
                observed_mutations = set()
                observed_prefix_suffix_mutations = set()

            if "DEBUG:Read:" in line:  # DEBUG:Read:[sequence]
                sequence = line[30 + 5:].strip()
            if "DEBUG:VisitedStates:" in line:  # DEBUG:VisitedStates:['state1', ...]
                visited = line[line.index('[') + 1:-2]
                split = visited.split(', ')
                split = [item[1:-1] for item in split]  # only take the state inside ' '
                visited_states = split

                # for visited_state in visited_states:
                #     if visited_state.startswith("I") or visited_state.startswith("D"):
                #         if visited_state in mutations:
                #             vid_to_aln_info[vid][visited_state].append((sequence, visited_states))

                ru_state_count = get_repeating_unit_state_count(visited_states)
                fully_observed_ru_count = len(ru_state_count)
                if 'partial_start' in ru_state_count:
                    fully_observed_ru_count -= 1
                if 'partial_end' in ru_state_count:
                    fully_observed_ru_count -= 1

                current_repeat = None
                is_valid_read = True
                reason_why_rejected = ""

                # Keep all mutations in a read, update them only if the read is valid
                mutation_count_temp = defaultdict(int)
                prefix_suffix_mutation_count_temp = defaultdict(int)

                prefix_match_count = 0
                prefix_mutation_count = 0
                suffix_match_count = 0
                suffix_mutation_count = 0

                for i in range(len(visited_states)):
                    current_state = visited_states[i]

                    if current_state.startswith('unit_start'):
                        if not is_valid_read:
                            break
                        if current_repeat is None:
                            current_repeat = 0
                        else:
                            current_repeat += 1

                    if current_state.endswith('fix'):  # Save all mutations observed in prefix or suffix
                        if current_state.startswith('I') or current_state.startswith('D'):
                            prefix_suffix_mutation_count_temp[current_state] += 1
                        if current_state.endswith('prefix'):
                            if current_state.startswith('M'):
                                prefix_match_count += 1
                            else:
                                prefix_mutation_count += 1
                        else:
                            if current_state.startswith('M'):
                                suffix_match_count += 1
                            else:
                                suffix_mutation_count += 1
                        continue

                    if not current_state.startswith('I') and not current_state.startswith('D'):
                        continue

                    # Reads starting with a partially observed repeat unit
                    if current_repeat is None:
                        if 'partial_start' in ru_state_count:
                            if ru_state_count['partial_start']['M'] < 5:
                                continue
                            if ru_state_count['partial_start']['I'] != ru_state_count['partial_start']['D']:
                                if current_state.startswith('I'):
                                    current_state += '_' + get_emitted_basepair_from_visited_states(current_state,
                                                                                                visited_states,
                                                                                                sequence)
                            mutation_count_temp[current_state] += 1
                            continue

                    # Reads ending with a partially observed repeat unit
                    if current_repeat >= fully_observed_ru_count:
                        if 'partial_end' in ru_state_count:
                            if ru_state_count['partial_end']['M'] < 5:
                                continue
                            if ru_state_count['partial_end']['I'] != ru_state_count['partial_end']['D']:
                                if current_state.startswith('I'):
                                    current_state += '_' + get_emitted_basepair_from_visited_states(current_state,
                                                                                                    visited_states,
                                                                                                    sequence)
                                mutation_count_temp[current_state] += 1
                            continue

                    pattern_index = current_state.split('_')[-1]

                    # ru_state_count is a dictionary of [repeat][M/I/D]
                    # This check is okay because insertion and deletion at a different position in a RU is very rare
                    if ru_state_count[current_repeat]['I'] == ru_state_count[current_repeat]['D']:
                        continue

                    pattern_length = len(pattern_clusters[int(pattern_index) - 1][0])
                    inserted_bp = abs(
                        ru_state_count[current_repeat]['M'] + ru_state_count[current_repeat]['I'] - pattern_length)

                    if inserted_bp > pattern_length / 2:
                        reason_why_rejected = "Rejected read: #M + #I - len(pattern) > {} bp in pattern {}, inserted {} bps".format(
                            pattern_length / 2, pattern_index, inserted_bp)
                        is_valid_read = False
                        break

                    indel_mutation_count = ru_state_count[current_repeat]['I'] + ru_state_count[current_repeat]['D']
                    if indel_mutation_count > pattern_length / 2:  # or indel_mutation_count > 4:  # Too many erroneous reads XXX
                        reason_why_rejected = "Rejected read: #I + #D > {} in pattern {}, diff {}".format(
                            pattern_length / 2, pattern_index, indel_mutation_count)
                        is_valid_read = False
                        break

                    # TODO If there are run of insertions, the sequence should be different
                    if current_state.startswith('I'):
                        current_state += '_' + get_emitted_basepair_from_visited_states(current_state, visited_states,
                                                                                        sequence)

                    mutation_count_temp[current_state] += 1

                if is_valid_read and len(mutation_count_temp) > 0:
                    # Check if the mutations are adjacent each other
                    # mutation_count_temp is the dictionary of mutations in a repeat unit
                    if len(mutation_count_temp) == 1:  # only one mutation
                        temp_mutation, count = mutation_count_temp.popitem()
                        if temp_mutation.startswith("D"):
                            if temp_mutation in target_mutations:
                                vid_to_aln_info[vid][temp_mutation].append((sequence, visited_states))
                        else:  # Insertion
                            temp_mutation = temp_mutation + "_LEN{}".format(count)
                            if temp_mutation in target_mutations:
                                vid_to_aln_info[vid][temp_mutation].append((sequence, visited_states))
                    else:
                        sorted_temp_mutations = sorted(mutation_count_temp.items(), key=lambda x: x[0])
                        prev_mutation = sorted_temp_mutations[0][0]
                        mutation_sequence = prev_mutation
                        if prev_mutation.startswith("I"):
                            mutation_sequence += "_LEN{}".format(sorted_temp_mutations[0][1])

                        for i in range(1, len(sorted_temp_mutations)):
                            temp_mutation = sorted_temp_mutations[i][0]
                            current_mutation_index = int(temp_mutation.split("_")[0][1:])
                            current_ru_index = int(temp_mutation.split("_")[1])
                            prev_mutation_index = int(prev_mutation.split("_")[0][1:])
                            prev_ru_index = int(prev_mutation.split("_")[1])

                            if current_ru_index != prev_ru_index:
                                if mutation_sequence is not None:  # Prev mutation was a deletion
                                    if mutation_sequence in target_mutations:
                                        vid_to_aln_info[vid][mutation_sequence].append((sequence, visited_states))
                                prev_mutation = temp_mutation
                                continue

                            if temp_mutation.startswith("D"):
                                # Case 1: D(i-1), D(i),
                                # In this case, the deletion is connected to the previous mutation sequence and skip
                                if prev_mutation_index + 1 == current_mutation_index:  # Only possible with D(i-1)
                                    mutation_sequence += '&' + temp_mutation
                                # Case 2: I/D(j), D(i), j < i-1
                                # In this case, they are not connected (This should be rare, two separated deletions in a RU)
                                else:
                                    # Save the previous mutation and initialize it
                                    if mutation_sequence is not None:  # Prev mutation was a deletion
                                        if mutation_sequence in target_mutations:
                                            vid_to_aln_info[vid][mutation_sequence].append((sequence, visited_states))
                                    mutation_sequence = temp_mutation

                            if temp_mutation.startswith("I"):
                                # Case 3: D(i), I(i)
                                if prev_mutation_index == current_mutation_index:  # Only possible with D(i)
                                    # Add the insertion and done
                                    mutation_sequence += "&{}_LEN{}".format(temp_mutation, mutation_count_temp[temp_mutation])
                                    if mutation_sequence in target_mutations:
                                        vid_to_aln_info[vid][mutation_sequence].append((sequence, visited_states))
                                    mutation_sequence = None
                                # Case 4: I/D(j), I(i), j < i-1
                                else:
                                    # Save the previous mutation and initialize it
                                    if mutation_sequence is not None:  # Prev mutation was a deletion
                                        if mutation_sequence in target_mutations:
                                            vid_to_aln_info[vid][mutation_sequence].append((sequence, visited_states))
                                    mutation_sequence = "{}_LEN{}".format(temp_mutation, mutation_count_temp[temp_mutation])
                                    if mutation_sequence in target_mutations:
                                        vid_to_aln_info[vid][mutation_sequence].append((sequence, visited_states))
                                    mutation_sequence = None

                            prev_mutation = temp_mutation

                        # Last check (e.g. D1-D2)
                        if mutation_sequence is not None:
                            if mutation_sequence in target_mutations:
                                vid_to_aln_info[vid][mutation_sequence].append((sequence, visited_states))

                    # for state in mutation_count_temp:
                    #     occurrence = mutation_count_temp[state]
                    #     if state.startswith('I'):
                    #         state += "_LEN{}".format(occurrence)  # Insertion length
                    #     mutations[state] += 1

                    # Update only when the pre-/suffix match rate is > 0.9
                    for state in prefix_suffix_mutation_count_temp.keys():
                        if state.endswith('prefix'):
                            if prefix_match_count != 0:
                                if prefix_mutation_count / float(prefix_match_count) < 0.9:
                                    continue
                        else:
                            if suffix_match_count != 0:
                                if suffix_mutation_count / float(suffix_match_count) < 0.9:
                                    continue

                        occurrence = prefix_suffix_mutation_count_temp[state]
                        if state.startswith('I'):
                            state += "_LEN{}".format(occurrence)  # Insertion length
                        if state in target_mutations:
                            vid_to_aln_info[vid][state].append((sequence, visited_states))

    for vid in vid_to_aln_info:
        for detected_mutation in vid_to_aln_info[vid]:
            out.write("Mutation {}\n".format(detected_mutation))
            ru_index = int(detected_mutation.split("&")[-1].split("_")[1])
            ru_sequence = sorted(list(set(ref_vntr_dict[vid].repeat_segments)))[ru_index-1]

            if is_multiple_mutation(detected_mutation):
                modified_base_count = get_modified_base_count(detected_mutation)
                if modified_base_count > 0:  # Insertion
                    mutation_position = 0
                    for mutation in detected_mutation.split("&"):
                        if mutation.startswith("I"):
                            mutation_position = int(mutation.split("_")[0][1:])
                            break
                    ru_sequence = ru_sequence[:mutation_position] + "-"*modified_base_count + ru_sequence[mutation_position:]
            else:
                if detected_mutation.startswith("I"):
                    mutation_position = int(detected_mutation.split("_")[0][1:])
                    modified_base_count = get_modified_base_count(detected_mutation)
                    ru_sequence = ru_sequence[:mutation_position] + "-"*modified_base_count + ru_sequence[mutation_position:]

            out.write("Reference repeat unit sequence\n")
            out.write("RefRU : {}{}".format(ru_sequence, "\n"))

            seq_state_pairs = vid_to_aln_info[vid][detected_mutation]
            aligned_reads = []
            for seq_state_pair in seq_state_pairs:
                aligned_reads.append(get_aligned_read(ru_index, ru_sequence, seq_state_pair[0], seq_state_pair[1], detected_mutation))

            for i, aligned_read in enumerate(sorted(aligned_reads)):
                out.write("Read{:>2}: {}{}".format(i, aligned_read, "\n"))


def get_samples_having_muc1_insertion(input_files):
    selected_samples = []
    for input_file in input_files:
        with open(input_file, "r") as f:
            for line in f:
                if "I22_2_G_LEN1" in line:
                    selected_samples.append(input_file)
                    break
    return selected_samples


if __name__ == "__main__":
    logfiles = glob("*.log")
    print("log file count", len(logfiles))

    reference_vntr_db = "hg38_VNTRs_by_TRF.db"
    ref_vntrs = load_unique_vntrs_data(reference_vntr_db)
    ref_vntr_dict = {ref_vntr.id: ref_vntr for ref_vntr in ref_vntrs}

    for logfile in logfiles:
        sample_id = logfile.split("/")[-1].split(".")[0].split("_")[1]
        print(sample_id)
        generate_aln(logfile, None, "", None, ref_vntr_dict)
