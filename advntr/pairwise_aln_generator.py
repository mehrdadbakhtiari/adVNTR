from advntr.models import load_unique_vntrs_data

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

from collections import defaultdict
from collections import Counter

import glob
import os.path
import argparse


def get_consensus_pattern(patterns):
    aligned_patterns = None
    if len(patterns) > 1:
        muscle_cline = MuscleCommandline('muscle', clwstrict=True)
        data = '\n'.join(['>%s\n' % str(i) + patterns[i] for i in range(len(patterns))])
        stdout, stderr = muscle_cline(stdin=data)
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_patterns = [str(aligned.seq) for aligned in alignment]
    else:
        aligned_patterns = patterns

    assert len(aligned_patterns[0]) == len(aligned_patterns[-1])
    consensus_seq = ""
    for seq_idx in range(len(aligned_patterns[0])):
        chr_counter = Counter()
        for pattern_idx in range(len(aligned_patterns)):
            chr_counter[aligned_patterns[pattern_idx][seq_idx]] += 1

        most_common_chr = chr_counter.most_common(1)[0][0]
        if most_common_chr == '-':
            if chr_counter.most_common(1)[0][1] < len(aligned_patterns) * 0.5:
                consensus_seq += chr_counter.most_common(2)[1][0]
        else:
            consensus_seq += most_common_chr

    return consensus_seq


def write_alignment(af, vntr_id, repeat_seq_dict, ref_vntr, read_length=151):
    query_id = "VID:" + str(vntr_id) + " " \
               + "REFRC:" + str(ref_vntr.estimated_repeats) + " " \
               + ref_vntr.chromosome + ":" \
               + str(ref_vntr.start_point) + "-" + str(ref_vntr.start_point + ref_vntr.get_length())

    left_flank = ref_vntr.left_flanking_region
    right_flank = ref_vntr.right_flanking_region

    patterns = ref_vntr.repeat_segments
    consensus_pattern = get_consensus_pattern(patterns)

    for repeat in sorted(repeat_seq_dict.keys()):
        seq_list = repeat_seq_dict[repeat]
        for idx, (sequence, visited_states) in enumerate(seq_list):
            # Alignment header
            af.write(">" + str(idx) + "_RC:" + str(repeat) + "_SEQLEN:" + str(len(sequence)) + " " + query_id + "\n")

            query_seq = ""
            ref_seq = ""
            match_line = ""

            left_flank_bp_count = 0
            right_flank_bp_count = 0
            mismatch_count_in_left_flanking = 0
            mismatch_count_in_right_flanking = 0
            mismatch_count_in_flanking_region = 0
            seq_index = 0

            max_hmm_index = -1
            prev_state = visited_states[0]  # initialization
            for state in visited_states:
                if 'suffix_end_suffix' in state:
                    max_hmm_index = int(prev_state.split("_")[0][1:])
                    break
                prev_state = state

            for state in visited_states:
                if 'start' in state:
                    if 'unit_start' in state:
                        query_seq += "|"
                        ref_seq += "|"
                        match_line += "+"
                    if 'Prefix Matcher HMM' in state:
                        query_seq += "*"
                        ref_seq += "*"
                        match_line += ">"
                    continue
                if 'end' in state:
                    if 'Suffix Matcher HMM' in state:
                        query_seq += "*"
                        ref_seq += "*"
                        match_line += "<"
                    if 'Repeat Matcher HMM' in state:
                        query_seq += "|"
                        ref_seq += "|"
                        match_line += "+"
                    continue

                split = state.split("_")  # [M/I/D][hmm_index]_[unit_index]
                state_name = split[0][0]  # [M/I/D]
                hmm_index = int(split[0][1:])
                if 'suffix' in state:
                    left_flank_bp_count += 1
                    if state_name == "M":
                        query_seq += sequence[seq_index]
                        ref_seq += left_flank[-(max_hmm_index - hmm_index + 1)]
                        if sequence[seq_index] == left_flank[-(max_hmm_index - hmm_index + 1)]:
                            match_line += "|"
                        else:
                            match_line += " "
                            mismatch_count_in_left_flanking += 1
                        seq_index += 1
                    if state_name == "I":
                        query_seq += sequence[seq_index]
                        ref_seq += '-'
                        match_line += " "
                        mismatch_count_in_left_flanking += 1
                        seq_index += 1
                    if state_name == "D":
                        query_seq += '-'
                        ref_seq += left_flank[-(max_hmm_index - hmm_index + 1)]
                        match_line += " "
                        mismatch_count_in_flanking_region += 1

                elif 'prefix' in state:
                    right_flank_bp_count += 1
                    if state_name == "M":
                        query_seq += sequence[seq_index]
                        ref_seq += right_flank[hmm_index - 1]
                        if sequence[seq_index] == right_flank[hmm_index - 1]:
                            match_line += "|"
                        else:
                            match_line += " "
                            mismatch_count_in_right_flanking += 1
                        seq_index += 1
                    if state_name == "I":
                        query_seq += sequence[seq_index]
                        ref_seq += '-'
                        match_line += " "
                        mismatch_count_in_right_flanking += 1
                        seq_index += 1
                    if state_name == "D":
                        query_seq += '-'
                        ref_seq += right_flank[hmm_index - 1]
                        match_line += " "
                        mismatch_count_in_right_flanking += 1

                else:  # Pattern matches
                    if state_name == "M":
                        query_seq += sequence[seq_index]
                        ref_seq += consensus_pattern[hmm_index - 1]
                        if sequence[seq_index] == consensus_pattern[hmm_index - 1]:
                            match_line += "|"
                        else:
                            match_line += " "
                        seq_index += 1
                    if state_name == "I":
                        query_seq += sequence[seq_index]
                        ref_seq += '-'
                        match_line += " "
                        seq_index += 1
                    if state_name == "D":
                        query_seq += '-'
                        ref_seq += consensus_pattern[hmm_index - 1]
                        match_line += " "

            af.write(query_seq + "\n")
            af.write(match_line + "\n")
            af.write(ref_seq + "\n")
            af.write("# Mismatch in flanking regions: {}/{} {:.2f}, L:{}/{} {:.2f}, R:{}/{} {:.2f}\n".format(\
                mismatch_count_in_left_flanking + mismatch_count_in_right_flanking,\
                left_flank_bp_count + right_flank_bp_count, \
                float(mismatch_count_in_left_flanking + mismatch_count_in_right_flanking)\
                / (left_flank_bp_count + right_flank_bp_count),\

                mismatch_count_in_left_flanking,\
                left_flank_bp_count,\
                float(mismatch_count_in_left_flanking)/left_flank_bp_count,\

                mismatch_count_in_right_flanking,\
                right_flank_bp_count,\
                float(mismatch_count_in_right_flanking)/right_flank_bp_count))


def _generate_pairwise_aln(log_file, aln_outfile, ref_vntrs, vid_list=None, sort_by_repeat=True):
    vid_to_aln_info = defaultdict(lambda: defaultdict(list))
    vid_read_length = defaultdict(int)

    is_target = False if vid_list is not None else True

    # Load adVNTR log files
    with open(log_file, "r") as f:
        for line in f:
            if "INFO:Using read length" in line:  # INFO:Using read length [read_length]
                read_length = int(line.split(" ")[-1])
                vid_read_length[vid] = read_length
            if "DEBUG:finding" in line:  # DEBUG:finding repeat count from alignment file for [vid]
                vid = int(line.split(" ")[-1])
                if vid_list is not None:
                    if vid in vid_list:
                        is_target = True
                    else:
                        is_target = False
            if is_target:
                if "DEBUG:repeats" in line:  # DEBUG:repeats: [repeat_count]
                    repeats = int(line.strip().split(" ")[-1])
                    vid_to_aln_info[vid][repeats].append((sequence, visited_states))
                if "DEBUG:spanning read visited states" in line or "DEBUG:[" in line:
                    visited = line[line.index('[') + 1:-2]
                    split = visited.split(', ')
                    split = [item[1:-1] for item in split]
                    visited_states = split
                else:
                    sequence = line[30:].strip()

    if sort_by_repeat:
        vid_to_aln_info = {k: v for k, v in sorted(vid_to_aln_info.items(), key=lambda (k, v): v)}

    with open(aln_outfile, "w") as af:
        for vid in sorted(vid_to_aln_info.keys()):
            af.write("#VID: {}\n".format(vid))
            if ref_vntrs[vid] is None:
                print("ERROR: The reference VNTR is not in the DB, VID: {}".format(vid))
                af.write("ERROR: The reference VNTR is not in the DB, VID: {}\n".format(vid))
                continue
            write_alignment(af, vid, vid_to_aln_info[vid], ref_vntrs[vid], vid_read_length[vid])


def generate_pairwise_aln(log_file, aln_file, ref_vntr_db=None, vntr_ids=None, sort_by_repeat=True):
    """
    Generate pairwise alignment for each spanning reads
    :param log_file: a log file or a directory
    :param aln_file: output file name
    :param ref_vntr_db: reference VNTR database
    :param vntr_ids: VNTR id list that you want to generate alignment
    :param sort_by_repeat: if True, the reads will be sorted by the number of repeats
    """
    # Load reference VNTRs
    reference_vntrs = load_unique_vntrs_data(ref_vntr_db)
    ref_vntrs = {ref_vntr.id: ref_vntr for ref_vntr in reference_vntrs}

    if os.path.isdir(log_file):
        log_files = glob.glob(log_file + "/log_*.log")
        for lf in log_files:
            if aln_file is not None:
                print("ERROR: If log file is given as a directory, output name should be None")
                exit(-1)
            if aln_file is None:
                aln_file = log_file.split("/")[-1].split(".")[0] + ".aln"
            _generate_pairwise_aln(lf, aln_file, ref_vntrs, vntr_ids, sort_by_repeat)
    else:
        if aln_file is None:
            aln_file = log_file.split("/")[-1].split(".")[0] + ".aln"
        _generate_pairwise_aln(log_file, aln_file, ref_vntrs, vntr_ids, sort_by_repeat)


def init_argparse():
    parser = argparse.ArgumentParser(
        usage="%(prog)s [-i logfile] [-o outfile] [-db vntr_db] [-vid vntr_ids]",
        description="Generate a pairwise alignment for all spanning reads identified by adVNTR."
    )
    parser.add_argument(
        "-i",
        action='store',
        help="A logfile or directory where a logfile or logfiles are located",
        required=True
    )
    parser.add_argument(
        "-o",
        action='store',
        default=None,
        help="Output file name. If not specified, .aln file is created"
    )
    parser.add_argument(
        "-db",
        action='store',
        default=None,
        help="Reference VNTRs DB. If not specified, the DB in the setting.py is used"
    )
    parser.add_argument(
        "-vid",
        action='store',
        nargs="*",
        default=None,
        help="VID list to generate pairwise alignment. If not specified, all VNTRs in the log file is the list"
    )
    return parser


def main():
    parser = init_argparse()
    args = parser.parse_args()
    arg_dict = (vars(args))
    target_vntrs_ids = None
    if arg_dict['vid'] is not None:
        target_vntrs_ids = [int(x) for x in arg_dict['vid']]
    generate_pairwise_aln(arg_dict['i'], arg_dict['o'], arg_dict['db'], vntr_ids=target_vntrs_ids, sort_by_repeat=True)


if __name__ == "__main__":
    main()