from repeat_finder import *
import pysam
from Bio import pairwise2, Seq, SeqRecord, SeqIO
from blast_parameters import *

def get_gc_content(s):
    res = 0
    for e in s:
        if e == 'G' or e == 'C':
            res += 1
    return float(res) / len(s)


def get_read_seq_from_samfile(read_name, read_file):
    result = None
    samfile = pysam.AlignmentFile(read_file, "r")
    for read in samfile.fetch():
        read_number = '/1'
        if read.is_read2:
            read_number = '/2'
        name = read.qname + read_number
        if name == read_name:
            result = read.query
            break
    return result


def get_related_reads_in_samfile(pattern, pattern_start, repeats, read_file):
    pattern_length = len(pattern)
    related_reads = []
    samfile = pysam.AlignmentFile(read_file, "r")
    for read in samfile.fetch():
        start = read.reference_start
        if pattern_start - 150 + pattern_length < int(start) < pattern_start + (repeats-1) * pattern_length:
            read_number = '/1'
            if read.is_read2:
                read_number = '/2'
            name = read.qname + read_number
            related_reads.append(name)
    related_reads = set(related_reads)
    return related_reads


def get_exact_number_of_repeats_from_sequence(pattern, pattern_start):
    file_name = 'chr15.fa'
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
    corresponding_region_in_seq = sequence[pattern_start:pattern_start + len(pattern) * 35].upper()
    repeats = get_occurrence_of_pattern_in_text(corresponding_region_in_seq, pattern, 0.66)
    return repeats


def find_sensitivity(pattern_num, pattern, pattern_start, min_length_for_pattern=None):
    repeats = get_exact_number_of_repeats_from_sequence(pattern, pattern_start)
    related_reads = get_related_reads_in_samfile(pattern, pattern_start, repeats, 'original_reads/paired_dat.sam')

    valid_parameters = read_blast_params_list()
    blast_pattern = pattern
    if min_length_for_pattern and len(pattern) < min_length_for_pattern:
        blast_pattern = pattern * int(round(50.0 / len(pattern) + 0.5))

    for word_size in [4, 7, 10]:
        word_size = str(word_size)
        for params in valid_parameters:
            for evalue in [0.1, 1.0, 10.0]:
                blast_selected_reads = get_blast_matched_ids(blast_pattern, 'original_reads/original_reads', max_seq='6000',
                                                             word_size=word_size, reward=params.reward,
                                                             penalty=params.penalty, gapopen=params.gapopen,
                                                             gapextend=params.gapextend, evalue=evalue)
                TP = [read for read in blast_selected_reads if read in related_reads]
                FP = [read for read in blast_selected_reads if read not in TP]
                FN = [read for read in related_reads if read not in blast_selected_reads]
                print('TP:', len(TP), 'FP:', len(FP), 'blast selected:', len(blast_selected_reads))
                print('FN:', len(FN))
                sensitivity = float(len(TP)) / len(related_reads) if len(related_reads) > 0 else 0
                precision = float(len(TP)) / (len(TP) + len(FP)) if len(TP) + len(FP) > 0 else 0
                min_len = min_length_for_pattern if min_length_for_pattern else 0
                with open('sensivity_over_precision_min_len%s_seq_%s.txt' % (min_len, pattern_num), 'a') as outfile:
                    outfile.write('%s %s\n' % (precision, sensitivity))

    # with open('0_size_related_reads.txt', 'a') as outfile: #0
    #     outfile.write('%s %s\n' % (len(pattern), len(related_reads)))
    # with open('1_size_sensitivity.txt', 'a') as outfile: #1
    #     outfile.write('%s %s\n' % (len(pattern), sensitivity))
    # with open('2_size_blast_selected.txt', 'a') as outfile: #2
    #     outfile.write('%s %s\n' % (len(pattern), len(blast_selected_reads)))
    # with open('3_sim_read_coverage__gc_content.txt', 'a') as outfile: #3
    #     outfile.write('%s %s\n' % (len(related_reads) * 150.0 / (len(pattern) * repeats), get_gc_content(pattern)))


def add_two_copy_to_all_patterns(patterns, start_points):
    file_name = 'chr15.fa'
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    sequence = ''
    record = SeqRecord.SeqRecord('')
    for fasta in fasta_sequences:
        record = fasta
        name, sequence = fasta.id, str(fasta.seq)
    total_added_chars = 0
    for i in range(len(patterns)):
        start_point = start_points[i] + total_added_chars
        pattern = patterns[i]
        sequence = sequence[:start_point] + pattern * 2 + sequence[start_point:]
        total_added_chars += len(pattern * 2)

    record.seq = Seq.Seq(sequence)
    output_name = 'edited_chr15_two_more_copies.fa'
    with open(output_name, 'w') as output_handle:
        SeqIO.write([record], output_handle, 'fasta')


def write_cn_over_true_cn_to_files(patterns, start_points):
    # read_files = [['paired_dat1.fasta', 'paired_dat2.fasta'],
    # ['paired_dat1.fasta', 'paired_dat2.fasta', 'edited_chr15_paired_dat1.fasta', 'edited_chr15_paired_dat2.fasta']]
    # out_files = ['original_computed_cn.txt', 'diploid_computed_cn.txt']
    # directory = ['original_reads/', 'diploid/']

    out_files = ['10X_ratio.txt', '20X_ratio.txt', '30X_ratio.txt']
    read_files = [['10X_paired_dat1.fasta', '10X_paired_dat2.fasta'],
                  ['paired_dat1.fasta', 'paired_dat2.fasta'],
                  ['30X_paired_dat1.fasta', '30X_paired_dat2.fasta']]
    directory = ['10X_reads/', 'original_reads/', '30X_reads/']
    true_cn = []
    for i in range(len(patterns)):
        true_cn.append(get_exact_number_of_repeats_from_sequence(patterns[i], start_points[i]))

    for k in range(len(out_files)):
        db_name = 'blast_db'
        if len(directory[k]):
            db_name = directory[k][:len(directory[k]) - 1]
        blast_db_name = directory[k] + db_name
        for t in range(len(read_files[k])):
            read_files[k][t] = directory[k] + read_files[k][t]
        # make_blast_database(read_files[k], blast_db_name)

        for i in range(len(patterns)):
            calculated_cn = get_copy_number_of_pattern(patterns[i], read_files[k], directory[k])
            with open(out_files[k], 'a') as outfile:
                outfile.write('%s %s\n' % (len(patterns[i]), calculated_cn / true_cn[i]))


with open('patterns.txt') as input:
    patterns = input.readlines()
    patterns = [pattern.strip() for pattern in patterns]
with open('start_points.txt') as input:
    lines = input.readlines()
    start_points = [int(num.strip())-1 for num in lines]

# write_cn_over_true_cn_to_files(patterns, start_points)

for i in range(len(patterns)):
    print(i)
    find_sensitivity(i+1, patterns[i], start_points[i])
