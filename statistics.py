from repeat_finder import *
import pysam


def get_gc_content(s):
    res = 0
    for e in s:
        if e == 'G' or e == 'C':
            res += 1
    return float(res) / len(s)


def get_related_reads(pattern, pattern_start, repeats, read_file):
    pattern_length = len(pattern)
    related_reads = []
    samfile = pysam.AlignmentFile(read_file, "r")
    for read in samfile.fetch():
        start = read.reference_start
        read_number = '/1'
        if read.is_read2:
            read_number = '/2'
        name = read.qname + read_number
        if pattern_start - 150 + pattern_length < int(start) < pattern_start + (repeats-1) * pattern_length:
            related_reads.append(name)
    related_reads = set(related_reads)
    return related_reads


def find_sensibility(pattern, pattern_start):
    # pattern_start = 89398653
    # pattern = 'CTCTGCCCCTGGAGTAGAGGACATCAGCGGGCTTCCTTCTGGAGAAGTTCTAGAGAC'
    file_name = 'chr15.fa'
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
    repeats = get_occurrence_of_pattern_in_text(sequence[pattern_start:pattern_start+len(pattern)*35].upper(), pattern, 0.66)

    related_reads = get_related_reads(pattern, pattern_start, repeats, 'paired_dat.sam')
    blast_selected_reads = get_blast_matched_ids(pattern)
    correctly_filtered_reads = [read for read in blast_selected_reads if read in related_reads]
    sensibility = float(len(correctly_filtered_reads)) / len(related_reads)
    with open('1_size_sensibility.txt', 'a') as outfile: #1
        outfile.write('%s %s\n' % (len(pattern), sensibility))
    with open('2_size_blast_selected.txt', 'a') as outfile: #2
        outfile.write('%s %s\n' % (len(pattern), len(blast_selected_reads)))
    with open('3_sim_read_coverage__gc_content.txt', 'a') as outfile: #3
        outfile.write('%s %s\n' % (len(related_reads) * 150.0 / (len(pattern) * repeats), get_gc_content(pattern)))


with open('patterns.txt') as input:
    patterns = input.readlines()
    patterns = [pattern.strip() for pattern in patterns]
with open('start_points.txt') as input:
    lines = input.readlines()
    start_points = [int(num.strip())-1 for num in lines]

for i in range(len(patterns)):
    print(i)
    find_sensibility(patterns[i], start_points[i])
