from Bio import Seq, SeqRecord, SeqIO
from utils import get_chromosome_reference_sequence


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


def create_reference_with_specific_repeats(reference_vntr, desired_repeats, output_name):
    record = SeqRecord.SeqRecord('')
    sequence = get_chromosome_reference_sequence(reference_vntr.chromosome)

    record.seq = Seq.Seq(sequence)
    with open(output_name, 'w') as output_handle:
        SeqIO.write([record], output_handle, 'fasta')
