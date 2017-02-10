from acgt_filter import *
from blast_wrapper import *
from Bio.Alphabet import IUPAC
from Bio import pairwise2, Seq, SeqRecord, SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
# import time


def find_exact_match_from_candid_reads(query, candid_reads, min_similarity_score):
    matched_reads = []
    for read_number, read_segment in candid_reads:
        # reward, mismatch, gapopen, gapextend
        alignment_score = pairwise2.align.localms(read_segment, query, 1, -1, -2, -2, score_only=True)
        if float(alignment_score) / len(query) > min_similarity_score:
            matched_reads.append(read_segment)
    return matched_reads


def find_reads_and_estimate_average_coverage(read_files, matched_read_ids, total_length):
    number_of_reads = 0
    read_length = 0
    matched_reads = []
    for read_file in read_files:
        reads = SeqIO.parse(read_file, 'fasta')
        for read_segment in reads:
            if number_of_reads == 0:
                read_length = len(str(read_segment.seq))
            if read_segment.id in matched_read_ids:
                matched_reads.append(str(read_segment.seq))
            number_of_reads += 1
    avg_coverage = float(number_of_reads * read_length) / total_length
    return matched_reads, avg_coverage


def get_occurrence_of_pattern_in_text(text, pattern, min_similarity):
    occurrence = 0
    while True:
        alignments = pairwise2.align.localms(text, pattern, 1, -0.5, -1, -1)
        if len(alignments) <= 0:
            break
        first_alignment = alignments[0]
        score = first_alignment[2]
        if float(score) / len(pattern) < min_similarity or len(text) < len(pattern) * 2 / 3:
            break
        occurrence += 1
        start_index = first_alignment[3]
        end_index = first_alignment[4]
        inserts = 0
        for i in range(start_index, end_index):
            if first_alignment[0][i] == '-':
                inserts += 1
        text = text[:start_index] + text[end_index-inserts:]

    return occurrence


def get_copy_number_of_pattern_in_reads(query, matched_reads, average_coverage=20.0, min_similarity=0.66):
    occurrence = 0
    for read in matched_reads:
        rc_query = str(Seq.Seq(query).reverse_complement())
        query_occurrence = get_occurrence_of_pattern_in_text(read, query, min_similarity)
        rc_query_occurrence = get_occurrence_of_pattern_in_text(read, rc_query, min_similarity)
        occurrence += max(query_occurrence, rc_query_occurrence)
    return float(occurrence) / average_coverage


def get_blast_matched_ids(query, word_size='7', blast_db_name='hg_19_chr_15_reads'):
    with open("query.fasta", "w") as output_handle:
        output_handle.write('query\n')
        output_handle.write(query)
    blastn_cline = NcbiblastnCommandline(query="query.fasta", db=blast_db_name, outfmt='"6 sallseqid"',
                                         out="result.txt", num_threads="4", word_size=word_size, gapopen='2',
                                         gapextend='2', penalty='-1', reward='1')
    blastn_cline()

    with open('result.txt') as result_input:
        ids = result_input.readlines()
        matched_ids = set([seq_id.strip() for seq_id in ids])

    return matched_ids


def get_copy_number_of_pattern(query, fasta_files):
    min_alignment_score = 0.66
    total_length = 100 * 1000 * 1000

    blast_db_name = 'hg_19_chr_15_reads'
    # make_blast_database(fasta_files, blast_db_name)

    matched_ids = get_blast_matched_ids(query, blast_db_name=blast_db_name)
    matched_reads, avg_coverage = find_reads_and_estimate_average_coverage(fasta_files, matched_ids, total_length)
    cn = get_copy_number_of_pattern_in_reads(query, matched_reads, avg_coverage, min_alignment_score)
    return cn


if __name__ == '__main__':
    repeating_unit = 'TAGAACAGAAGGACAAGGCCCTGGAACCAAAAGATAAAGACT'
    paired_end_files = ['paired_dat1.fasta', 'paired_dat2.fasta']
    print(get_copy_number_of_pattern(repeating_unit, paired_end_files))

