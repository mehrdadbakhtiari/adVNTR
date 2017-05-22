from Bio import pairwise2, Seq, SeqIO
from blast_wrapper import get_blast_matched_ids
import settings


def find_exact_match_from_candid_reads(query, candid_reads, min_similarity_score):
    matched_reads = []
    for read_number, read_segment in candid_reads:
        # reward, mismatch, gapopen, gapextend
        penalty = -1
        if len(query) < 25:
            penalty = -3
        alignment_score = pairwise2.align.localms(read_segment, query, 1, penalty, -2, -2, score_only=True)
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


def get_number_of_occurrence_of_pattern_in_text(text, pattern, min_similarity, get_sequences=False):
    occurrence = 0
    sequences = []
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
        if get_sequences:
            sequences.append(text[start_index:end_index-inserts])

    if get_sequences:
        return occurrence, sequences
    else:
        return occurrence


def get_copy_number_of_pattern_in_reads(query, matched_reads, average_coverage=20.0, min_similarity=0.66):
    occurrence = 0
    for read in matched_reads:
        rc_query = str(Seq.Seq(query).reverse_complement())
        query_occurrence = get_number_of_occurrence_of_pattern_in_text(read, query, min_similarity)
        rc_query_occurrence = get_number_of_occurrence_of_pattern_in_text(read, rc_query, min_similarity)
        occurrence += max(query_occurrence, rc_query_occurrence)
    return float(occurrence) / average_coverage


def get_copy_number_of_pattern(query, fasta_files, directory='', min_len=None):
    min_alignment_score = 0.66
    total_length = settings.GENOME_LENGTH

    db_name = 'blast_db'
    if len(directory):
        db_name = directory[:len(directory)-1]
    blast_db_name = directory + db_name
    # make_blast_database(fasta_files, blast_db_name)

    if min_len and len(query) < min_len:
        blast_pattern = query * int(round(50.0 / len(query) + 0.5))
    else:
        blast_pattern = query
    matched_ids = get_blast_matched_ids(blast_pattern, blast_db_name=blast_db_name, evalue=1.0)
    matched_reads, avg_coverage = find_reads_and_estimate_average_coverage(fasta_files, matched_ids, total_length)
    cn = get_copy_number_of_pattern_in_reads(query, matched_reads, avg_coverage, min_alignment_score)
    return cn


if __name__ == '__main__':
    repeating_unit = 'TAGAACAGAAGGACAAGGCCCTGGAACCAAAAGATAAAGACT'
    paired_end_files = ['paired_dat1.fasta', 'paired_dat2.fasta']
    files_dir = 'original_reads/'
    for i in range(len(paired_end_files)):
        paired_end_files[i] = files_dir + paired_end_files[i]
    print(get_copy_number_of_pattern(repeating_unit, paired_end_files, files_dir))

