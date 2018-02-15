from distance import *
from Bio import SeqIO, Seq


def match_query_by_sliding_windows(query, query_acgt_content, rc_query_acgt_content, number_of_copies, read_segment):
    read_segment_content = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for i in range(len(read_segment)):
        if i >= len(query):
            read_segment_content[read_segment[i - len(query)].upper()] -= 1
        read_segment_content[read_segment[i].upper()] += 1
        if nucleotide_dist(query_acgt_content, read_segment_content) < 3:
            return 1
        if nucleotide_dist(rc_query_acgt_content, read_segment_content) < 3:
            return 2
    return 0


def get_candid_reads_by_sliding_window_method(query, number_of_copies, fastq_files):
    candid_reads = []
    reversed_complement_query = str(Seq.Seq(query).reverse_complement())
    query_acgt_content = get_nucleotide_map(query)
    rc_query_acgt_content = get_nucleotide_map(reversed_complement_query)
    for fastq_file in fastq_files:
        reads = SeqIO.parse(fastq_file, 'fasta')
        read_counter = 0
        for read_segment in reads:
            match_result = match_query_by_sliding_windows(query, query_acgt_content, rc_query_acgt_content, number_of_copies, str(read_segment.seq))
            if match_result > 0:
                candid_reads.append((read_counter, str(read_segment.seq)))
            read_counter += 1
        print('total reads: ', read_counter)
    return candid_reads


def get_kmers(query, k):
    doubled_query = query + query
    return [doubled_query[i:i+k] for i in range(len(doubled_query)-k+1)]


def has_kmer(kmers, read_segment):
    for kmer in kmers:
        if kmer in read_segment:
            return True
    return False


def get_candid_reads_by_kmer_method(query, number_of_copies, k, fastq_files):
    candid_reads = []
    kmers = get_kmers(query, k)
    for fastq_file in fastq_files:
        reads = SeqIO.parse(fastq_file, 'fastq')
        read_counter = 0
        for read_segment in reads:
            if has_kmer(kmers, str(read_segment.seq)):
                candid_reads.append((read_counter, str(read_segment.seq)))
            read_counter += 1
        print('total reads: ', read_counter)
    return candid_reads
