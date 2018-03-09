import logging

from Bio import SeqIO

from advntr.settings import *


def get_min_number_of_copies_to_span_read(pattern, read_length=150):
    return int(round(float(read_length) / len(pattern) + 0.499))


def get_gc_content(s):
    res = 0
    for e in s:
        if e == 'G' or e == 'C':
            res += 1
    return float(res) / len(s)


def is_low_quality_read(read):
    if read.mapq <= MAPQ_CUTOFF:
        logging.debug('Rejecting read for poor mapping quality')
        return True
    low_quality_base_pairs = [i for i, q in enumerate(read.query_qualities) if q < QUALITY_SCORE_CUTOFF]
    if len(low_quality_base_pairs) >= LOW_QUALITY_BP_TO_DISCARD_READ * len(read.query_qualities):
        logging.debug('Rejecting read for having so many low quality base pairs')
        return True
    maximum_low_quality_run = int(LOW_QUALITY_BP_TO_DISCARD_READ * len(read.query_qualities) / 4)
    for i in low_quality_base_pairs:
        passed = False
        for j in range(i+1, i+maximum_low_quality_run):
            if j not in low_quality_base_pairs:
                passed = True
                break
        if not passed:
            logging.debug('Rejecting read for having long run of low quality base pairs')
            return True
    return False


def get_chromosome_reference_sequence(chromosome):
    ref_file_name = HG19_DIR + chromosome + '.fa'
    fasta_sequences = SeqIO.parse(open(ref_file_name), 'fasta')
    ref_sequence = ''
    for fasta in fasta_sequences:
        name, ref_sequence = fasta.id, str(fasta.seq)
    return ref_sequence
