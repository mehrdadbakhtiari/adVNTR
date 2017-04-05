import pysam


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


def get_related_reads_and_read_count_in_samfile(pattern, pattern_start, repeats=None, read_file='', pattern_end=None):
    pattern_length = len(pattern)
    related_reads = []
    samfile = pysam.AlignmentFile(read_file, "r")
    read_count = 0
    if pattern_end is None:
        pattern_end = pattern_start + (repeats-1) * pattern_length
    for read in samfile.fetch():
        read_count += 1
        start = read.reference_start
        if pattern_start - 150 + pattern_length < int(start) < pattern_end:
            read_number = '/1'
            if read.is_read2:
                read_number = '/2'
            name = read.qname + read_number
            related_reads.append(name)
    related_reads = set(related_reads)
    return related_reads, read_count
