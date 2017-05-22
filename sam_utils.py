import pysam
import settings


def get_read_seq_from_samfile(read_name, read_file='original_reads/paired_dat.sam'):
    read = get_read_from_samfile(read_name, read_file)
    return read.query


def get_reads_seq_from_samfile(read_names, read_file='original_reads/paired_dat.sam'):
    reads = get_reads_from_samfile(read_names, read_file)
    result = [(name, read.query) for name, read in reads]
    return result


def get_read_from_samfile(read_name, read_file='original_reads/paired_dat.sam'):
    result = None
    samfile = pysam.AlignmentFile(read_file, "r")
    for read in samfile.fetch():
        read_number = '/1'
        if read.is_read2:
            read_number = '/2'
        name = read.qname + read_number
        if name == read_name:
            result = read
            break
    return result


def get_reads_from_samfile(read_names, read_file='original_reads/paired_dat.sam'):
    result = []
    samfile = pysam.AlignmentFile(read_file, "r")
    for read in samfile.fetch():
        read_number = '/1'
        if read.is_read2:
            read_number = '/2'
        name = read.qname + read_number
        for read_name in read_names:
            if name == read_name:
                result.append((name, read))
    return result


def get_related_reads_and_read_count_in_samfile(pattern, pattern_start, repeats=None, read_file='', pattern_end=None):
    pattern_length = len(pattern)
    related_reads = []
    samfile = pysam.AlignmentFile(read_file, "r")
    read_count = 0
    if pattern_end is None:
        pattern_end = pattern_start + repeats * pattern_length
    for read in samfile.fetch():
        read_count += 1
        start = read.reference_start
        offset = 0
        if pattern_start - 150 + offset < int(start) < pattern_end - offset:
            read_number = '/1'
            if read.is_read2:
                read_number = '/2'
            name = read.qname + read_number
            related_reads.append(name)
    related_reads = set(related_reads)
    return related_reads, read_count


def get_VNTR_coverage_over_total_coverage(start_point, end_point, read_file='original_reads/paired_dat.sam'):
    samfile = pysam.AlignmentFile(read_file, "r")
    read_count = 0
    avg_read_size = 0
    vntr_bp_in_reads = 0
    for read in samfile.fetch():
        avg_read_size = (len(read.seq) + avg_read_size * read_count) / (read_count + 1)
        read_count += 1
        if start_point <= read.reference_start < end_point or start_point < read.reference_end <= end_point:
            end = min(read.reference_end, end_point)
            start = max(read.reference_start, start_point)
            vntr_bp_in_reads += end - start
    total_coverage = float(read_count * avg_read_size) / settings.GENOME_LENGTH
    vntr_coverage = vntr_bp_in_reads / float(end_point - start_point)
    return vntr_coverage / total_coverage
