import pysam
import os
import settings


def extract_unmapped_reads_to_fasta_file(alignment_file, working_directory='./', use_existing_computed_files=True):
    base_name = os.path.basename(alignment_file).rsplit('.', 1)[0]
    unmapped_bam_file = working_directory + base_name + '.unmapped.bam'
    unmapped_read_file = working_directory + base_name + '.unmapped.fasta'
    fastq_to_fasta_command = "paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n'"
    if not use_existing_computed_files or not os.path.exists(unmapped_bam_file):
        os.system('samtools view -b -f4 %s > %s' % (alignment_file, unmapped_bam_file))
    if not use_existing_computed_files or not os.path.exists(unmapped_read_file):
        os.system('samtools bam2fq %s | %s > %s' % (unmapped_bam_file, fastq_to_fasta_command, unmapped_read_file))
    return unmapped_read_file


def get_reference_genome_of_alignment_file(samfile):
    result = None
    if '1' in samfile.references:
        result = 'GRCh37'
    for reference in samfile.references:
        if reference.startswith('chr'):
            result = 'HG19'
    return result


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


def get_vntr_coverage_over_total_coverage(start_point, end_point, read_file='original_reads/paired_dat.sam'):
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
