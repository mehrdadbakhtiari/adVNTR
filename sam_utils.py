import pysam
import os
import settings

from profiler import time_usage


@time_usage
def extract_unmapped_reads_to_fasta_file(alignment_file, working_directory='./', use_existing_computed_files=True):
    base_name = os.path.basename(alignment_file).rsplit('.', 1)[0]
    unmapped_bam_file = working_directory + base_name + '.unmapped.bam'
    unmapped_read_file = working_directory + base_name + '.unmapped.fasta'
    fastq_to_fasta_command = "paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n'"
    if not use_existing_computed_files or not os.path.exists(unmapped_read_file):
        if not use_existing_computed_files or not os.path.exists(unmapped_bam_file):
            os.system('samtools view -b -f4 %s | samtools rmdup -S - %s' % (alignment_file, unmapped_bam_file))
        os.system('samtools bam2fq %s | %s > %s' % (unmapped_bam_file, fastq_to_fasta_command, unmapped_read_file))
        os.system('rm -f %s' % unmapped_bam_file)
    return unmapped_read_file


def make_bam_and_index(samfile):
    bamfile = samfile[:-4]
    os.system('samtools view -bS %s | samtools sort - %s' % (samfile, bamfile))
    os.system('samtools index %s %s' % (bamfile + '.bam', bamfile + '.bai'))


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


def get_id_of_reads_mapped_to_vntr_in_bamfile(bam_file, reference_vntr):
    alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    start = reference_vntr.start_point
    end = reference_vntr.start_point + reference_vntr.get_length()
    reads = []
    for read in alignment_file.fetch(reference_vntr.chromosome, start, end):
        if read.is_secondary or read.is_supplementary:
            continue
        reads.append(read.qname)
    return reads


def get_reads_mapped_to_vntr_in_samfile(sam_file, reference_vntr, read_length=150, region=None):
    alignment_file = pysam.AlignmentFile(sam_file, 'r')
    if not region:
        start = reference_vntr.start_point
        end = reference_vntr.start_point + reference_vntr.get_length()
    else:
        start, end = region
    reads = []
    for read in alignment_file.fetch():
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if reference_vntr.chromosome == read.reference_name:
            if start - read_length < read.reference_start < end:
                reads.append(read)
    return reads


def get_id_of_reads_mapped_to_vntr_in_samfile(sam_file, reference_vntr, read_length=150, region=None):
    reads = get_reads_mapped_to_vntr_in_samfile(sam_file, reference_vntr, read_length, region=region)
    return [read.qname for read in reads]


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
