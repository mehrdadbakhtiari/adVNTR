import os
import glob
import pysam

from reference_vntr import load_unique_vntrs_data


def count_reads(bam_file):
    print bam_file
    alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    return alignment_file.unmapped + alignment_file.mapped


def count_reads_mapped_to_vntr(bam_file, reference_vntr):
    alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    start = reference_vntr.start_point
    end = reference_vntr.start_point + reference_vntr.get_length()
    reads = []
    for read in alignment_file.fetch(reference_vntr.chromosome, start, end):
        if read.is_secondary or read.is_supplementary:
            continue
        reads.append(read)
    return len(reads)


def make_bam_and_index(samfile):
    bamfile = samfile[:-4]
    os.system('samtools view -bS %s | samtools sort - %s' % (samfile, bamfile))
    os.system('samtools index %s %s' % (bamfile + '.bam', bamfile + '.bai'))


def bowtie_alignment(fq_file):
    bowtie_alignment_file = fq_file[:-3] + '_bowtie_aln.sam'
    os.system('bowtie2 -x hg19_chromosomes/hg19_bt2_idx -U %s -S %s --threads 20' % (fq_file, bowtie_alignment_file))
    make_bam_and_index(bowtie_alignment_file)
    return bowtie_alignment_file[:-4] + '.bam'


def bwamem_alignment(fq_file):
    bwa_alignment_file = fq_file[:-3] + '_bwa_aln.sam'
    os.system('bwa mem -t 20 hg19_chromosomes/CombinedHG19_Reference.fa %s > %s' % (fq_file, bwa_alignment_file))
    make_bam_and_index(bwa_alignment_file)
    return bwa_alignment_file[:-4] + '.bam'


reference_vntrs = load_unique_vntrs_data()
id_to_gene = {119: 'DRD4', 1220: 'GP1BA', 1221: 'CSTB', 1214: 'MAOA', 1219: 'IL1RN'}
genes = glob.glob('../Illumina_copy_number/*')
for gene_dir in genes:
    print(gene_dir)
    files = glob.glob(gene_dir + '/*30x.sam')
    print(len(files))
    mapped_reads = {}
    for file_name in files:
        copies = file_name.split('_')[-2]
        make_bam_and_index(file_name)
        base_name = file_name[:-4]
        bowtie_bam = bowtie_alignment(base_name + '.fq')
        bwa_bam = bwamem_alignment(base_name + '.fq')
        original_bam = file_name[:-4] + '.bam'

        original = count_reads(original_bam)
        bwa = count_reads_mapped_to_vntr(bwa_bam, reference_vntrs[1221])
        bowtie = count_reads_mapped_to_vntr(bowtie_bam, reference_vntrs[1221])
        mapped_reads[copies] = [original, bwa, bowtie]
    with open(gene_dir + '/result.txt', 'w') as out:
        for copies in mapped_reads.keys():
            original, bwa, bowtie = mapped_reads[copies]
            out.write('%s %s %s %s\n' % (copies, original, bwa, bowtie))
    break
