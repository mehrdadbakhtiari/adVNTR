import os
import glob
import pysam

from Bio import SeqIO

from reference_vntr import load_unique_vntrs_data
from vntr_finder import VNTRFinder


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


def get_our_selected_reads_count(fq_file, vntr_finder):
    hmm = vntr_finder.get_vntr_matcher_hmm(read_length=150)
    min_score_to_count_read = vntr_finder.get_min_score_to_select_a_read(hmm, None, 150)
    selected_reads = []
    fasta_sequences = SeqIO.parse(open(fq_file), 'fastq')
    for fasta in fasta_sequences:
        logp, vpath = hmm.viterbi(fasta.seq)
        rev_logp, rev_vpath = hmm.viterbi(str(fasta.seq.reverse_complement()))
        logp = max(rev_logp, logp)
        if logp > min_score_to_count_read:
            selected_reads.append(fasta)
    return len(selected_reads)


def get_our_filtered_reads_count(fq_file, vntr_finder):
    unmapped_reads_file = fq_file[:-3] + '.fa'
    os.system("cat %s | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > %s" % (fq_file, unmapped_reads_file))
    filtered_ids = vntr_finder.filter_reads_with_keyword_matching('working_directory/', unmapped_reads_file)
    return len(filtered_ids)


reference_vntrs = load_unique_vntrs_data()
id_to_gene = {119: 'DRD4', 1220: 'GP1BA', 1221: 'CSTB', 1214: 'MAOA', 1219: 'IL1RN'}
genes = glob.glob('../Illumina_copy_number/*')
for gene_dir in genes:
    print(gene_dir)
    files = glob.glob(gene_dir + '/*30x.sam')
    gene_name = gene_dir.split('/')[-1]
    print(len(files))
    mapped_reads = {}
    for file_name in files:
        copies = file_name.split('_')[-2]
        make_bam_and_index(file_name)
        base_name = file_name[:-4]
        bowtie_bam = bowtie_alignment(base_name + '.fq')
        bwa_bam = bwamem_alignment(base_name + '.fq')
        original_bam = file_name[:-4] + '.bam'

        vntr_id = None
        for vid, gname in id_to_gene.items():
            if gname == gene_name:
                vntr_id = vid

        vntr_finder = VNTRFinder(reference_vntrs[vntr_id])
        original = count_reads(original_bam)
        our_selection = get_our_selected_reads_count(base_name + '.fq', vntr_finder)
        our_filtering = get_our_filtered_reads_count(base_name + '.fq', vntr_finder)
        bwa = count_reads_mapped_to_vntr(bwa_bam, reference_vntrs[vntr_id])
        bowtie = count_reads_mapped_to_vntr(bowtie_bam, reference_vntrs[vntr_id])
        mapped_reads[copies] = [original, our_filtering, our_selection, bwa, bowtie]
    with open(gene_dir + '/result.txt', 'w') as out:
        for copies in mapped_reads.keys():
            original, our_filtering, our_selection, bwa, bowtie = mapped_reads[copies]
            out.write('%s %s %s %s %s %s\n' % (copies, original, our_filtering, our_selection, bwa, bowtie))
