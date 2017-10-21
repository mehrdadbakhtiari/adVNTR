import glob
import os
import pysam
import sys

from reference_vntr import load_unique_vntrs_data
from sam_utils import get_id_of_reads_mapped_to_vntr_in_samfile


def clean_up_tmp():
    os.system('rm -rf /tmp/*.sam')
    os.system('rm -rf /tmp/*.fasta')


def bowtie_alignment(fasta_file, output, param):
    os.system('bowtie2 -x hg19_chromosomes/hg19_bt2_idx --end-to-end -f %s -S %s --threads 24 --score-min L,-0.6,%s' % (fasta_file, output, param))


def bwa_alignment(fasta_file, output, param):
    os.system('bwa mem -T %s -t 24 hg19_chromosomes/CombinedHG19_Reference.fa %s > %s' % (param, fasta_file, output))


def save_reads_stat(file_name, reads):
    with open(file_name, 'w') as out:
        for read in reads:
            alignment_score = None
            for key, value in read.tags:
                if key == 'AS':
                    alignment_score = value
            out.write('%s %s\n' % (read.qname, alignment_score))


def get_positive_and_fn_reads_from_samfile(sam_file, reference_vntr, true_reads, read_length=150):
    alignment_file = pysam.AlignmentFile(sam_file, 'r')
    start = reference_vntr.start_point
    end = reference_vntr.start_point + reference_vntr.get_length()
    positive_reads = []
    false_negative_reads = []
    for read in alignment_file.fetch():
        if read.is_unmapped:
            if read.qname in true_reads:
                false_negative_reads.append(read)
            continue
        if read.is_secondary or read.is_supplementary:
            continue
        if reference_vntr.chromosome == read.reference_name:
            if start - read_length < read.reference_start < end:
                positive_reads.append(read)
                continue
        if read.qname in true_reads:
            false_negative_reads.append(read)
    return positive_reads, false_negative_reads


def find_info_by_mapping(sim_dir='simulation_data/', dir_index=0):
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {119: 'DRD4', 1220: 'GP1BA', 1221: 'CSTB', 1214: 'MAOA', 1219: 'IL1RN'}
    clean_up_tmp()
    dirs = glob.glob(sim_dir+'/*')
    simulation_dir = dirs[dir_index]
    files = glob.glob(simulation_dir + '/*')
    for fasta_file in files:
        if fasta_file.endswith('WGS_30x.fasta'):
            gene_name = simulation_dir.split('/')[-1].split('_')[0]
            vntr_id = None
            for vid, gname in id_to_gene.items():
                if gname == gene_name:
                    vntr_id = vid
            ref_vntr = reference_vntrs[vntr_id]

            true_reads_file = fasta_file[:-6] + '_true_reads.txt'
            if not os.path.exists(true_reads_file):
                simulated_sam_file = fasta_file[:-6] + '.sam'
                true_reads = get_id_of_reads_mapped_to_vntr_in_samfile(simulated_sam_file, ref_vntr)
                with open(true_reads_file) as out:
                    for true_read in true_reads:
                        out.write('%s\n' % true_read)
            else:
                with open(true_reads_file) as input:
                    lines = input.readlines()
                    true_reads = [line.strip() for line in lines if line.strip() != '']

            for i, parameter in enumerate([30]):
                positive_file = fasta_file[:-6] + '_bwa_%s_positive_reads.txt' % abs(parameter)
                false_negative_file = fasta_file[:-6] + '_bwa_%s_fp_reads.txt' % abs(parameter)
                if os.path.exists(positive_file) and os.path.exists(false_negative_file):
                    continue
                bwa_alignment_file = '/tmp/_gene%s_' % dir_index + 'bwa_alignment_%s.sam' % i
                bwa_alignment(fasta_file, bwa_alignment_file, parameter)
                positive_reads, fn_reads = get_positive_and_fn_reads_from_samfile(bwa_alignment_file, ref_vntr, true_reads)
                save_reads_stat(positive_file, positive_reads)
                save_reads_stat(false_negative_file, fn_reads)

                clean_up_tmp()

            for i, parameter in enumerate([-0.6]):
                positive_file = fasta_file[:-6] + '_bowtie_%s_positive_reads.txt' % abs(parameter)
                false_negative_file = fasta_file[:-6] + '_bowtie_%s_fp_reads.txt' % abs(parameter)
                if os.path.exists(positive_file) and os.path.exists(false_negative_file):
                    continue
                bowtie_alignment_file = '/tmp/_gene%s_' % dir_index + 'bowtie_alignment_%s.sam' % i
                bowtie_alignment(fasta_file, bowtie_alignment_file, parameter)
                positive_reads, fn_reads = get_positive_and_fn_reads_from_samfile(bowtie_alignment_file, ref_vntr, true_reads)
                save_reads_stat(positive_file, positive_reads)
                save_reads_stat(false_negative_file, fn_reads)

                clean_up_tmp()


find_info_by_mapping('simulation_data/', int(sys.argv[1]))
