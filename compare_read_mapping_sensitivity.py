import glob
import os
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

            if not os.path.exists(fasta_file[:-6] + '_true_reads.txt'):
                simulated_sam_file = fasta_file[:-6] + '.sam'
                true_reads = get_id_of_reads_mapped_to_vntr_in_samfile(simulated_sam_file, ref_vntr)
                with open(fasta_file[:-6] + '_true_reads.txt', 'w') as out:
                    for true_read in true_reads:
                        out.write('%s\n' % true_read)

            for i, parameter in enumerate([30, 20, 40]):
                result_file = fasta_file[:-6] + '_bwa_%s_mapped_reads.txt' % abs(parameter)
                if os.path.exists(result_file):
                    continue
                bwa_alignment_file = '/tmp/_gene%s_' % dir_index + 'bwa_alignment_%s.sam' % i
                bwa_alignment(fasta_file, bwa_alignment_file, parameter)
                bwa_ids = get_id_of_reads_mapped_to_vntr_in_samfile(bwa_alignment_file, ref_vntr)
                with open(result_file, 'w') as out:
                    for bwa_id in bwa_ids:
                        out.write('%s\n' % bwa_id)
                clean_up_tmp()

            for i, parameter in enumerate([-0.6, -0.3, -1]):
                result_file = fasta_file[:-6] + '_%s_mapped_reads.txt' % abs(parameter)
                if os.path.exists(result_file):
                    continue
                bowtie_alignment_file = '/tmp/_gene%s_' % dir_index + 'bowtie_alignment_%s.sam' % i
                bowtie_alignment(fasta_file, bowtie_alignment_file, parameter)
                bowtie_ids = get_id_of_reads_mapped_to_vntr_in_samfile(bowtie_alignment_file, ref_vntr)
                with open(result_file, 'w') as out:
                    for bowtie_id in bowtie_ids:
                        out.write('%s\n' % bowtie_id)
                clean_up_tmp()


find_info_by_mapping('simulation_data/', int(sys.argv[1]))
