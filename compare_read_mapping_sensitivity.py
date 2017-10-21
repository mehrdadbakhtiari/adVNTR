import glob
import os
import sys

from sam_utils import get_id_of_reads_mapped_to_vntr_in_samfile


def clean_up_tmp():
    os.system('rm -rf /tmp/*.sam')
    os.system('rm -rf /tmp/*.fasta')


def bowtie_alignment(fasta_file, output, param):
    os.system('bowtie2 -x hg19_chromosomes/hg19_bt2_idx --end-to-end -f %s -S %s --threads 24 --score-min L,-0.6,%s' % (fasta_file, output, param))


def find_info_by_mapping(sim_dir='simulation_data/', dir_index=0):
    clean_up_tmp()
    dirs = glob.glob(sim_dir+'/*')
    simulation_dir = dirs[dir_index]
    files = glob.glob(simulation_dir + '/*')
    for fasta_file in files:
        if fasta_file.endswith('WGS_30x.fasta'):
            simulated_sam_file = fasta_file[:-6] + '.sam'
            true_reads = get_id_of_reads_mapped_to_vntr_in_samfile(simulated_sam_file)
            with open(fasta_file[:-6] + '_true_reads.txt', 'w') as out:
                for true_read in true_reads:
                    out.write('%s\n' % true_read)
            for i, parameter in enumerate([-0.6, -0.3, -1]):
                bowtie_alignment_file = '/mnt/_gene%s_' % dir_index + 'bowtie_alignment_%s.sam' % i
                bowtie_alignment(fasta_file, bowtie_alignment_file, parameter)
                bowtie_ids = get_id_of_reads_mapped_to_vntr_in_samfile(bowtie_alignment_file)
                with open(fasta_file[:-6] + '_%s_mapped_reads.txt' % i, 'w') as out:
                    for bowtie_id in bowtie_ids:
                        out.write('%s\n' % bowtie_id)
                clean_up_tmp()


find_info_by_mapping(int(sys.argv[1]))
