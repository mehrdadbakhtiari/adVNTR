from multiprocessing import Process, Manager, Value, Semaphore
import os
import glob

def get_our_result(file_name):
    print('VeNTeR for %s' % file_name)
    copies = file_name.split('_')[-4]
    command = 'python main.py --alignment_file %s --frameshift --working_directory ./working_directory/ -t 20 > out_temp_%s.txt' % (file_name, copies)
    print(command)
    os.system(command)
    with open('out_temp_%s.txt' % copies) as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            if len(line) < 1:
                continue
            if line == 1123:
                continue
            if line == 'None':
                return False
            else:
                return True

def get_samtools_result(file_name):
    os.system('samtools mpileup -uf hg19_chromosomes/chr9.fa %s | bcftools view -bvcg - > var.raw.bcf' % (file_name))
    os.system('bcftools view var.raw.bcf | /usr/share/samtools/vcfutils.pl varFilter -D100 > var.flt.vcf')
    with open('var.flt.vcf') as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            if len(line) < 1:
                continue
            if line[0] == '#':
                continue
            print('True %s' % line)
            return True
    return False


def get_gatk_result(bowtie_bam_file):
    rg_bam_file = bowtie_bam_file[:-4] + '.rg.bam'
    copies = bowtie_bam_file.split('_')[-4]
    os.system('java -jar ~/picard.jar AddOrReplaceReadGroups I=%s O=%s RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20' % (bowtie_bam_file, rg_bam_file))
    os.system('samtools index %s %s' % (rg_bam_file, rg_bam_file[:-4] + '.bai'))
    os.system('java -jar ~/GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19_chromosomes/chr9.fa -I %s --genotyping_mode DISCOVERY -o gatk_30x_results/raw_variants_%s.vcf' % (rg_bam_file, copies))
    waiting_for_indel = False
    with open('gatk_30x_results/raw_variants_%s.vcf' % copies) as infile:
        lines = infile.readlines()
        for line in lines:
            line = line.strip()
            if len(line) < 1:
                continue
            if line.startswith('#CHROM'):
                waiting_for_indel = True
            if line[0] == '#':
                continue
            if waiting_for_indel:
                print('True %s' % line)
                return True
    return False


def get_freebayes_result(bowtie_bam_file):
    os.system('freebayes -f hg19_chromosomes/chr9.fa %s > %s' % (bowtie_bam_file, ))
    return A


def create_mapped_bam_file_and_index(base_name):
    os.system('bowtie2 -x data/CEL_frameshift_simulation/bowtie_reference/chr9_bowtie_index -1 %s1.fq -2 %s2.fq -S %s_bowtie_aln.sam --threads 20' % (base_name, base_name, base_name))
    os.system('samtools view -bS %s_bowtie_aln.sam | samtools sort - %s_bowtie_aln' % (base_name, base_name))
    os.system('samtools index %s_bowtie_aln.bam %s_bowtie_aln.bai' % (base_name, base_name))


files = glob.glob('data/CEL_frameshift_simulation/30x/*_30x.sam')
print(files)

ours = 0
result = False
sams = 0
sam_result = False
sams_correct = []
gatks = 0
for file_name in files:
    print(str(file_name))
    base_name = file_name[:-4]
#    create_mapped_bam_file_and_index(base_name)
    bowtie_bam_file = '%s_bowtie_aln.bam' % base_name
#    sam_result = get_samtools_result(bowtie_bam_file)
    if sam_result:
        sams += 1
        sams_correct.append(file_name)
#    result = get_our_result(bowtie_bam_file)
    if result:
        ours += 1
    gatk_result = get_gatk_result(bowtie_bam_file)
    if gatk_result:
        gatks += 1

#print(ours)
#print(sams)
print(gatks)
if len(sams_correct):
    print(sams_correct)

