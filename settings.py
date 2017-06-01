
HG19_DIR = './hg19_chromosomes/'
CHROMOSOMES = ['chr' + str(chr_number) for chr_number in list(range(1, 23))] + ['chrX', 'chrY']
# GENOME_LENGTH = 102531392
GENOME_LENGTH = 3100000000

USE_TRAINED_HMMS = True
TRAINED_HMMS_DIR = 'trained_HMMs/'

BLAST_TMP_DIR = 'blast_tmp/'

GC_CONTENT_WINDOW_SIZE = 100
OUTLIER_COVERAGE = 200

SEQTK_DIR = '../seqtk/seqtk'
