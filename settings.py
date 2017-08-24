import logging

HG19_DIR = './hg19_chromosomes/'
CHROMOSOMES = ['chr' + str(chr_number) for chr_number in list(range(1, 23))] + ['chrX', 'chrY']
GENOME_LENGTH = 3100000000
MAX_INSERT_SIZE = 1000

USE_TRAINED_HMMS = True
TRAINED_HMMS_DIR = 'trained_HMMs/'
SCORE_FINDING_READS_FRACTION = 0.0001

BLAST_TMP_DIR = 'blast_tmp/'

GC_CONTENT_WINDOW_SIZE = 100
GC_CONTENT_BINS = 10
OUTLIER_COVERAGE = 200

CORES = 96

SEQTK_DIR = '../seqtk/seqtk'
MUSCLE_DIR = 'tools/muscle3.8.31_i86linux64'

LOG_FILE = 'log.log'
logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', filename=LOG_FILE, level=logging.DEBUG)
