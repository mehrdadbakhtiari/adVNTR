import logging
import os
import sys

from Bio import SeqIO

from advntr.genome_analyzer import GenomeAnalyzer
from advntr.models import load_unique_vntrs_data, get_largest_id_in_database, save_reference_vntr_to_database
from advntr.reference_vntr import ReferenceVNTR
from advntr.vntr_finder import VNTRFinder
from advntr import settings


def valid_vntr_for_frameshift(target_vntrs):
    for vntr_id in target_vntrs:
        if vntr_id not in settings.FRAMESHIFT_VNTRS:
            return False
    return True


def print_error(subparser, msg):
    subparser.print_help()
    sys.exit("\nERROR: %s" % msg)


def genotype(args, genotype_parser):
    if args.alignment_file is None and args.fasta is None:
        print_error(genotype_parser, 'No input specified. Please specify alignment file or fasta file')

    if args.nanopore:
        settings.MAX_ERROR_RATE = 0.3
    elif args.pacbio:
        settings.MAX_ERROR_RATE = 0.3
    else:
        settings.MAX_ERROR_RATE = 0.05

    if args.threads < 1:
        print_error(genotype_parser, 'threads cannot be less than 1')
    settings.CORES = args.threads

    if args.expansion and args.coverage is None:
        print_error(genotype_parser, 'Please specify the average coverage to identify the expansion')
    average_coverage = args.coverage if args.expansion else None

    input_file = args.alignment_file if args.alignment_file else args.fasta
    input_is_alignment_file = input_file.endswith('bam') or input_file.endswith('sam')
    working_directory = args.working_directory + '/' if args.working_directory else os.path.dirname(input_file) + '/'

    settings.BLAST_TMP_DIR = working_directory + settings.BLAST_TMP_RELATIVE_DIR
    log_file = working_directory + 'log_%s.log' % os.path.basename(input_file)
    log_format = '%(asctime)s %(levelname)s:%(message)s'
    logging.basicConfig(format=log_format, filename=log_file, level=logging.DEBUG, filemode='w')

    settings.TRAINED_MODELS_DB = args.models
    settings.TRAINED_HMMS_DIR = os.path.dirname(os.path.realpath(settings.TRAINED_MODELS_DB)) + '/'
    reference_vntrs = load_unique_vntrs_data()
    # reference_vntrs = identify_homologous_vntrs(reference_vntrs, 'chr15')
    illumina_targets = [532789, 188871, 301645, 600000]

    target_vntrs = []
    for i in range(len(reference_vntrs)):
        if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
            continue
        target_vntrs.append(reference_vntrs[i].id)

    if args.vntr_id is not None:
        target_vntrs = [int(vid) for vid in args.vntr_id.split(',')]
    else:
        target_vntrs = illumina_targets
    genome_analyzier = GenomeAnalyzer(reference_vntrs, target_vntrs, working_directory)
    if args.pacbio:
        if input_is_alignment_file:
            genome_analyzier.find_repeat_counts_from_pacbio_alignment_file(input_file)
        else:
            genome_analyzier.find_repeat_counts_from_pacbio_reads(input_file, args.naive)
    else:
        if args.frameshift:
            if valid_vntr_for_frameshift(target_vntrs):
                genome_analyzier.find_frameshift_from_alignment_file(input_file)
            else:
                genotype_parser.error('--frameshift is only available for these IDs: %s' % settings.FRAMESHIFT_VNTRS)
        elif input_is_alignment_file:
            genome_analyzier.find_repeat_counts_from_alignment_file(input_file, average_coverage)
        else:
            genome_analyzier.find_repeat_counts_from_short_reads(input_file)


def print_models(reference_vntrs):
    print('VNTR ID\t| Chr\t| Gene\t| Start Position | Pattern')
    print('--------------------------------------------------')
    for ref_vntr in reference_vntrs:
        gene_name = ref_vntr.gene_name
        if len(gene_name) < 7:
            gene_name += '\t'
        print('%s\t| %s\t|%s| %s\t | %s' % (ref_vntr.id, ref_vntr.chromosome, gene_name,
                                            ref_vntr.start_point, ref_vntr.pattern))


def view_model(args, viewmodel_parser):
    valid_characters = {'A', 'C', 'G', 'T'}
    if args.pattern:
        for element in set(args.pattern.upper()):
            if element not in valid_characters:
                print_error(viewmodel_parser, 'Pattern should only contain A, C, G, T')

    genes = [gene.upper() for gene in args.gene.split(',') if gene]
    reference_vntrs = load_unique_vntrs_data()
    results = []
    for ref_vntr in reference_vntrs:
        if len(genes) and ref_vntr.gene_name not in genes:
            continue
        if args.pattern and ref_vntr.pattern != args.pattern.upper():
            continue
        # if ref_vntr.get_length() > 130:
        #     continue
        results.append(ref_vntr)
    print_models(results)


def add_model(args, addmodel_parser):
    if not args.reference:
        print_error(addmodel_parser, '--reference is required')
    if not args.chromosome:
        print_error(addmodel_parser, '--chromosome is required')
    if not args.pattern:
        print_error(addmodel_parser, '--pattern is required')
    if not args.start:
        print_error(addmodel_parser, '--start is required')
    if not args.end:
        print_error(addmodel_parser, '--end is required')

    chromosome = args.chromosome
    chr_sequence = ''

    fasta_sequences = SeqIO.parse(open(args.reference), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name == chromosome:
            chr_sequence = sequence
    vntr_id = get_largest_id_in_database() + 1
    estimated_repeats = (args.end - args.start) / len(args.pattern) + 5
    ref_vntr = ReferenceVNTR(vntr_id, args.pattern, args.start, chromosome, None, None, estimated_repeats, chr_sequence)
    ref_vntr.init_from_vntrseek_data()
    vntr_finder = VNTRFinder(ref_vntr)

    print('Searching reference genome for regions with shared kmers with VNTR. It takes a few hours for human genome')
    alphabet = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    m = 4194301

    def get_hash(string):
        result = 0
        for k in range(len(string)):
            result = (result + alphabet[string[k].upper()] * (4 ** (keyword_size - k - 1))) % m
        return result

    false_filtered_reads = []
    read_size = 150
    keyword_size = 11
    keywords = vntr_finder.get_keywords_for_filtering(True, keyword_size)
    hashed_keywords = set([get_hash(keyword) for keyword in keywords])
    match_positions = []
    fasta_sequences = SeqIO.parse(open(args.reference), 'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        if name != args.chromosome:
            continue
        window_hash = None
        for i in range(0, len(sequence) - keyword_size):
            if sequence[i].upper() == 'N' or sequence[i - 1 + keyword_size].upper() == 'N':
                continue
            if window_hash is None or sequence[i - 1].upper() == 'N':
                if 'N' in sequence[i:i + keyword_size].upper():
                    window_hash = None
                    continue
                window_hash = get_hash(sequence[i:i + keyword_size])
                continue
            window_hash -= alphabet[sequence[i - 1].upper()] * (4 ** (keyword_size - 1))
            window_hash = (window_hash * 4 + alphabet[sequence[i - 1 + keyword_size].upper()]) % m
            if window_hash in hashed_keywords:
                if name == args.chromosome and args.start - read_size < i < args.end:
                    continue
                if sequence[i:i + keyword_size].upper() in keywords:
                    match_positions.append(i)
                    if len(match_positions) > 3 and match_positions[-1] - match_positions[-3] < read_size:
                        for j in range(match_positions[-1] - read_size, match_positions[-3], 5):
                            false_filtered_reads.append(sequence[j:j + read_size])

    scaled_recruitment_score = vntr_finder.train_classifier_threshold(false_filtered_reads)
    ref_vntr.scaled_score = scaled_recruitment_score
    save_reference_vntr_to_database(ref_vntr)
    print('Training completed. VNTR saved with ID: %s to the database' % vntr_id)


def not_implemented_command(parser, command):
    parser.error('%s command has not been implemented yet. Sorry for inconvenience.' % command)
