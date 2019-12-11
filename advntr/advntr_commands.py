import logging
import os
import sys

from Bio import SeqIO

from advntr.genome_analyzer import GenomeAnalyzer
from advntr.models import load_unique_vntrs_data, get_largest_id_in_database, save_reference_vntr_to_database
from advntr.models import delete_vntr_from_database, create_vntrs_database
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


def get_default_vntrs(reference_vntrs, is_pacbio=False):
    pacbio_results = []
    illumina_results = []

    for ref_vntr in reference_vntrs:
        if not ref_vntr.is_non_overlapping() or ref_vntr.has_homologous_vntr():
            continue
        if 'N' in ref_vntr.left_flanking_region[-100:] or 'N' in ref_vntr.right_flanking_region[:100]:
            continue

        if ref_vntr.get_length() < 140 and ref_vntr.annotation in ['Coding', 'UTR', 'Promoter']:
            illumina = True
        else:
            illumina = False
        if ref_vntr.id in [532789, 188871, 301645, 468671, 503431]:
            illumina = True

        pacbio = True if illumina or ref_vntr.annotation in ['Coding', 'UTR', 'Promoter'] else False

        if ref_vntr.id in [3056, 25561, 69212, 415277, 519759, 379159, 532789, 70186, 188143, 193369, 193364, 258405,
                           188871, 301645, 400825, 468671]:
            pacbio = True

        if pacbio:
            pacbio_results.append(ref_vntr.id)
        if illumina:
            illumina_results.append(ref_vntr.id)

    if is_pacbio:
        return pacbio_results
    else:
        return illumina_results


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
    input_is_alignment_file = input_file.endswith('bam') or input_file.endswith('sam') or input_file.endswith('cram')
    if not input_is_alignment_file:
        print_error(genotype_parser, "The input file format is not supported. Please use BAM/CRAM files.")
    if args.working_directory is None:
        print_error(genotype_parser, 'Please specify working directory by -wd or --working_directory')
    working_directory = args.working_directory + '/' if args.working_directory else os.path.dirname(input_file) + '/'

    log_file = working_directory + 'log_%s.log' % os.path.basename(input_file)
    log_format = '%(asctime)s %(levelname)s:%(message)s'
    logging.basicConfig(format=log_format, filename=log_file, level=logging.DEBUG, filemode='w')

    if args.outfile:
        sys.stdout = open(args.outfile, 'w')

    models_file = args.models
    if models_file is None:
        models_file = settings.ILLUMINA_DEFAULT_MODELS_FILE
        if args.pacbio:
            models_file = settings.PACBIO_DEFAULT_MODELS_FILE
    settings.TRAINED_MODELS_DB = models_file
    settings.TRAINED_HMMS_DIR = os.path.dirname(os.path.realpath(settings.TRAINED_MODELS_DB)) + '/'

    reference_vntrs = load_unique_vntrs_data()
    target_vntrs = [ref_vntr.id for ref_vntr in reference_vntrs]
    if args.vntr_id is not None:
        target_vntrs = [int(vid) for vid in args.vntr_id.split(',')]
    logging.info('Running adVNTR for %s VNTRs' % len(target_vntrs))
    genome_analyzier = GenomeAnalyzer(reference_vntrs, target_vntrs, working_directory, args.outfmt, args.haploid,
                                      args.reference_filename, input_file)
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
                print_error(genotype_parser, '--frameshift is not available for these VNTRs')
        elif input_is_alignment_file:
            genome_analyzier.find_repeat_counts_from_alignment_file(input_file, average_coverage, args.update)
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
    if args.models is not None:
        settings.TRAINED_MODELS_DB = args.models

    genes = [gene.upper() for gene in args.gene.split(',') if gene]
    reference_vntrs = load_unique_vntrs_data()
    results = []
    for ref_vntr in reference_vntrs:
        if len(genes) and ref_vntr.gene_name.upper() not in genes:
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

    if args.models is not None:
        settings.TRAINED_MODELS_DB = args.models
    if not os.path.exists(settings.TRAINED_MODELS_DB):
        create_vntrs_database(settings.TRAINED_MODELS_DB)

    vntr_id = get_largest_id_in_database() + 1
    estimated_repeats = int((args.end - args.start) / len(args.pattern) + 5)
    ref_vntr = ReferenceVNTR(vntr_id, args.pattern, args.start, chromosome, None, None, estimated_repeats, chr_sequence)
    ref_vntr.init_from_vntrseek_data()
    vntr_finder = VNTRFinder(ref_vntr)

    print('Searching reference genome for regions with shared kmers with VNTR. It takes a few hours for human genome')
    scaled_recruitment_score = vntr_finder.train_classifier_threshold(args.reference)
    ref_vntr.scaled_score = scaled_recruitment_score
    save_reference_vntr_to_database(ref_vntr)
    print('Training completed. VNTR saved with ID: %s to the database' % vntr_id)


def del_model(args, delmodel_parser):
    if not args.vntr_id:
        print_error(delmodel_parser, '--vntr_id is required')
    if args.models is not None:
        settings.TRAINED_MODELS_DB = args.models
    delete_vntr_from_database(args.vntr_id)
