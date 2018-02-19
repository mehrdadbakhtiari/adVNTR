import logging
import os
import sys

from genome_analyzer import GenomeAnalyzer
from reference_vntr import load_unique_vntrs_data
import settings


def valid_vntr_for_frameshift(target_vntrs):
    for vntr_id in target_vntrs:
        if vntr_id not in settings.FRAMESHIFT_VNTRS:
            return False
    return True


def genotype(args, genotype_parser):
    if args.alignment_file is None and args.fasta is None:
        genotype_parser.print_help()
        sys.exit("ERROR: No input specified. Please specify alignment file or fasta file")

    if args.nanopore:
        settings.MAX_ERROR_RATE = 0.3
    elif args.pacbio:
        settings.MAX_ERROR_RATE = 0.3
    else:
        settings.MAX_ERROR_RATE = 0.05

    if args.threads < 1:
        genotype_parser.error('threads cannot be less than 1')
    settings.CORES = args.threads

    input_file = args.alignment_file if args.alignment_file else args.fasta
    input_is_alignment_file = input_file.endswith('bam') or input_file.endswith('sam')
    working_directory = args.working_directory + '/' if args.working_directory else os.path.dirname(input_file) + '/'

    log_file = working_directory + 'log_%s.log' % os.path.basename(input_file)
    log_format = '%(asctime)s %(levelname)s:%(message)s'
    logging.basicConfig(format=log_format, filename=log_file, level=logging.DEBUG, filemode='w')

    reference_vntrs = load_unique_vntrs_data()
    # reference_vntrs = identify_homologous_vntrs(reference_vntrs, 'chr15')
    illumina_targets = [1214, 1220, 1221, 1222, 1223, 1224, 377, 378, 809]

    target_vntrs = []
    for i in range(len(reference_vntrs)):
        if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
            continue
        target_vntrs.append(i)

    if args.vntr_id is not None:
        target_vntrs = [args.vntr_id]
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
            genome_analyzier.find_repeat_counts_from_alignment_file(input_file)
        else:
            genome_analyzier.find_repeat_counts_from_short_reads(input_file)


def not_implemented_command(parser):
    parser.error('Command has not been implemented yet. Sorry for inconvenience.')
