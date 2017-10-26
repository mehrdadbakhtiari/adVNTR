import argparse
import logging
import os
import sys

from genome_analyzer import GenomeAnalyzer
from reference_vntr import identify_homologous_vntrs, load_unique_vntrs_data
import settings
# from vntr_graph import plot_graph_components, get_nodes_and_edges_of_vntr_graph

parser = argparse.ArgumentParser(description='VNTRFinder 1.0.0')
parser.add_argument('-a', '--alignment_file', type=str, help='Alignment file in BAM format or SAM format',
                    metavar='FILE')
parser.add_argument('-f', '--fasta', type=str, help='Fasta file containing raw reads', metavar='FILE')
parser.add_argument('-fs', '--frameshift', action='store_true',
                    help='Search for a frameshift in VNTR instead of copy number')
parser.add_argument('-p', '--pacbio', action='store_true',
                    help='Input file contains PacBio reads instead of Illumina reads')
parser.add_argument('-n', '--nanopore', action='store_true',
                    help='Input file contains Nanopore MinION reads instead of Illumina reads')
parser.add_argument('-wd', '--working_directory', type=str, metavar='DIRECTORY',
                    help='Working directory for creating temporary files needed for computation')
parser.add_argument('-t', '--threads', type=int, metavar='<nthreads>', default=4,
                    help='Run the tool on <nthreads> parallel threads which will run on separate processors/cores')
parser.add_argument('-vid', '--vntr_id', type=int, metavar='<VNTR ID>', default=None,
                    help='ID of the VNTR')
parser.add_argument('-naive', '--naive', action='store_true', default=False,
                    help='Use naive approach for PacBio reads')
args = parser.parse_args()

if args.alignment_file is None and args.fasta is None:
    parser.print_help()
    sys.exit("ERROR: No input specified. Please specify alignment file or fasta file")

if args.nanopore:
    settings.MAX_ERROR_RATE = 0.3
elif args.pacbio:
    settings.MAX_ERROR_RATE = 0.2
else:
    settings.MAX_ERROR_RATE = 0.05

if args.threads < 1:
    parser.error('threads cannot be less than 1')
settings.CORES = args.threads

input_file = args.alignment_file if args.alignment_file else args.fasta
input_is_alignment_file = input_file.endswith('bam') or input_file.endswith('sam')
working_directory = args.working_directory if args.working_directory else os.path.dirname(input_file) + '/'

LOGFILE = working_directory + 'log_%s.log' % os.path.basename(input_file)
logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', filename=LOGFILE, level=logging.DEBUG, filemode='w')

reference_vntrs = load_unique_vntrs_data()

# reference_vntrs = identify_homologous_vntrs(reference_vntrs, 'chr15')
accurate_vntr_list = [7, 69, 119, 970, 1123, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 809, 377, 378]
accurate_vntr_list += [7, 119, 1214, 1218, 1220, 1221, 377, 378, 809] # short VNTRs
accurate_vntr_list += [1123, 1214, 1220, 1221, 1222] # grant
#accurate_vntr_list = [377, 378, 809, 69, 1123] # frameshift
#accurate_vntr_list = [69, 1123]
accurate_vntr_list += [1215] # INS
accurate_vntr_list += [970, 1213, 1215, 1216, 1217, 1219, 1222]

target_vntrs = []
for i in range(len(reference_vntrs)):
    if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
        continue
    if reference_vntrs[i].id in accurate_vntr_list:
        continue
    target_vntrs.append(i)

if args.vntr_id:
    target_vntrs = [args.vntr_id]
genome_analyzier = GenomeAnalyzer(reference_vntrs, target_vntrs, working_directory)
if args.pacbio:
    if input_is_alignment_file:
        copy_number = genome_analyzier.find_repeat_counts_from_pacbio_alignment_file(input_file)
    else:
        copy_number = genome_analyzier.find_repeat_counts_from_pacbio_reads(input_file, args.naive)
else:
    if args.frameshift:
        genome_analyzier.find_frameshift_from_alignment_file(input_file)
    elif input_is_alignment_file:
        copy_number = genome_analyzier.find_repeat_counts_from_alignment_file(input_file)
    else:
        copy_number = genome_analyzier.find_repeat_counts_from_short_reads(input_file)

# print(len(reference_vntrs))
# nodes, edges = get_nodes_and_edges_of_vntr_graph(reference_vntrs)
# plot_graph_components(nodes, edges)
