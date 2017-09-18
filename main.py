import argparse
import logging
import os
import sys

from vntr_finder import VNTRFinder
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
parser.add_argument('-wd', '--working_directory', type=str, metavar='DIRECTORY',
                    help='Working directory for creating temporary files needed for computation')
parser.add_argument('-t', '--threads', type=int, metavar='<nthreads>', default=1,
                    help='Run the tool on <nthreads> parallel threads which will run on separate processors/cores')
args = parser.parse_args()

if args.alignment_file is None and args.fasta is None:
    parser.print_help()
    sys.exit("ERROR: No input specified. Please specify alignment file or fasta file")

if args.threads < 1:
    parser.error('threads cannot be less than 1')
settings.CORES = args.threads

input_file = args.alignment_file if args.alignment_file else args.fasta
input_is_alignment_file = input_file.endswith('bam') or input_file.endswith('sam')
working_directory = args.working_directory if args.working_directory else os.path.dirname(input_file) + '/'

LOGFILE = 'log_%s.log' % os.path.basename(input_file)
logging.basicConfig(format='%(asctime)s %(levelname)s:%(message)s', filename=LOGFILE, level=logging.DEBUG, filemode='w')

reference_vntrs = load_unique_vntrs_data()

# reference_vntrs = identify_homologous_vntrs(reference_vntrs, 'chr15')
accurate_vntr_list = [7, 69, 119, 970, 1123, 1213, 1214, 1215, 1216, 1217, 1218, 1219, 1220, 809, 377, 378]
accurate_vntr_list = [7, 119, 1214, 1218, 1220, 1221, 377, 378, 809]

for i in range(len(reference_vntrs)):
    if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
        continue
    if reference_vntrs[i].id not in accurate_vntr_list:
        continue
    print(i, len(reference_vntrs[i].get_repeat_segments()))
    vntr_finder = VNTRFinder(reference_vntrs[i])
    if args.pacbio:
        if input_is_alignment_file:
            copy_number = vntr_finder.find_repeat_count_from_pacbio_alignment_file(input_file, working_directory)
        else:
            copy_number = vntr_finder.find_repeat_count_from_pacbio_reads(input_file, working_directory)
    else:
        if input_is_alignment_file:
            copy_number = vntr_finder.find_repeat_count_from_alignment_file(input_file, working_directory)
        else:
            copy_number = vntr_finder.find_repeat_count_from_short_reads(input_file)

# print(len(reference_vntrs))
# nodes, edges = get_nodes_and_edges_of_vntr_graph(reference_vntrs)
# plot_graph_components(nodes, edges)
