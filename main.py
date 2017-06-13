from vntr_finder import VNTRFinder
from reference_vntr import identify_homologous_vntrs, load_unique_vntrs_data
# from vntr_graph import plot_graph_components, get_nodes_and_edges_of_vntr_graph


read_files = ['original_reads/paired_dat1.fasta', 'original_reads/paired_dat2.fasta']
alignment_file = '12878_reads_1/CEUTrio.HiSeq.WGS.b37_decoy.NA12878.clean.dedup.recal.20120117.bam'
reference_vntrs = load_unique_vntrs_data()

# reference_vntrs = identify_homologous_vntrs(reference_vntrs, 'chr15')
accurate_vntr_list = [271, 281, 283, 287, 288, 325, 327, 328, 329]

for i in range(len(reference_vntrs)):
    if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
        continue
    if reference_vntrs[i].id not in accurate_vntr_list:
        continue
    print(i, len(reference_vntrs[i].get_repeat_segments()))
    vntr_finder = VNTRFinder(reference_vntrs[i])
    # copy_number = vntr_finder.find_repeat_count_from_short_reads(read_files)
    copy_number = vntr_finder.find_repeat_count_from_alignment_file(alignment_file)
    # vntr_finder.find_accuracy()

    with open('hmm_repeat_count.txt', 'a') as output:
        output.write('%s %s\n' % (i, copy_number / len(reference_vntrs[i].get_repeat_segments())))
    # end_point = start_points[i] + sum([len(e) for e in repeat_segments])
    # VNTR_coverage_ratio = get_VNTR_coverage_over_total_coverage(start_points[i], end_point)
    # with open('vntr_coverage_ratio.txt', 'a') as output:
    #     output.write('%s %s\n' % (i, VNTR_coverage_ratio))

# print(len(reference_vntrs))
# nodes, edges = get_nodes_and_edges_of_vntr_graph(reference_vntrs)
# plot_graph_components(nodes, edges)
