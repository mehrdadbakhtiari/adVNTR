from Bio import Seq, SeqRecord, SeqIO

from advntr.utils import get_chromosome_reference_sequence
from advntr.models import load_unique_vntrs_data


def add_two_copy_to_all_patterns(patterns, start_points):
    file_name = 'chr15.fa'
    fasta_sequences = SeqIO.parse(open(file_name), 'fasta')
    sequence = ''
    record = SeqRecord.SeqRecord('')
    for fasta in fasta_sequences:
        record = fasta
        name, sequence = fasta.id, str(fasta.seq)
    total_added_chars = 0
    for i in range(len(patterns)):
        start_point = start_points[i] + total_added_chars
        pattern = patterns[i]
        sequence = sequence[:start_point] + pattern * 2 + sequence[start_point:]
        total_added_chars += len(pattern * 2)

    record.seq = Seq.Seq(sequence)
    output_name = 'edited_chr15_two_more_copies.fa'
    with open(output_name, 'w') as output_handle:
        SeqIO.write([record], output_handle, 'fasta')


def create_reference_with_indel(ref_vntr, output_name, position, insertion=True, inserted_bp='C'):
    sequence = get_chromosome_reference_sequence(ref_vntr.chromosome)
    left_flank = sequence[ref_vntr.start_point-1000:ref_vntr.start_point]
    vntr_end = ref_vntr.start_point + ref_vntr.get_length()
    vntr = sequence[ref_vntr.start_point:vntr_end]
    right_flank = sequence[vntr_end:vntr_end+1000]

    if insertion:
        sequence = left_flank + vntr[:position] + inserted_bp + vntr[position:] + right_flank
    else:
        sequence = left_flank + vntr[:position] + vntr[position+1:] + right_flank
    record = SeqRecord.SeqRecord('')
    record.seq = Seq.Seq(sequence)
    with open(output_name, 'w') as output_handle:
        SeqIO.write([record], output_handle, 'fasta')


def create_cel_frameshifts(cel_vntr):
    sequence = get_chromosome_reference_sequence(cel_vntr.chromosome)
    left_flank = sequence[cel_vntr.start_point-3000:cel_vntr.start_point]
    vntr_end = cel_vntr.start_point + cel_vntr.get_length()
    vntr = sequence[cel_vntr.start_point:vntr_end]
    right_flank = sequence[vntr_end:vntr_end+3000]

    for frameshift in xrange(160, 360, 10):
        deletion = left_flank + vntr[:frameshift] + vntr[frameshift+1:] + right_flank
        insertion = left_flank + vntr[:frameshift+5] + 'C' + vntr[frameshift+5:] + right_flank
        record = SeqRecord.SeqRecord('')
        record.seq = Seq.Seq(deletion)
        deletion_output_name = 'chr9_cel_deletion_%s.fa' % str(frameshift)
        insertion_output_name = 'chr9_cel_insertion_%s.fa' % str(frameshift + 5)
        with open(deletion_output_name, 'w') as output_handle:
            SeqIO.write([record], output_handle, 'fasta')
        record.seq = Seq.Seq(insertion)
        with open(insertion_output_name, 'w') as output_handle:
            SeqIO.write([record], output_handle, 'fasta')


def create_reference_region_with_specific_repeats(reference_vntr, desired_repeats_count, output_name, flanks=30000, repeat_patterns=None):
    record = SeqRecord.SeqRecord('')
    sequence = get_chromosome_reference_sequence(reference_vntr.chromosome)
    vntr_end = reference_vntr.start_point + reference_vntr.get_length()
    if flanks is None:
        region_start = 0
        region_end = len(sequence)
    else:
        region_start = reference_vntr.start_point - flanks
        region_end = vntr_end + flanks
    new_sequence = sequence[region_start:reference_vntr.start_point]
    if repeat_patterns is None:
        repeats = reference_vntr.get_repeat_segments()
    else:
        repeats = repeat_patterns
    for i in range(desired_repeats_count):
        new_sequence += repeats[i % len(repeats)]
    new_sequence += sequence[vntr_end:region_end]

    record.seq = Seq.Seq(new_sequence)
    with open(output_name, 'w') as output_handle:
        SeqIO.write([record], output_handle, 'fasta')


def create_illumina_genotyping_references(illumina_read_dir='../Genotyping/'):
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {1220: 'GP1BA', 1221: 'CSTB', 1214: 'MAOA'}
    repeats = {'GP1BA': range(1, 5), 'CSTB': range(1, 16), 'MAOA': range(1, 6)}
    repeats_patterns = {'GP1BA': ['AGCCCGACCACCCCAGAGCCCACCTCAGAGCCCGCCCCC', 'AGCCCGACCACCCCGGAGCCCACCTCAGAGCCCGCCCCC', 'AGCCCGACCACCCCGGAGCCCACCCCAATCCCGACCATCGCCA'],
                        'CSTB': ['CGCGGGGCGGGG', 'CGCGGGGCGGGG', 'CGCGGGGCGGGG', 'CGGCGGGCGGGG'],
                        'MAOA': ['ACCGGCACCGGCACCAGTACCCGCACCAGT', 'ACCGGCACCGGCACCGAGCGCAAGGCGGAG', 'ACCGGCACCGGCACCAGTACCCGCACCAGT']}

    for vntr_id in id_to_gene.keys():
        # if vntr_id != 1221:
        #     continue
        for repeat in repeats[id_to_gene[vntr_id]]:
            outfile = illumina_read_dir + id_to_gene[vntr_id] + '/' + str(repeat) + '.pacfa'
            create_reference_region_with_specific_repeats(reference_vntrs[vntr_id], repeat, outfile, 3000, repeats_patterns[id_to_gene[vntr_id]])


def create_illumina_copy_number_variation_references(illumina_read_dir='../Illumina_copy_number/'):
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {119: 'DRD4', 1220: 'GP1BA', 1221: 'CSTB', 1214: 'MAOA', 1219: 'IL1RN'}
    repeats = {'DRD4': range(1, 12), 'GP1BA': range(1, 6), 'CSTB': range(1, 16), 'MAOA': range(1, 6),
               'IL1RN': range(1, 10)}

    for vntr_id in id_to_gene.keys():
        for repeat in repeats[id_to_gene[vntr_id]]:
            outfile = illumina_read_dir + id_to_gene[vntr_id] + '/' + str(repeat) + '.fa'
            create_reference_region_with_specific_repeats(reference_vntrs[vntr_id], repeat, outfile, 149)


def create_pacbio_copy_number_variation_references(pacbio_read_dir='../pacbio_recruitment/set1/'):
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {1221: 'CSTB', 1216: 'HIC1', 1215: 'INS'}
    repeats = {'CSTB': range(1, 69), 'HIC1': range(2, 36), 'INS': range(10, 171)}

    for vntr_id in id_to_gene.keys():
        for repeat in repeats[id_to_gene[vntr_id]]:
            if id_to_gene[vntr_id] == 'INS' and repeat % 5 != 0:
                continue
            if id_to_gene[vntr_id] == 'CSTB' and repeat % 2 != 0:
                continue
            outfile = pacbio_read_dir + id_to_gene[vntr_id] + '/' + str(repeat) + '.fa'
            create_reference_region_with_specific_repeats(reference_vntrs[vntr_id], repeat, outfile, 1000)


def create_pacbio_coverage_data_for_3_genes_and_10_cn(pacbio_read_dir='../pacbio_coverage_experiment/'):
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {1221: 'CSTB', 1216: 'HIC1', 1215: 'INS'}
    repeats = {'CSTB': range(2, 42), 'HIC1': range(2, 22), 'INS': range(10, 110)}

    for vntr_id in id_to_gene.keys():
        for repeat in repeats[id_to_gene[vntr_id]]:
            if id_to_gene[vntr_id] == 'INS' and repeat % 5 != 0:
                continue
            if id_to_gene[vntr_id] == 'CSTB' and repeat % 2 != 0:
                continue
            if id_to_gene[vntr_id] != 'INS':
                continue
            outfile = pacbio_read_dir + id_to_gene[vntr_id] + '/' + str(repeat) + '.fa'
            create_reference_region_with_specific_repeats(reference_vntrs[vntr_id], repeat, outfile, 3000)


def create_pacbio_ru_length_data_for_all_vntrs(pacbio_read_dir='../pacbio_ru_data_for_all_vntrs/'):
    reference_vntrs = load_unique_vntrs_data()

    with open('vntr_complex.txt') as infile:
        lines = infile.readlines()
        complex_vntrs = [int(r.strip().split()[0]) for r in lines] + [0]

    repeat_units = {}
    for vntr_id in range(len(reference_vntrs)):
        if vntr_id in complex_vntrs:
            continue
        ru = len(reference_vntrs[vntr_id].pattern)
        if ru not in repeat_units.keys():
            repeat_units[ru] = []
        if len(repeat_units[ru]) >= 4:
            continue
        repeat_units[ru].append(vntr_id)

    import os
    for ru in repeat_units.keys():
        if len(repeat_units[ru]) < 2:
            continue
        for vntr_id in repeat_units[ru]:
            original_repeats = len(reference_vntrs[vntr_id].get_repeat_segments())
            start = max(3, original_repeats - 10)
            for repeat in range(start, start + 21):
                if repeat % 5 != 0:
                    continue
                outfile = pacbio_read_dir + str(ru) + '/vntr_id_' + str(vntr_id) + '_' + str(repeat) + '.fa'
                if not os.path.exists(os.path.dirname(outfile)):
                    os.makedirs(os.path.dirname(outfile))
                create_reference_region_with_specific_repeats(reference_vntrs[vntr_id], repeat, outfile, 1000)

if __name__ == '__main__':
    create_pacbio_ru_length_data_for_all_vntrs()
