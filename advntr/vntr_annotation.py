
ANNOTATION_DIR = 'results/hg38_annotation/'
EXONS = ANNOTATION_DIR + '%s_gene_coding_exons.bed'
INTRONS = ANNOTATION_DIR + '%s_gene_introns.bed'
UTR5 = ANNOTATION_DIR + '%s_gene_5utr.bed'
UTR3 = ANNOTATION_DIR + '%s_gene_3utr.bed'
GENES = ANNOTATION_DIR + '%s_genes.bed'
ENSEMBL_TO_GENE = ANNOTATION_DIR + 'ensemblToGeneName.txt'
UCSC_TO_ENSMBL = ANNOTATION_DIR + 'knownToEnsembl.txt'
REFSEQ_TO_GENE = ANNOTATION_DIR + 'Refseq2Gene.txt'
PROMOTER_RANGE = 500


def intersect(s1, e1, s2, e2):
    return s1 <= e2 and s2 <= e1

def include(s1, e1, vntr_s, vntr_e):
    return s1 <= vntr_s <= vntr_e <= e1

def get_gene_name_from_ensmbl(known_ensmbl):
    with open(ENSEMBL_TO_GENE) as infile:
        lines = infile.readlines()
        for line in lines:
            ensmbl, gene_name = line.strip().split()
            if ensmbl == known_ensmbl:
                return gene_name
    return 'None'


def get_gene_name_from_ucsc_id(known_ucsc_id):
    with open(UCSC_TO_ENSMBL) as infile:
        lines = infile.readlines()
        for line in lines:
            ucsc_id, ensmbl = line.strip().split()
            if ucsc_id == known_ucsc_id:
                return get_gene_name_from_ensmbl(ensmbl)
    return 'None'


def get_refseq_id_to_gene_name_map():
    result = {}
    with open(REFSEQ_TO_GENE) as infile:
        lines = infile.readlines()
        for line in lines:
            refseq_id, gene_name = line.strip().split()
            result[refseq_id] = gene_name
    return result


def get_gene_name_from_refseq_id(query_refseq_id, mapping):
    if query_refseq_id in mapping:
        return mapping[query_refseq_id]
    return 'None'


def get_gene_name_and_annotation_of_vntr(vntr_chromosome, vntr_start, vntr_end, genes, exons, introns, utr3, utr5, name_mapping=None, gene_reference='refseq'):
    if name_mapping is None:
        name_mapping = get_refseq_id_to_gene_name_map()

    def get_annotation(vntr_start, vntr_end, regions, region_name='Coding'):
        gene_name, annotation = 'None', 'None'
        for start, end, identifier, direction in regions:
            if intersect(start, end, vntr_start, vntr_end):
                if gene_reference == 'ucsc':
                    gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
                else:
                    gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0], name_mapping)
                annotation = region_name
                break
            if start > vntr_end:
                break
        return gene_name, annotation

    gene_name, annotation = get_annotation(vntr_start, vntr_end, exons[vntr_chromosome])

    if gene_name == 'None':
        gene_name, annotation = get_annotation(vntr_start, vntr_end, utr5[vntr_chromosome], 'UTR')
    if gene_name == 'None':
        gene_name, annotation = get_annotation(vntr_start, vntr_end, utr3[vntr_chromosome], 'UTR')
    if gene_name == 'None':
        gene_name, annotation = get_annotation(vntr_start, vntr_end, introns[vntr_chromosome], 'Intron')

    if gene_name == 'None':
        for start, end, identifier, direction in genes[vntr_chromosome]:
            if direction == '+':
                start = end
                end += PROMOTER_RANGE
            else:
                end = start
                start -= PROMOTER_RANGE
            if intersect(start, end, vntr_start, vntr_end):
                if gene_reference == 'ucsc':
                    gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
                else:
                    gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0], name_mapping)
                annotation = 'Promoter'
                break
            if start - PROMOTER_RANGE > vntr_end:
                break
    return gene_name, annotation


def is_vntr_close_to_gene(genes_info, vntr_chromosome, vntr_start, vntr_end):
    for start, end, _, _direction in genes_info[vntr_chromosome]:
        if intersect(start-1000, end+1000, vntr_start, vntr_end):
            return True
        if start > vntr_end:
            return False
    return False


def get_translate_ranges(exons, gene_reference='refseq'):
    translate_ranges = {}
    mapping = get_refseq_id_to_gene_name_map()

    for chromosome in exons.keys():
        for start, end, identifier, direction in exons[chromosome]:
            if gene_reference == 'ucsc':
                gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
            else:
                gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0], mapping)
            if gene_name not in translate_ranges.keys():
                translate_ranges[gene_name] = (start, end)
            else:
                current_start, current_end = translate_ranges[gene_name]
                translate_ranges[gene_name] = (min(start, current_start), max(end, current_end))
    return translate_ranges


def get_exons_info(annotation_file=EXONS, gene_reference='refseq'):
    exons_info = {}
    number_of_segments = {}
    with open(annotation_file % gene_reference) as infile:
        lines = infile.readlines()
        for line in lines:
            chromosome, start, end, identifier, _, direction = line.strip().split()
            start = int(start)
            end = int(end)
            if chromosome not in exons_info.keys():
                exons_info[chromosome] = []
            segment_number = int(identifier.split('.')[1].split('_')[2])
            exons_info[chromosome].append((start, end, identifier, direction, segment_number))
            number_of_segments[identifier.split('.')[0]] = int(segment_number)
    results = {}
    for chromosome, coordinates in exons_info.items():
        results[chromosome] = sorted(coordinates)
    return results, number_of_segments

def is_within_coding_exon(vntr_chromosome, vntr_start, vntr_end, exons):
    for start, end, _, _, _ in exons[vntr_chromosome]:
        if start > vntr_end:
            break
        if include(start, end, vntr_start, vntr_end):
            return True
    return False

def intersects_with_coding_exon(vntr_chromosme, vntr_start, vntr_end, exons):
    for start, end, _, _, _ in exons[vntr_chromosome]:
        if start > vntr_end:
            break
        if intersect(start, end, vntr_start, vntr_end):
            return True
    return False

def get_RepeatMasker_info(repeat_ref_file):
    # Repeating Elements by RepeatMasker
    # http://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=rep&hgta_track=rmsk&hgta_table=rmsk&hgta_doSchema=describe+table+schema

    from collections import defaultdict
    repeat_info = defaultdict(list)
    with open(repeat_ref_file, "r") as infile:
        for line in infile:
            split_line = line.strip().split()
            # Columns
            # bin   swScore	milliDiv	milliDel	milliIns
            # genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily
            # repStart	repEnd	repLeft	id

            # repName: (TAACCC)n, TAR1, L1MC5a, MER5B, (TGG)n, L3...
            # repClass: Simple_repeat, Satellite, LINE, DNA, SINE, ...
            # repFamily: Simple_repeat, telo, L1, MIR,
            _, _, _, _, _, chromosome, genoStart, genoEnd, _, strand, repName, repClass, repFamily, _, _, _, _ = split_line
            repeat_info[chromosome].append((int(genoStart), int(genoEnd), strand, repName, repClass, repFamily))

    results = {}
    for chromosome, coordinates in repeat_info.items():
        results[chromosome] = sorted(coordinates)
    return results

def is_within_alu_elements(vntr_chromosome, vntr_start, vntr_end, repeat_masker_info):
    for start, end, strand, repName, repClass, repFamily, in repeat_masker_info[vntr_chromosome]:
        if repClass == "LINE" or repClass =="SINE":
            if start > vntr_end:
                break
            if include(start, end, vntr_start, vntr_end):
                return True
    return False

def get_genes_info(gene_reference='refseq'):
    genes_info = {}
    with open(GENES % gene_reference) as infile:
        genes_lines = infile.readlines()
        for line in genes_lines:
            line = line.strip().split()[:6]
            chromosome, start, end, identifier, _, direction = line
            start = int(start)
            end = int(end)
            if chromosome not in genes_info.keys():
                genes_info[chromosome] = []
            genes_info[chromosome].append((start, end, identifier, direction))
    results = {}
    for chromosome, coordinates in genes_info.items():
        results[chromosome] = sorted(coordinates)
    return results


def sort_file(filename):
    with open(filename) as infile:
        lines = infile.readlines()
    lines = [line.split() for line in lines]
    for i in range(len(lines)):
        lines[i][1] = int(lines[i][1])
    lines.sort(key=lambda tup: tup[1])
    lines.sort(key=lambda tup: tup[0])
    with open(filename, 'w') as outfile:
        for line in lines:
            for i in range(len(line)):
                outfile.write('%s\t' % line[i])
            outfile.write('\n')


def get_introns_count(introns_info):
    introns_count = {}
    found_genes = set()
    for c in introns_info.keys():
        for start, end, id, dir in introns_info[c]:
            id = id.split('.')[0]
            if id not in found_genes:
                found_genes.add(id)
                introns_count[id] = 0
            introns_count[id] += 1
    return introns_count

def get_introns(identifier, chromosome, introns_count=None):
    return introns_count[identifier]

def get_intron_count(vntr_start, vntr_end, chromosome, regions):
    index = 0
    res = None
    for start, end, identifier, direction in regions[chromosome]:
        if intersect(start, end, vntr_start, vntr_end):
            if gene_reference == 'ucsc':
                gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
            else:
                gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0], name_mapping)
            if direction == '+':
                res = index + 1
            else:
                res = get_introns(identifier.split('.')[0], chromosome) - index
            break
        if start > vntr_end:
            break
    return res

if __name__ == '__main__':
    from models import update_gene_name_and_annotation_in_database, load_unique_vntrs_data
    genes_info = get_genes_info()
    exons_info = get_exons_info()
    introns_info = get_exons_info(INTRONS)
    utr5_info = get_exons_info(UTR5)
    utr3_info = get_exons_info(UTR3)
    name_mapping = get_refseq_id_to_gene_name_map()
    # translate_ranges = get_translate_ranges(exons_info)
    db_file = 'vntr_data/hg38_genic_VNTRs.db'
    reference_vntrs = load_unique_vntrs_data(db_file)
    for ref_vntr in reference_vntrs:
        end = ref_vntr.start_point + ref_vntr.get_length()
        new_gene, new_annotation = get_gene_name_and_annotation_of_vntr(ref_vntr.chromosome, ref_vntr.start_point, end, genes_info, exons_info, introns_info, utr3_info, utr5_info, name_mapping)
        if new_gene == 'None' and ref_vntr.gene_name != 'None':
            new_gene = ref_vntr.gene_name
        print(ref_vntr.id)
        update_gene_name_and_annotation_in_database(ref_vntr.id, new_gene, new_annotation, db_file)
