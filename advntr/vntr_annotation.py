
ANNOTATION_DIR = 'results/annotation/'
EXONS = ANNOTATION_DIR + '%s_gene_exons.bed'
GENES = ANNOTATION_DIR + '%s_genes.bed'
ENSEMBL_TO_GENE = ANNOTATION_DIR + 'ensemblToGeneName.txt'
UCSC_TO_ENSMBL = ANNOTATION_DIR + 'knownToEnsembl.txt'
REFSEQ_TO_GENE = ANNOTATION_DIR + 'Refseq2Gene.txt'
PROMOTER_RANGE = 500


def intersect(s1, e1, s2, e2):
    return s1 <= e2 and s2 <= e1


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


def get_gene_name_from_refseq_id(query_refseq_id):
    with open(REFSEQ_TO_GENE) as infile:
        lines = infile.readlines()
        for line in lines:
            refseq_id, gene_name = line.strip().split()
            if refseq_id == query_refseq_id:
                return gene_name
    return 'None'


def get_gene_name_and_annotation_of_vntr(vntr_chromosome, vntr_start, vntr_end, genes, exons, gene_reference='refseq'):
    gene_name = 'None'
    annotation = 'Noncoding'

    for start, end, identifier in exons[vntr_chromosome]:
        if intersect(start, end, vntr_start, vntr_end):
            if gene_reference == 'ucsc':
                gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
            else:
                gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0])
            annotation = 'Coding'
            break
        if start > vntr_end:
            break

    if gene_name == 'None':
        for start, end, identifier in genes[vntr_chromosome]:
            if intersect(start, end, vntr_start, vntr_end):
                if gene_reference == 'ucsc':
                    gene_name = get_gene_name_from_ucsc_id(identifier.split('_')[0])
                else:
                    gene_name = get_gene_name_from_refseq_id(identifier.split('.')[0])
                annotation = 'Noncoding'
                break
            if start > vntr_end:
                break
    return gene_name, annotation


def is_vntr_close_to_gene(genes_info, vntr_chromosome, vntr_start, vntr_end):
    for start, end, _ in genes_info[vntr_chromosome]:
        if intersect(start, end, vntr_start, vntr_end):
            return True
        if start > vntr_end:
            return False
    return False


def get_exons_info(gene_reference='refseq'):
    exons_info = {}
    with open(EXONS % gene_reference) as infile:
        lines = infile.readlines()
        for line in lines:
            chromosome, start, end, identifier, _, _ = line.strip().split()
            start = int(start)
            end = int(end)
            if chromosome not in exons_info.keys():
                exons_info[chromosome] = []
            exons_info[chromosome].append((start, end, identifier))
    results = {}
    for chromosome, coordinates in exons_info.items():
        results[chromosome] = sorted(coordinates)
    return results


def get_genes_info(gene_reference='refseq'):
    genes_info = {}
    with open(GENES % gene_reference) as infile:
        genes_lines = infile.readlines()
        for line in genes_lines:
            line = line.strip().split()[:4]
            chromosome, start, end, identifier = line
            start = int(start)
            end = int(end)
            start -= PROMOTER_RANGE
            end += PROMOTER_RANGE
            if chromosome not in genes_info.keys():
                genes_info[chromosome] = []
            genes_info[chromosome].append((start, end, identifier))
    results = {}
    for chromosome, coordinates in genes_info.items():
        results[chromosome] = sorted(coordinates)
    return results
