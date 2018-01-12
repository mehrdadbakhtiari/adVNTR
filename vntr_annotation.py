
ANNOTATION_DIR = 'annotation/'
EXONS = ANNOTATION_DIR + 'ucsc_gene_exons.bed'
GENES = ANNOTATION_DIR + 'ucsc_genes.bed'
ENSEMBL_TO_GENE = ANNOTATION_DIR + 'ensemblToGeneName.txt'
UCSC_TO_ENSMBL = ANNOTATION_DIR + 'knownToEnsembl.txt'
PROMOTER_RANGE = 300


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


def get_gene_name_and_annotation_of_vntr(vntr_chromosome, vntr_start, vntr_end):
    gene_name = 'None'
    annotation = 'Noncoding'
    with open(EXONS) as infile:
        lines = infile.readlines()
        for line in lines:
            chromosome, start, end, ucsc_id, _, _ = line.strip().split()
            start = int(start)
            end = int(end)
            if chromosome == vntr_chromosome and intersect(start, end, vntr_start, vntr_end):
                gene_name = get_gene_name_from_ucsc_id(ucsc_id.split('_')[0])
                annotation = 'Coding'
                break
    if gene_name == 'None':
        with open(GENES) as infile:
            lines = infile.readlines()
            for line in lines:
                line = line.strip().split()[:4]
                chromosome, start, end, ucsc_id = line
                start = int(start)
                end = int(end)
                start -= PROMOTER_RANGE
                end += PROMOTER_RANGE
                if chromosome == vntr_chromosome and intersect(start, end, vntr_start, vntr_end):
                    gene_name = get_gene_name_from_ucsc_id(ucsc_id.split('_')[0])
                    annotation = 'Noncoding'
                    break
    return gene_name, annotation
