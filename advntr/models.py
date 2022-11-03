import os
import sqlite3

from multiprocessing import Process, Semaphore, Manager

from Bio import Seq, SeqRecord, SeqIO

from advntr.reference_vntr import ReferenceVNTR

from advntr.vntr_annotation import get_gene_name_and_annotation_of_vntr, is_vntr_close_to_gene, get_genes_info
from advntr import settings
from advntr.utils import get_chromosome_reference_sequence

from threading import Lock
lock = Lock()


PROCESSING_DIR = 'whole_genome_loc/'


def load_unprocessed_vntrseek_data(vntrseek_output, chromosome=None, only_genic=False):
    vntrs = []
    if only_genic:
        genes_info = get_genes_info()
    else:
        genes_info = {}
    chr_seq = get_chromosome_reference_sequence(chromosome)
    with open(vntrseek_output) as input_file:
        input_lines = [line.strip() for line in input_file.readlines() if line.strip() != '']
        for vntr_id, line in enumerate(input_lines):
            vntrseek_repeat, _, pattern, chromosome_number, start = line.split()
            if len(pattern) > 100 or len(pattern) < 6:
                continue
            start = int(start) - 1
            estimated_repeats = int(float(vntrseek_repeat) + 2)
            if chromosome is not None and chromosome_number != chromosome:
                continue
            estimated_end = estimated_repeats * len(pattern) + start
            if only_genic and not is_vntr_close_to_gene(genes_info, chromosome_number, start, estimated_end):
                continue
            vntrs.append(ReferenceVNTR(vntr_id, pattern, start, chromosome_number, None, None, estimated_repeats, chromosome_sequence=chr_seq))
    print('%s VNTRs are close to a gene' % len(vntrs))
    return vntrs


def find_non_overlapping_vntrs(manager, vntrs, result, index, chrom=None, output_file=None, sema=None):
    skipped_vntrs = set([])
    for i in range(len(vntrs)):
        if chrom is not None and vntrs[i].chromosome != chrom:
            continue
        estimated_end = len(vntrs[i].pattern) * vntrs[i].estimated_repeats + vntrs[i].start_point
        print(i, estimated_end - vntrs[i].start_point)
        if i < len(vntrs)-1 and vntrs[i].chromosome == vntrs[i+1].chromosome and estimated_end > vntrs[i+1].start_point:
            vntrs[i].estimated_repeats += vntrs[i+1].estimated_repeats
        if len(vntrs[i].pattern) * vntrs[i].estimated_repeats > 1000:
            vntrs[i].non_overlapping = False
            continue
        vntrs[i].init_from_vntrseek_data()
        repeat_segments = vntrs[i].get_repeat_segments()
        if i in skipped_vntrs:
            vntrs[i].non_overlapping = False
        else:
            j = i + 1
            end_point = len(vntrs[i].pattern) * len(repeat_segments) + vntrs[i].start_point
            while j < len(vntrs) and vntrs[i].chromosome == vntrs[j].chromosome and end_point > vntrs[j].start_point:
                skipped_vntrs.add(j)
                j += 1
    with lock:
        print('writing %s for %s' % (len(vntrs), chrom))
        for vntr in vntrs:
            if not vntr.is_non_overlapping():
                continue
            repeat_segments = ','.join(vntr.get_repeat_segments())
            with open(output_file, 'a') as out:
                end_point = vntr.start_point + vntr.get_length()
                gene_name, annotation = None, None
                out.write('%s %s %s %s %s %s %s %s %s %s\n' % (vntr.id, vntr.is_non_overlapping(), vntr.chromosome,
                                                               vntr.start_point, gene_name, annotation, vntr.pattern,
                                                               vntr.left_flanking_region, vntr.right_flanking_region,
                                                               repeat_segments,))
    res = manager.list(result[index])
    result[index] = res
    vntrs = None
    if sema is not None:
        sema.release()


def process_vntrseek_data(unprocessed_vntrs_file, output_file='vntr_data/VNTRs.txt', chrom=None):
    process_list = []
    unprocessed_vntrs = load_unprocessed_vntrseek_data(unprocessed_vntrs_file, chrom)
    sema = Semaphore(settings.CORES)
    manager = Manager()
    partial_vntrs = manager.list([])
    for i in range(settings.CORES):
        sema.acquire()
        partial_vntrs.append(manager.list())
        q = len(unprocessed_vntrs) / settings.CORES
        start = i * q
        end = (i+1) * q if i + 1 < settings.CORES else len(unprocessed_vntrs)
        partial_input = unprocessed_vntrs[start:end]
        p = Process(target=find_non_overlapping_vntrs, args=(manager, partial_input, partial_vntrs, i, chrom, output_file, sema))
        process_list.append(p)
        p.start()
    for p in process_list:
        p.join()
    print(chrom, 'Done')


def identify_homologous_vntrs(vntrs, chromosome=None):
    for i in range(len(vntrs)):
        for j in range(i + 1, len(vntrs)):
            if chromosome and (chromosome != vntrs[i].chromosome and chromosome != vntrs[j].chromosome):
                continue
            if vntrs[i].is_homologous_vntr(vntrs[j]):
                vntrs[i].has_homologous = True
                vntrs[j].has_homologous = True
    return vntrs


def create_vntrs_database(db_file):
    if not os.path.exists(os.path.dirname(db_file)):
        os.makedirs(os.path.dirname(db_file))
    db = sqlite3.connect(db_file)
    cursor = db.cursor()
    cursor.execute('''
    CREATE TABLE vntrs(id INTEGER PRIMARY KEY, nonoverlapping TEXT, chromosome TEXT, ref_start INTEGER, gene_name TEXT,
    annotation TEXT, pattern TEXT, left_flanking TEXT, right_flanking TEXT, repeats TEXT, scaled_score REAL default 0)
    ''')

    db.commit()
    db.close()


def load_unique_vntrs_data(db_file=None):
    vntrs = []
    if db_file is None:
        db_file = settings.TRAINED_MODELS_DB
    db = sqlite3.connect(db_file)
    cursor = db.cursor()
    cursor.execute('''SELECT id, nonoverlapping, chromosome, ref_start, gene_name, annotation, pattern, left_flanking,
    right_flanking, repeats, scaled_score FROM vntrs''')

    for row in cursor:
        new_row = []
        for element in row:
            if type(element) != int and type(element) != float:
                new_row.append(str(element))
            else:
                new_row.append(element)
        vntr_id, overlap, chrom, start, gene, annotation, pattern, left_flank, right_flank, segments, score = new_row
        if "," in segments:
            repeat_segments = segments.split(',')
        else:
            repeat_segments = []
        repeats = len(repeat_segments)
        vntr = ReferenceVNTR(int(vntr_id), pattern, int(start), chrom, gene, annotation, repeats, scaled_score=score)
        vntr.init_from_xml(repeat_segments, left_flank, right_flank)
        vntr.non_overlapping = True if overlap == 'True' else False
        vntrs.append(vntr)

    return vntrs


def save_vntrs_to_database(processed_vntrs, db_file):
    with open(processed_vntrs) as input_file:
        lines = input_file.readlines()
    db = sqlite3.connect(db_file)
    cursor = db.cursor()
    singles = 0
    paired = 0
    for line in lines:
        line = line.strip()
        vntr_id, overlap, chromosome, start, gene, annotation, pattern, left_flank, right_flank, segments = line.split()
        cursor.execute('''INSERT INTO vntrs(id, nonoverlapping, chromosome, ref_start, gene_name, annotation, pattern,
                       left_flanking, right_flanking, repeats, scaled_score) VALUES(?,?,?,?,?,?,?,?,?,?,?)''',
                       (vntr_id, overlap, chromosome, start, gene, annotation, pattern, left_flank, right_flank,
                        segments, 0))
        if len(segments) < 145:
            singles += 1
        if len(segments) < 290:
            paired += 1
    db.commit()
    db.close()
    print('%s %s %s' % (processed_vntrs, singles, paired))


def update_trained_score_in_database(vntr_id, scaled_recruitment_score):
    db = sqlite3.connect(settings.TRAINED_MODELS_DB)
    cursor = db.cursor()
    cursor.execute('''UPDATE vntrs SET scaled_score=? WHERE id=?''', (scaled_recruitment_score, vntr_id))
    db.commit()
    db.close()


def update_gene_name_and_annotation_in_database(vntr_id, gene_name, annotation, db_file=None):
    if db_file is None:
        db_file = settings.TRAINED_MODELS_DB
    db = sqlite3.connect(db_file)
    cursor = db.cursor()
    cursor.execute('''UPDATE vntrs SET gene_name=?, annotation=? WHERE id=?''', (gene_name, annotation, vntr_id))
    db.commit()
    db.close()


def save_reference_vntr_to_database(ref_vntr, db_file=None):
    if db_file is None:
        db_file = settings.TRAINED_MODELS_DB
    db = sqlite3.connect(db_file)
    cursor = db.cursor()
    segments = ','.join(ref_vntr.get_repeat_segments())
    non_overlapping = "True" if ref_vntr.non_overlapping else "False"
    cursor.execute('''INSERT INTO vntrs(id, nonoverlapping, chromosome, ref_start, gene_name, annotation, pattern,
                   left_flanking, right_flanking, repeats, scaled_score) VALUES(?,?,?,?,?,?,?,?,?,?,?)''',
                   (ref_vntr.id, non_overlapping, ref_vntr.chromosome, ref_vntr.start_point,
                    ref_vntr.gene_name, ref_vntr.annotation, ref_vntr.pattern, ref_vntr.left_flanking_region,
                    ref_vntr.right_flanking_region, segments, ref_vntr.scaled_score))
    db.commit()
    db.close()


def get_largest_id_in_database():
    db = sqlite3.connect(settings.TRAINED_MODELS_DB)
    cursor = db.cursor()
    cursor.execute('''SELECT MAX(id) FROM vntrs''')
    result = 0
    for row in cursor:
        if row[0] is not None:
            result = row[0]
    return result


def delete_vntr_from_database(vntr_id):
    db = sqlite3.connect(settings.TRAINED_MODELS_DB)
    cursor = db.cursor()
    cursor.execute('DELETE FROM vntrs WHERE id=%s' % vntr_id)
    db.commit()
    db.close()


def is_false_vntr_hit(qresult, ref_vntr):
    for hit in qresult:
        for hsp in hit:
            if int(hsp.hit_id) == ref_vntr.id:
                continue
            score = hsp.match_num - hsp.mismatch_num - hsp.hit_gapopen_num
            length = len(ref_vntr.pattern) + 60
            if score / float(length) > 0.75:
                print(hsp.hit_id, hsp.hit_start)
                return True
    return False


def find_similar_region_for_vntr(sema, reference_vntr, ref_file, result_list):
    from Bio import SearchIO
    vntr_id = reference_vntr.id
    q = reference_vntr.left_flanking_region[-30:] + reference_vntr.pattern + reference_vntr.right_flanking_region[:30]
    search_index = vntr_id
    qfile = str(vntr_id) + '_' + str(search_index) + '_query.fasta'
    with open(qfile, "w") as output_handle:
        my_rec = SeqRecord.SeqRecord(seq=Seq.Seq(q), id='query', description='')
        SeqIO.write([my_rec], output_handle, 'fasta')
    output = 'blat_out/output_%s_%s.psl' % (vntr_id, search_index)
    command = 'blat -q=dna -oneOff=1 -tileSize=8 -stepSize=3 -minIdentity=75 %s %s %s' % (ref_file, qfile, output)
    os.system(command)
    os.system('rm %s' % qfile)
    try:
        qresult = SearchIO.read(output, 'blat-psl')
        if is_false_vntr_hit(qresult, reference_vntr):
            print('there is similar sequence for %s' % vntr_id)
            result_list.append(vntr_id)
    except ValueError:
        pass
    sema.release()


def identify_similar_regions_for_vntrs_using_blat():
    from multiprocessing import Process, Semaphore, Manager
    reference_vntrs = load_unique_vntrs_data()

    records = []
    for ref_vntr in reference_vntrs:
        record = SeqRecord.SeqRecord('')
        sequence = ref_vntr.left_flanking_region[-30:] + ref_vntr.pattern + ref_vntr.right_flanking_region[:30]
        record.seq = Seq.Seq(sequence)
        record.id = str(ref_vntr.id)
        records.append(record)
    vntr_structures_file = 'reference_vntr_structures.fa'
    with open(vntr_structures_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

    sema = Semaphore(7)
    manager = Manager()
    result_list = manager.list()
    process_list = []
    for ref_vntr in reference_vntrs:
        sema.acquire()
        p = Process(target=find_similar_region_for_vntr, args=(sema, ref_vntr, vntr_structures_file, result_list))
        process_list.append(p)
        p.start()

    for p in process_list:
        p.join()
    result_list = list(result_list)
    with open('similar_vntrs.txt', 'a') as out:
        for vntr_id in result_list:
            out.write('%s\n' % vntr_id)


def extend_flanking_regions_in_processed_vntrs(flanking_size=500, output_file='vntr_data/repeats_and_segments2.txt'):
    vntrs = load_unique_vntrs_data()
    reference_genomes = {}
    for vntr in vntrs:
        comma_separated_segments = ','.join(vntr.get_repeat_segments())
        if vntr.chromosome not in reference_genomes.keys():
            reference_genomes[vntr.chromosome] = get_chromosome_reference_sequence(vntr.chromosome)
        start = vntr.start_point
        left_flanking_region = reference_genomes[vntr.chromosome][start-flanking_size:start].upper()
        end = vntr.start_point + vntr.get_length()
        right_flanking_region = reference_genomes[vntr.chromosome][end:end+flanking_size].upper()
        with open(output_file, 'a') as out:
            out.write('%s %s %s %s %s\n' % (vntr.id, vntr.is_non_overlapping(), left_flanking_region,
                                            right_flanking_region, comma_separated_segments))


def fill_vntr_database():
    for chrom in settings.CHROMOSOMES:
        processed_vntrs = PROCESSING_DIR + 'VNTRs_%s.txt' % chrom
        database_file = PROCESSING_DIR + 'GRCh38_VNTRs_%s.db' % chrom
        create_vntrs_database(database_file)
        save_vntrs_to_database(processed_vntrs, database_file)

def test_load_save_vntrs_basic():
    from pprint import pprint

    db_file = os.path.join(os.getcwd(), "test_db_file.db")
    vntrs = []
    vntr_0 = ReferenceVNTR(vntr_id=0, pattern="ACC", start_point=10,
                         chromosome="chr1", gene_name="test_gene",
                         annotation="promoter", estimated_repeats=0,
                         scaled_score=0.2)
    vntrs.append(vntr_0)
    vntr_1 = ReferenceVNTR(vntr_id=1, pattern="ACCTTTGG", start_point=5,
                         chromosome="chr22", gene_name="test_gene",
                         annotation="promoter", estimated_repeats=0,
                         scaled_score=0.1)
    vntr_1.non_overlapping = False
    vntrs.append(vntr_1)
    vntr_2 = ReferenceVNTR(vntr_id=3, pattern="CA", start_point=300,
                         chromosome="chrX", gene_name="test_gene",
                         annotation="enhancer", estimated_repeats=2,
                         scaled_score=0.6)
    vntr_2.repeat_segments = ["CA", "CAC"]
    vntrs.append(vntr_2)
    vntr_3 = ReferenceVNTR(vntr_id=4, pattern="CA", start_point=105,
                         chromosome="chrX", gene_name="test_gene",
                         annotation="enhancer", estimated_repeats=0,
                         scaled_score=0.6)

    vntr_3.left_flanking_region = "CAAA"
    vntr_3.right_flanking_region = "CCCC"
    vntrs.append(vntr_3)
    # Empty the db_file and create a fresh database.
    open(db_file, "w").close()
    create_vntrs_database(db_file=db_file)
    # Save the test vntrs on the db_file.
    for vntr in vntrs:
        save_reference_vntr_to_database(ref_vntr=vntr, db_file=db_file)
    # Load the vntrs from the db_file.
    loaded_vntrs = load_unique_vntrs_data(db_file=db_file)

    # Compare the list of vntrs, not caring about the order.
    assert(len(vntrs) == len(loaded_vntrs))
    for i in range(len(vntrs)):
        try:
            assert(loaded_vntrs[i] == vntrs[i])
        except AssertionError:
            print("Assertion failed for index {}".format(i))
            print("loaded_vntrs[{}]: ".format(i))
            pprint(vars(loaded_vntrs[i]))
            print("vntrs[{}]: ".format(i))
            pprint(vars(vntrs[i]))

if __name__ == "__main__":
    import sys
    for i, chrom in enumerate(settings.CHROMOSOMES):
        if i != int(sys.argv[1]):
            continue
        sorted_vntrseek_output = PROCESSING_DIR + '/repeats_length_patterns_chromosomes_starts.txt'
        processed_vntrs = PROCESSING_DIR + '/VNTRs_%s.txt' % chrom
        process_vntrseek_data(sorted_vntrseek_output, processed_vntrs, chrom)
        break
