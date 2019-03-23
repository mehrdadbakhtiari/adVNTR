import os
import sqlite3

from multiprocessing import Process, Semaphore, Manager

from Bio import Seq, SeqRecord, SeqIO

from advntr.reference_vntr import ReferenceVNTR
from advntr.vntr_annotation import get_gene_name_and_annotation_of_vntr, is_vntr_close_to_gene, get_genes_info
from advntr import settings
from advntr.utils import get_chromosome_reference_sequence


def load_unprocessed_vntrseek_data(vntrseek_output, chromosome=None):
    vntrs = []
    genes_info = get_genes_info()
    with open(vntrseek_output) as input_file:
        input_lines = [line.strip() for line in input_file.readlines() if line.strip() != '']
        for vntr_id, line in enumerate(input_lines):
            vntrseek_repeat, _, pattern, chromosome_number, start = line.split()
            if len(pattern) > 100:
                continue
            start = int(start) - 1
            estimated_repeats = int(float(vntrseek_repeat) + 5)
            if chromosome is not None and chromosome_number != chromosome:
                continue
            estimated_end = estimated_repeats * len(pattern) + start
            if not is_vntr_close_to_gene(genes_info, chromosome_number, start, estimated_end):
                continue
            vntrs.append(ReferenceVNTR(vntr_id, pattern, start, chromosome_number, None, None, estimated_repeats))
    print('%s VNTRs are close to a gene' % len(vntrs))
    return vntrs


def find_non_overlapping_vntrs(vntrs, result, chrom=None, sema=None):
    skipped_vntrs = []
    for i in range(len(vntrs)):
        if chrom is not None and vntrs[i].chromosome != chrom:
            continue
        estimated_end = len(vntrs[i].pattern) * vntrs[i].estimated_repeats + vntrs[i].start_point
        print(i, estimated_end - vntrs[i].start_point)
        if i < len(vntrs)-1 and vntrs[i].chromosome == vntrs[i+1].chromosome and estimated_end > vntrs[i+1].start_point:
            vntrs[i].estimated_repeats += vntrs[i+1].estimated_repeats
        vntrs[i].init_from_vntrseek_data()
        repeat_segments = vntrs[i].get_repeat_segments()
        if i in skipped_vntrs:
            vntrs[i].non_overlapping = False
        else:
            j = i + 1
            end_point = len(vntrs[i].pattern) * len(repeat_segments) + vntrs[i].start_point
            while j < len(vntrs) and vntrs[i].chromosome == vntrs[j].chromosome and end_point > vntrs[j].start_point:
                skipped_vntrs.append(j)
                j += 1
    for vntr in vntrs:
        result.append(vntr)
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
        p = Process(target=find_non_overlapping_vntrs, args=(partial_input, partial_vntrs[i], chrom, sema))
        process_list.append(p)
        p.start()
    for p in process_list:
        p.join()
    vntrs = []
    for partial_list in partial_vntrs:
        vntrs.extend(list(partial_list))
    print(chrom, len(vntrs))

    for vntr in vntrs:
        if not vntr.is_non_overlapping():
            continue
        repeat_segments = ','.join(vntr.get_repeat_segments())
        with open(output_file, 'a') as out:
            end_point = vntr.start_point + vntr.get_length()
            gene_name, annotation = get_gene_name_and_annotation_of_vntr(vntr.chromosome, vntr.start_point, end_point)
            out.write('%s %s %s %s %s %s %s %s %s %s\n' % (vntr.id, vntr.is_non_overlapping(), vntr.chromosome,
                                                           vntr.start_point, gene_name, annotation, vntr.pattern,
                                                           vntr.left_flanking_region, vntr.right_flanking_region,
                                                           repeat_segments,))


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
    if not os.path.exists(os.path.basename(db_file)):
        os.makedirs(os.path.basename(db_file))
    db = sqlite3.connect(db_file)
    cursor = db.cursor()
    cursor.execute('''
    CREATE TABLE vntrs(id INTEGER PRIMARY KEY, nonoverlapping TEXT, chromosome TEXT, ref_start INTEGER, gene_name TEXT,
    annotation TEXT, pattern TEXT, left_flanking TEXT, right_flanking TEXT, repeats TEXT, scaled_score REAL default 0)
    ''')

    db.commit()
    db.close()


def load_unique_vntrs_data():
    vntrs = []
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
        repeat_segments = segments.split(',')
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


def update_gene_name_and_annotation_in_database(vntr_id, gene_name, annotation):
    db = sqlite3.connect(settings.TRAINED_MODELS_DB)
    cursor = db.cursor()
    cursor.execute('''UPDATE vntrs SET gene_name=?, annotation=? WHERE id=?''', (gene_name, annotation, vntr_id))
    db.commit()
    db.close()


def save_reference_vntr_to_database(ref_vntr):
    db = sqlite3.connect(settings.TRAINED_MODELS_DB)
    cursor = db.cursor()
    segments = ','.join(ref_vntr.get_repeat_segments())
    cursor.execute('''INSERT INTO vntrs(id, nonoverlapping, chromosome, ref_start, gene_name, annotation, pattern,
                   left_flanking, right_flanking, repeats, scaled_score) VALUES(?,?,?,?,?,?,?,?,?,?,?)''',
                   (ref_vntr.id, ref_vntr.non_overlapping, ref_vntr.chromosome, ref_vntr.start_point,
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
        processed_vntrs = 'vntr_data/VNTRs_%s.txt' % chrom
        database_file = 'vntr_data/hg19_VNTRs.db'
        save_vntrs_to_database(processed_vntrs, database_file)

if __name__ == "__main__":
    pass
    # process_vntrseek_data(sorted_vntrseek_output, processed_vntrs, chrom)
