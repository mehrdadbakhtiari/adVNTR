from subprocess import call
import shlex
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import Seq, SeqRecord, SeqIO
import settings


def make_blast_database(fasta_files, db_name):
    for file_name in fasta_files:
        make_db_args = '-in %s -parse_seqids -dbtype nucl -title "%s"' % (file_name, file_name)
        make_db_args = shlex.split(make_db_args)
        call(['makeblastdb'] + make_db_args)

    file_names = ' '.join(fasta_files)
    merge_db_args = '-dblist "%s" -dbtype nucl -out %s -title "%s"' % (file_names, db_name, db_name)
    merge_db_args = shlex.split(merge_db_args)

    call(['blastdb_aliastool'] + merge_db_args)


def get_blast_matched_ids(query, blast_db_name, word_size=None, max_seq='6000', evalue=10.0, search_id=''):
    query_file = settings.BLAST_TMP_DIR + search_id + '_query.fasta'
    result_file = settings.BLAST_TMP_DIR + search_id + '_blast_result.txt'
    with open(query_file, "w") as output_handle:
        my_rec = SeqRecord.SeqRecord(seq=Seq.Seq(query), id='query', description='')
        SeqIO.write([my_rec], output_handle, 'fasta')

    if word_size is None:
        if len(query) < 30:
            word_size = '7'
        else:
            word_size = '11'
    if len(query) <= 15:
        task = 'blastn-short'
    else:
        task = 'blastn'
    blastn_cline = NcbiblastnCommandline(query=query_file, db=blast_db_name, outfmt='"6 sallseqid"', dust='no',
                                         out=result_file, num_threads="4", word_size=word_size, max_target_seqs=max_seq,
                                         evalue=evalue, task=task)
    blastn_cline()

    with open(result_file) as result_input:
        ids = result_input.readlines()
        matched_ids = set([seq_id.strip() for seq_id in ids])
    os.remove(result_file) if os.path.exists(result_file) else None
    os.remove(query_file) if os.path.exists(query_file) else None

    return matched_ids
