from subprocess import call
import shlex


def make_blast_database(fasta_files, db_name):
    for file_name in fasta_files:
        make_db_args = '-in %s -parse_seqids -dbtype nucl -title "%s"' % (file_name, file_name)
        make_db_args = shlex.split(make_db_args)
        call(['makeblastdb'] + make_db_args)

    file_names = ' '.join(fasta_files)
    merge_db_args = '-dblist %s -dbtype nucl -out %s -title "%s"' % (file_names, db_name, db_name)
    merge_db_args = shlex.split(merge_db_args)

    call(['blastdb_aliastool'] + merge_db_args)
