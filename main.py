from distance import *


def read_fasta_seq(file_name='chr15.fa'):
    from Bio import SeqIO
    fasta_sequences = SeqIO.parse(open(file_name),'fasta')
    sequence = ''
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
    return sequence

sequence = ""
with open('chr15_mehrdad') as ifile:
    sequence = ifile.readline()

res['C'] = res['G'] = res['T'] = res['A'] = 0
found = 0
last_found = 0
for i in range(len(sequence)):
    if i > len(query):
        res[sequence[i-len(query)].upper()] -= 1
    if i % 200 == 0:
        found = 0
        last_found = 0
    res[sequence[i].upper()] += 1
    if nucleotide_dist(res, m) < 2:
        found = 1
        last_found = i