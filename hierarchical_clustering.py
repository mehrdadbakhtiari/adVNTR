from Bio import SeqIO
from distance import hamming


def clusters_dist(c1, c2, distance_matrix):
    dist = 0
    for seq1 in c1:
        for seq2 in c2:
            dist += distance_matrix[seq1][seq2]
    return dist / (len(c1) * len(c2))


def find_closest_clusters(clusters, distance_matrix):
    min_dist = 1e10
    closest = (0, 0)
    for i in range(len(clusters)):
        for j in range(len(clusters)):
            if i == j:
                continue
            if clusters_dist(clusters[i], clusters[j], distance_matrix) < min_dist:
                min_dist = clusters_dist(clusters[i], clusters[j], distance_matrix)
                closest = (i, j)
    if closest[0] > closest[1]:
        closest = (closest[1], closest[0])
    return closest


def hierarchical_clustering(k, distance_matrix):
    clusters = [[i] for i in range(len(distance_matrix))]
    while len(clusters) > k:
        closest = find_closest_clusters(clusters, distance_matrix)
        new_cluster = clusters[closest[0]] + clusters[closest[1]]
        clusters.append(new_cluster)
        clusters = clusters[:closest[1]] + clusters[closest[1]+1:]
        clusters = clusters[:closest[0]] + clusters[closest[0]+1:]
    return clusters

fasta_sequences = SeqIO.parse(open('24149_VNTR1217_alignment_of_reads.fa'), 'fasta')
mat = []
seqs = []
for fasta in fasta_sequences:
    seqs.append(str(fasta.seq))

seqs = get_informative_columns(seqs)

for row in seqs:
    print(row)

for i in range(len(seqs)):
    mat.append([])
    for seq in seqs:
        mat[i].append(hamming(seq, seqs[i]))
# for line in mat:
#     print line
# print('')
clusters = hierarchical_clustering(2, mat)
print(clusters)
