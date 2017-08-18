from StringIO import StringIO

from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO

from distance import hamming
from hierarchical_clustering import hierarchical_clustering
from settings import MUSCLE_DIR


class PacBioHaplotyper:
    def __init__(self, reads):
        self.reads = [read.upper() for read in reads]

    def get_error_corrected_haplotypes(self):
        haplotypes = []
        clusters = self.get_read_clusters()
        for cluster in clusters:
            muscle_cline = MuscleCommandline(MUSCLE_DIR, clwstrict=True)
            stdout, stderr = muscle_cline(stdin=data)
            data = '\n'.join(['>%s\n' % str(i) + cluster[i] for i in range(len(cluster))])
            alignment = AlignIO.read(StringIO(stdout), "clustal")
            aligned_reads = [str(aligned.seq) for aligned in alignment]
            seq = self.get_consensus_sequence_from_multiple_alignment(aligned_reads)
            haplotypes.append(seq)
        return haplotypes

    @staticmethod
    def get_consensus_sequence_from_multiple_alignment(aligned_reads):
        """
        It finds the most frequent element in each column using counting sort.
        Then it add the element to the result if it is not a gap element.
        """
        seq = ''
        for i in range(len(aligned_reads[0])):
            bins = {}
            for row in aligned_reads:
                if row[i] in bins.keys():
                    bins[row[i]] += 1
                else:
                    bins[row[i]] = 0
            sorted_frequencies = sorted(bins.items(), key=lambda x: x[1])
            most_frequent_element = sorted_frequencies[-1][0]
            if most_frequent_element != '-':
                seq += most_frequent_element

        return seq

    def get_read_clusters(self):
        muscle_cline = MuscleCommandline(MUSCLE_DIR, clwstrict=True)
        data = '\n'.join(['>%s\n' % str(i) + self.reads[i] for i in range(len(self.reads))])
        stdout, stderr = muscle_cline(stdin=data)
        alignment = AlignIO.read(StringIO(stdout), "clustal")
        aligned_reads = [str(aligned.seq) for aligned in alignment]
        aligned_read_ids = [str(aligned.id) for aligned in alignment]

        seqs = self.get_informative_columns(aligned_reads)
        distance_matrix = []
        for i in range(len(seqs)):
            distance_matrix.append([])
            for seq in seqs:
                distance_matrix[i].append(hamming(seq, seqs[i]))
        clusters = hierarchical_clustering(2, distance_matrix)

        result = [[self.reads[int(aligned_read_ids[i])] for i in cluster] for cluster in clusters]
        return result

    def get_informative_columns(self, aligned_reads):
        result = ['' for _ in aligned_reads]
        for col in range(len(aligned_reads[0])-1):
            current = '$'
            counter = 0
            for row in aligned_reads:
                if row[col] == current:
                    counter += 1
                else:
                    counter -= 1
                if counter < 0:
                    counter = 0
                    current = row[col]
            if counter <= 3:
                for i in range(len(aligned_reads)):
                    result[i] += aligned_reads[i][col]
        return result
