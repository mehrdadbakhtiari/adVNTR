def hamming(s1, s2):
    result = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            result += 1
    return result


def get_similarity(s1, s2):
    similarity = 0
    for i in range(min(len(s1), len(s2))):
        if s1[i] == s2[i]:
            similarity += 1
    return similarity


def get_nucleotide_map(sequence):
    n_map = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    for s in sequence:
        n_map[s] += 1
    return n_map


def nucleotide_dist(map1, map2):
    nucleotides = 'ACTG'
    res = 0
    for n in nucleotides:
        res += abs(map1[n] - map2[n])
    return res
