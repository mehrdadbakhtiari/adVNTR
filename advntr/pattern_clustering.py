
from Bio import pairwise2
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import numpy.matlib as matlib


def get_sequence_distance(s, t):
    max_length = max(len(s), len(t))
    return max_length - pairwise2.align.globalxx(s, t, score_only=True)


def get_distance_matrix(patterns):
    distance_matrix = []
    for p1 in patterns:
        dist_row = []
        for p2 in patterns:
            dist_row.append(get_sequence_distance(p1, p2))
        distance_matrix.append(dist_row)
    return distance_matrix


def get_cluster_similarities(clusters, distance_matrix):
    result = []
    for cluster in clusters:
        sim = 0.0
        for item1 in cluster:
            for item2 in cluster:
                sim += distance_matrix[item1][item2]
        sim /= len(cluster) ** 2
        result.append(sim)
    return result


def get_elbow_point_index(wcss):
    curve = wcss
    number_of_points = len(curve)
    all_coordinates = np.vstack((range(number_of_points), curve)).T
    np.array([range(number_of_points), curve])
    first_point = all_coordinates[0]
    line_vector = all_coordinates[-1] - all_coordinates[0]
    line_vector_norm = line_vector / np.sqrt(np.sum(line_vector**2))
    vec_from_first = all_coordinates - first_point
    scalar_product = np.sum(vec_from_first * matlib.repmat(line_vector_norm, number_of_points, 1), axis=1)
    vec_from_first_parallel = np.outer(scalar_product, line_vector_norm)
    vectors_to_line = vec_from_first - vec_from_first_parallel
    dists_to_line = np.sqrt(np.sum(vectors_to_line ** 2, axis=1))
    index_of_best_point = np.argmax(dists_to_line)
    return index_of_best_point


def get_pattern_clusters(patterns):
    if len(patterns) == 1:
        return patterns
    distance_matrix = get_distance_matrix(patterns)
    distortions = []
    clusterings = []
    for k in range(1, len(patterns)+1):
        f = AgglomerativeClustering(affinity='precomputed', linkage='complete', n_clusters=k).fit(distance_matrix)
        clusters = [[] for _ in range(k)]
        for pattern_index, label in enumerate(f.labels_):
            clusters[label].append(pattern_index)
        cluster_similarities = get_cluster_similarities(clusters, distance_matrix)
        clustering_quality = sum(cluster_similarities) / float(len(cluster_similarities))

        distortions.append(clustering_quality)
        clusterings.append(clusters)

    distortions.reverse()
    clusterings.reverse()

    best_clustering_index = get_elbow_point_index(distortions)
    best_clustering = clusterings[best_clustering_index]
    result = []
    for cluster in best_clustering:
        result.append([patterns[i] for i in cluster])
    return result
