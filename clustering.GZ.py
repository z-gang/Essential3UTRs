## TODO: Rewrite using numpy
import numpy as np

import sys


import serialization


chrom_names = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII', 'XIV', 'XV', 'XVI', 'Mito']


def count_samples(X, i, eps, sample_weight):
    left = right = i

    while right < len(X) - 1:
        if X[right + 1] <= X[i] + eps:
            right += 1
        else:
            break

    while left > 0:
        if X[left - 1] >= X[i] + eps:
            left -= 1
        else:
            break

    return sum(sample_weight[left: right + 1]), left, right

"""
## TODO: Implement proper DBSCAN
def linear_dbscan(X, eps=0.5, min_samples=5, metric='eucledean', sample_weight=None):
    if not eps > 0.0:
        raise ValueError("eps must be positive.")

    if sample_weight is None:
        sample_weight = [1 for x in X]

    labels = [-1 for x in X]

    i = 0
    last_label = -1
    cur_label = -1

    while i < len(X):
        samples_num, left, right = count_samples(X, i, eps, sample_weight)

        if samples_num < min_samples:
            if cur_label == -1:
                i += 1
                continue
            else:
                raise Exception('Wrong order of samples 1')

        if samples_num >= min_sample:
            if cur_label == -1:
                j = i + 1
                next_counter, next_left, next_rigth = count_samples(X, j, eps, sample_weight)
                prev_counter, prev_left, prev_right = samples_num, left, right


                if right == j:
                    cluster


                while X[j] - X[j - 1] <= eps and next_counter >= min_samples:
                    prev_counter, prev_left, prev_right = next_counter, next_left, next_right
                    j += 1
                    next_counter, next_left, next_rigth = count_samples(X, j, eps, sample_weight)

                if X[j] - X[j - 1] > eps:


                    for k in range(left, next_right + 1):
                        if labels[k] != -1:
                            raise ValueError('sample has been already labeled as a different cluster')
                        labels[k] = last_label + 1



            else:
                raise Exception('Wrong order of samples 2')
"""
 
## TODO: Add ability to use different metrics
def find_dbscan_core_points(X, eps=1, min_samples=5, metric='eucledean', sample_weight=None):
    if not eps > 0.0:
        raise ValueError("eps must be positive.")

    if sample_weight is None:
        sample_weight = [1 for x in X]

    labels = [-1 for x in X]

    i = 0
    last_label = -1

    while i < len(X):
        samples_num, left, right = count_samples(X, i, eps, sample_weight)

        if samples_num < min_samples:
            i += 1
            continue

        if right == i:
            i += 1
            continue

        prev_counter, prev_left, prev_right = samples_num, left, right
        j = i + 1
        next_counter, next_left, next_right = count_samples(X, j, eps, sample_weight)

        while j < next_right and next_counter >= min_samples:
            prev_counter, prev_left, prev_right = next_counter, next_left, next_right
            j += 1
            next_counter, next_left, next_right = count_samples(X, j, eps, sample_weight)

    
        if j == next_right:
            if next_counter >= min_samples:
                if X[j]- X[i] > 5:
                    for k in range(i, j + 1):
                        labels[k] = last_label + 1
                    last_label += 1

            else:
                if X[j - 1]- X[i] > 5:
                    for k in range(i, j):
                        labels[k] = last_label + 1
                    last_label += 1

        else:
            if X[j - 1] - X[i] > 5:
                for k in range(i, j):
                    labels[k] = last_label + 1
                last_label += 1
            
        i = j + 1

    return labels


def convert_labels_to_clusters(X, labels):
    clusters = {}

    labels_ = np.array(labels)

    for k in set(labels_):
        if k != -1:
            class_members = [index[0] for index in np.argwhere(labels_ == k)]
            start = min([X[t] for t in class_members])
            end = max([X[t] for t in class_members])

            if end < start:
                raise Exception('End of the cluster is smaller than the start!')

            clusters[start] = (start, end)

    return_list = []

    for key in sorted(clusters):
        return_list.append(clusters[key])

    return return_list

"""
def find_dbscan_borders(tn_clusters, chromosomes_length):
    borders = dict()

    for chrom in chrom_names:
        borders[chrom] = []

        start = -1

        for i in range(len(tn_clusters[chrom])):

            if start == -1:
                borders[chrom].append((1, tn_clusters[chrom][i][0] - 1))
                start = tn_clusters[chrom][i][1] + 1

            else:
                borders[chrom].append((start, tn_clusters[chrom][i][0] - 1))
                start = tn_clusters[chrom][i][1] + 1


            if i == len(tn_clusters[chrom]) - 1:
                if start < chromomosomes_length[chrom]:
                    borders[chrom].append((start, chromosomes_length[chrom]))

    return borders
"""

def dump_tn_clustering_labels(aligner, tn_clustering_labels, method, eps, min_samples):
    dump_file_name = '../analysis/clusters/{}/dump_{}_tn_clustering_labels_{}_{}.txt'.format(aligner, method,
                                                                                             eps, min_samples)
    with open(dump_file_name, 'w') as dumpfile:
        for chrom in tn_clustering_labels:
            for pos in tn_clustering_labels[chrom]:
                dumpfile.write('{}\t{}\t{}\n'.format(chrom, pos, tn_clustering_labels[chrom][pos]))


def dump_tn_clusters(aligner, tn_clusters, method, eps, min_samples):
    dump_file_name = '../analysis/clusters/{}/dump_{}_tn_clusters_{}_{}.txt'.format(aligner, method, eps, min_samples)
    with open(dump_file_name, 'w') as dumpfile:
        for chrom in tn_clusters:
            for cluster in tn_clusters[chrom]:
                dumpfile.write('{}\t{}\t{}\n'.format(chrom, cluster[0], cluster[1]))


def clustering(aligner, tn_position_scores, method, eps_limit, min_samples_limit):
    tn_clustering_labels = dict()
    tn_clusters = dict()
    for chrom in chrom_names:
        tn_clustering_labels[chrom] = {}
        tn_clusters[chrom] = {}

    for chrom in chrom_names:
        weights = []

        X = sorted(tn_position_scores[chrom])
        for x in X:
            weights.append(tn_position_scores[chrom][x])
    
        labels = find_dbscan_core_points(X, eps=eps_limit, min_samples=min_samples_limit, sample_weight=weights)

        for i in range(len(X)):
            tn_clustering_labels[chrom][X[i]] = labels[i]

        tn_clusters[chrom] = convert_labels_to_clusters(X, labels)

    dump_tn_clustering_labels(aligner, tn_clustering_labels, method, eps_limit, min_samples_limit)
    dump_tn_clusters(aligner, tn_clusters, method, eps_limit, min_samples_limit)



def main(aligner):
    #merging_methods = ['sum', 'linear', 'squareroot', 'log']
    merging_methods = ['linear']
    
    for method in merging_methods:
        file_name = '../analysis/dump_{}_merged_tn_position_scores_{}_{}.txt'.format(method, aligner, 'two_sided')
        tn_position_scores = serialization.read_tn_position_scores_from_dump(file_name)

        #for eps in range(15, 81, 5):
        #    for min_samples in range(10, 201, 2):
        for eps in range(15, 52, 2):
            for min_samples in range(6, 15, 2):
                clustering(aligner, tn_position_scores, method, eps, min_samples)


if __name__ == '__main__':
    aligner = sys.argv[1]
    main(aligner)
