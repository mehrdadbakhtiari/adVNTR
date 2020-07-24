from time import time
import os
from random import gauss
import sys

import numpy as np
from keras.models import Sequential, load_model
from keras.layers import Dense, Activation

from Bio.SeqIO import SeqRecord
from Bio import SeqIO, Seq

from advntr.sam_utils import get_id_of_reads_mapped_to_vntr_in_bamfile, make_bam_and_index
from advntr.models import load_unique_vntrs_data
from advntr import settings


hg38 = True
PREFIX_DIR = '/mnt/hg38_dnn/' if hg38 else '/mnt/'
OUTPUT_DIR_PREFIX = '../advntr2_recruitment_comparison_hg38/' if hg38 else '../advntr2_recruitment_comparison/'

read_length = 150
kmer_length = 6
if __name__ == '__main__':
    if len(sys.argv) > 1:
        kmer_length = int(sys.argv[1])
# input_dim = 4 ** kmer_length * position_partition
input_dim = 4 ** kmer_length
# input_dim = read_length * 4
reduced_dimensions = 150 * (2 * kmer_length - 6)

losses = ['binary_crossentropy', 'mean_squared_error', 'mean_absolute_error', 'mean_squared_logarithmic_error', 'hinge', 'squared_hinge']
loss_to_activatio = {'binary_crossentropy': 'sigmoid',
                     'mean_squared_error': 'linear',
                     'mean_absolute_error': 'linear',
                     'mean_squared_logarithmic_error': 'linear',
                     'hinge': 'tanh',
                     'squared_hinge': 'tanh'}
loss_index = 0
if __name__ == '__main__':
    if len(sys.argv) > 4:
        loss_index = int(sys.argv[4])
loss_function = losses[loss_index]
loss_suffix = '_%s' % loss_index if loss_index > 0 else ''

result_dir = OUTPUT_DIR_PREFIX + 'hmm_dnn_comparison_%s/' % (str(kmer_length) + loss_suffix)
bowtie_result_dir = OUTPUT_DIR_PREFIX + 'hmm_dnn_comparison_bowtie/'
bowtie_working_dir = PREFIX_DIR + '/bowtie_recruitment/'
dnn_models_dir = PREFIX_DIR + '/dnn_models_%s/' % (str(kmer_length)) + loss_suffix


def align_with_bowtie(fq_file):
    bowtie_alignment = fq_file[:-3] + '_bowtie_aln.sam'
    if not os.path.exists(bowtie_alignment[:-4] + '.bam'):
        os.system('bowtie2 -x /mnt/hg19_chromosomes/hg19_bt2_idx -f %s -S %s --threads 7' % (fq_file, bowtie_alignment))
        make_bam_and_index(bowtie_alignment)
    return bowtie_alignment[:-4] + '.bam'

def get_embedding_of_string(sequence, kmer_length=6):
    input_dim = 4 ** kmer_length
    sequence = sequence.upper()
    num = 0
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for c in 'ASDFGHJKLPOIUYTREWQZXCVBNM':
        if c not in mapping.keys():
            mapping[c] = 0
    for i in range(len(sequence[:kmer_length])):
        num += mapping[sequence[i]] * (4 ** (kmer_length - i - 1))
    result = [0] * input_dim
    result[num] = 1
    # result = set()
    # result.add(num)
    highest_position = 4 ** (kmer_length-1)
    for i in range(kmer_length, len(sequence)):
        num -= highest_position * mapping[sequence[i-kmer_length]]
        num *= 4
        num += mapping[sequence[i]]
        result[num] = 1
        # result.add(num)
    return result


def get_google_embedding_of_string(sequence):
    sequence = sequence.upper()
    result = [0] * input_dim
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i, c in enumerate(sequence):
        if c not in mapping.keys():
            mapping[c] = 0
        result[i * 4 + mapping[c]] = 1
    return result


def make_random_unit_vector(dims):
    vec = [gauss(0, 1) for i in range(dims)]
    mag = sum(x ** 2 for x in vec) ** .5
    return [x / mag for x in vec]


def get_random_vector_set(seed=10):
    import random
    random.seed(seed)

    random_vectors = []
    for _ in range(reduced_dimensions):
        random_vectors.append(make_random_unit_vector(input_dim))
    return random_vectors


def get_hashed_embedding_of_string(sequence, random_vectors):
    original_embedding = get_embedding_of_string(sequence)
    hashed_embedding = []
    for i, random_vector in enumerate(random_vectors):
        hashed_embedding.append(0)
        for position in list(original_embedding):
            hashed_embedding[i] += random_vector[position]
    # hashed_embedding = np.array(random_vectors).dot(np.array(original_embedding))
    return hashed_embedding


def generate_d_neighborhood(pattern, d):
    neigh = set([pattern])
    for i in range(d):
        addition = set([])
        for p in neigh:
            for j in range(len(p)):
                insertion = [p[j] + inserted for inserted in 'ACGT']
                for sub in [''] + ['A', 'C', 'G', 'T'] + insertion:
                    new_str = p[:j] + sub + p[j+1:]
                    addition.add(new_str)
        neigh |= addition
    return neigh


def get_blast_keywords(reference_vntr, keyword_size=11):
    vntr = ''.join(reference_vntr.get_repeat_segments())
    if len(vntr) < keyword_size:
        min_copies = int(keyword_size / len(vntr)) + 1
        vntr = str(vntr) * min_copies
    locus = reference_vntr.left_flanking_region[-15:] + vntr + reference_vntr.right_flanking_region[:15]
    queries = []
    step_size = 5 if len(reference_vntr.pattern) != 5 else 6
    for i in range(0, len(locus) - keyword_size + 1, step_size):
        queries.append(locus[i:i+keyword_size])
    return queries


def get_hmm_accuracy(vntr_finder, simulated_true_reads, simulated_false_filtered_reads):
    output_dir = result_dir + '/%s/' % vntr_finder.reference_vntr.id

    print('running BLAST')
    from blast_wrapper import get_blast_matched_ids, make_blast_database
    blast_dir = output_dir + 'blast_dir/'
    if not os.path.exists(blast_dir):
        os.makedirs(blast_dir)
    vntr_id = vntr_finder.reference_vntr.id
    fasta_file = blast_dir + 'reads.fasta'
    records = []
    for i, read in enumerate(simulated_false_filtered_reads):
        records.append(SeqRecord(seq=Seq.Seq(read), id='fasle_%s' % i))
    for i, read in enumerate(simulated_true_reads):
        records.append(SeqRecord(seq=Seq.Seq(read), id='true_%s' % i))
    with open(fasta_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

    make_blast_database(fasta_file, blast_dir + 'blast_db_%s' % vntr_id)
    query = '@'.join(get_blast_keywords(vntr_finder.reference_vntr))
    search_id = 'search_id'
    search_results = get_blast_matched_ids(query, blast_dir + 'blast_db_%s' % vntr_id, max_seq='100000', word_size='7',
                                           evalue=sys.maxsize, search_id=search_id, identity_cutoff='100', blast_tmp_dir=blast_dir)
    from collections import Counter
    res = Counter(search_results)
    filtered = [item for item, occur in res.items() if occur >= 2]
    print('BLAST results computed')

    print(len(filtered))
    print(len(simulated_true_reads))
    print(len(simulated_false_filtered_reads))
    tp = float(len([e for e in filtered if e.startswith('true')]))
    fp = float(len([e for e in filtered if e.startswith('false')]))
    fn = float(len(simulated_true_reads) - tp)
    tn = float(len(simulated_false_filtered_reads) - fp)
    train_time = 0
    passed_time = 0

    precision = tp / (tp + fp) if tp > 0 else 0
    recall = tp / (tp + fn)
    accuracy = (100 * (tp + tn) / (fp + fn + tp + tn))
    print('BLAST:')
    print(tp, fp, fn, tn)
    print('Precision:', precision)
    print('Recall:', recall)
    print('acc: %s' % accuracy)

    with open(output_dir + '/blast.txt', 'w') as outfile:
        outfile.write('%s\n' % train_time)
        outfile.write('%s\n' % passed_time)
        outfile.write('%s\n' % precision)
        outfile.write('%s\n' % recall)
        outfile.write('%s\n' % accuracy)
        outfile.write('%s,%s,%s,%s\n' % (tp, fn, fp, tn))
    return passed_time

    output_dir = result_dir + '/%s/' % vntr_finder.reference_vntr.id
    if os.path.exists(output_dir + '/hmm.txt') and os.path.getsize(output_dir + '/hmm.txt') > 0:
        if sum(1 for _ in open(output_dir + 'hmm.txt')) > 5:
            print('HMM info is already calculated')
            with open(output_dir + 'hmm.txt') as infile:
                lines = infile.readlines()
                return float(lines[1])

    train_true_reads = [read for i, read in enumerate(simulated_true_reads) if i % 2 == 0]
    train_false_reads = [read for i, read in enumerate(simulated_false_filtered_reads) if i % 2 == 0]
    test_true_reads = [read for i, read in enumerate(simulated_true_reads) if i % 2 == 1]
    test_false_reads = [read for i, read in enumerate(simulated_false_filtered_reads) if i % 2 == 1]

    start_time = time()
    hmm = vntr_finder.get_vntr_matcher_hmm(read_length=read_length)

    processed_true_reads = vntr_finder.find_hmm_score_of_simulated_reads(hmm, train_true_reads)
    processed_false_reads = vntr_finder.find_hmm_score_of_simulated_reads(hmm, train_false_reads)

    recruitment_score = vntr_finder.find_recruitment_score_threshold(processed_true_reads, processed_false_reads)
    train_time = time() - start_time
    print('HMM train time: %s' % train_time)

    tp = 0.0
    fn = 0.0
    tn = 0.0
    fp = 0.0
    start_time = time()
    true_reads = vntr_finder.find_hmm_score_of_simulated_reads(hmm, test_true_reads)
    false_reads = vntr_finder.find_hmm_score_of_simulated_reads(hmm, test_false_reads)
    passed_time = time() - start_time
    for read in true_reads:
        if read.logp > recruitment_score:
            tp += 1
        else:
            fn += 1
    for read in false_reads:
        if read.logp > recruitment_score:
            fp += 1
        else:
            tn += 1
    precision = tp / (tp + fp) if tp > 0 else 0
    recall = tp / (tp + fn)
    accuracy = (100 * (tp + tn) / (fp + fn + tp + tn))
    print('HMM: %s' % passed_time)
    print(tp, fp, fn, tn)
    print('Precision:', precision)
    print('Recall:', recall)
    print('acc: %s' % accuracy)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + '/hmm.txt', 'w') as outfile:
        outfile.write('%s\n' % train_time)
        outfile.write('%s\n' % passed_time)
        outfile.write('%s\n' % precision)
        outfile.write('%s\n' % recall)
        outfile.write('%s\n' % accuracy)
        outfile.write('%s,%s,%s,%s\n' % (tp, fn, fp, tn))
    return passed_time


def simulate_and_store_false_reads(vntr_finder, false_reads_file, min_matches=4):
    simulated_false_filtered_reads = []
    reference_files = []
    for chromosome in settings.CHROMOSOMES:
        reference_file = settings.HG19_DIR + chromosome + '.fa'
        reference_files.append(reference_file)

    if hg38:
        reference_files = ['/mnt/hg38_chromosomes/hg38.fa']
    for reference_file in reference_files:
        simulated_false_filtered_reads += vntr_finder.simulate_false_filtered_reads(reference_file, min_matches)
        print(len(simulated_false_filtered_reads))
        if len(simulated_false_filtered_reads) > 41000:
            break
    with open(false_reads_file, 'w') as outfile:
        for read in simulated_false_filtered_reads:
            outfile.write('%s\n' % read)


def get_true_reads_and_false_reads(vntr_finder, vntr_id):
    simulated_true_reads = vntr_finder.simulate_true_reads(read_length)
    print('true reads: %s' % len(simulated_true_reads))
    false_reads_file = PREFIX_DIR + '/false_reads/false_reads_%s.txt' % vntr_id
    if not os.path.exists(false_reads_file) or os.path.getsize(false_reads_file) == 0:
        if os.path.exists(false_reads_file) and os.path.getsize(false_reads_file) == 0:
            print('There is no false read in the file')
            no_false_read = True
        else:
            no_false_read = False
        min_matches = 1 if no_false_read else 4
        simulate_and_store_false_reads(vntr_finder, false_reads_file, min_matches)

    min_matches = 6
    while True:
        with open(false_reads_file) as infile:
            lines = infile.readlines()
        if len(lines) > 40000:
            print('There are more than %s reads in the file. Trying min_matches = %s' % (len(lines), min_matches))
            simulate_and_store_false_reads(vntr_finder, false_reads_file, min_matches)
            min_matches += 2
        else:
            break
    simulated_false_filtered_reads = [read.strip() for read in lines if 'N' not in read.upper()]
    print('false reads: %s' % len(simulated_false_filtered_reads))
    print('true reads: %s' % len(simulated_true_reads))

    return simulated_true_reads, simulated_false_filtered_reads


def get_nn_model(train, three_hidden_layers=False, model_function='relu', first_layer=100, second_layer=0):
    model = Sequential()
    model.add(Dense(first_layer, input_dim=input_dim, kernel_initializer="uniform", activation=model_function))
    # if three_hidden_layers:
    if second_layer > 0:
        model.add(Dense(second_layer, activation=model_function, kernel_initializer="uniform"))
    model.add(Dense(2))
    model.add(Activation("softmax"))

    model.compile(loss=loss_function, optimizer='adam', metrics=['accuracy'])
    model.fit(train[0], train[1], epochs=3, batch_size=10)
    return model


def is_true(result_class):
    return result_class[0] > result_class[1]


def select_positive_and_negative_reads_with_bowtie(reads, vntr_finder, label):
    working_dir = bowtie_working_dir + '/%s/' % vntr_finder.reference_vntr.id
    if not os.path.exists(working_dir):
        os.makedirs(working_dir)
    fq_file = working_dir + label + '.fa'

    records = []
    for i, read in enumerate(reads):
        record = SeqRecord('')
        record.seq = Seq.Seq(read)
        record.id = 'read_%s/1' % str(i)
        records.append(record)
    with open(fq_file, 'w') as output_handle:
        SeqIO.write(records, output_handle, 'fasta')

    passed_time = time()
    bowtie_bamfile = align_with_bowtie(fq_file)
    bowtie_selected = len(get_id_of_reads_mapped_to_vntr_in_bamfile(bowtie_bamfile, vntr_finder.reference_vntr))
    return float(bowtie_selected), float(len(reads) - bowtie_selected), time() - passed_time


def run_bowtie2(true_reads, false_reads, vntr_finder):
    output_dir = bowtie_result_dir + '%s/' % vntr_finder.reference_vntr.id
    if os.path.exists(output_dir + '/bowtie.txt') and os.path.getsize(output_dir + '/bowtie.txt') > 0:
        if sum(1 for _ in open(output_dir + 'bowtie.txt')) > 5:
            print('bowtie results is already computed')
            return

    train_time = 0
    tp, fn, t1 = select_positive_and_negative_reads_with_bowtie(true_reads, vntr_finder, 'true')
    fp, tn, t2 = select_positive_and_negative_reads_with_bowtie(false_reads, vntr_finder, 'false')
    passed_time = t1 + t2
    precision = tp / (tp + fp) if tp > 0 else 0
    recall = tp / (tp + fn)
    accuracy = float(tp + tn) / (tp + tn + fp + fn)

    print('Bowtie2: %s' % passed_time)
    print('Precision:', precision)
    print('Recall:', recall)
    print('acc: %s' % accuracy)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + '/bowtie.txt', 'w') as outfile:
        outfile.write('%s\n' % train_time)
        outfile.write('%s\n' % passed_time)
        outfile.write('%s\n' % precision)
        outfile.write('%s\n' % recall)
        outfile.write('%s\n' % accuracy)
        outfile.write('%s,%s,%s,%s\n' % (tp, fn, fp, tn))


def run_simulation(vntr_map, vntr_id):
    print('vntr:', vntr_id)
    ref_vntr = vntr_map[vntr_id]
    from advntr.vntr_finder import VNTRFinder
    vntr_finder = VNTRFinder(ref_vntr)
    simulated_true_reads, simulated_false_filtered_reads = get_true_reads_and_false_reads(vntr_finder, vntr_id)
    if len(simulated_false_filtered_reads) > 30000 or len(simulated_false_filtered_reads) == 0:
        print('skipping VNTR', vntr_id)
        return

    # map with bowtie2
#    run_bowtie2(simulated_true_reads, simulated_false_filtered_reads, vntr_finder)

#    hmm_time = get_hmm_accuracy(vntr_finder, simulated_true_reads, simulated_false_filtered_reads)

    if not os.path.exists(dnn_models_dir):
        os.makedirs(dnn_models_dir)
    output_dir = result_dir + '/%s/' % vntr_finder.reference_vntr.id
    model_dir = dnn_models_dir + '%s.hd5' % vntr_finder.reference_vntr.id
    if os.path.exists(output_dir + 'dnn.txt') and os.path.getsize(output_dir + 'dnn.txt') > 0:
        if sum(1 for _ in open(output_dir + 'dnn.txt')) > 5:
            print('dnn information is already calculated')
            return

    true_embeddings = [get_embedding_of_string(seq) for seq in simulated_true_reads]
    false_embeddings = [get_embedding_of_string(seq) for seq in simulated_false_filtered_reads]
    train_x = [embedding for i, embedding in enumerate(true_embeddings) if i % 2 == 0]
    start_time = time()
    test_x = [embedding for i, embedding in enumerate(true_embeddings) if i % 2 == 1]
    embedding_time = time() - start_time
    train_y = [[1, 0]] * len(train_x)
    test_y = [[1, 0]] * len(test_x)
    train_x += [embedding for i, embedding in enumerate(false_embeddings) if i % 2 == 0]
    start_time = time()
    test_x += [embedding for i, embedding in enumerate(false_embeddings) if i % 2 == 1]
    embedding_time += time() - start_time
    train_y += [[0, 1]] * (len(train_x) - len(train_y))
    test_y += [[0, 1]] * (len(test_x) - len(test_y))
    train = [np.array(train_x), np.array(train_y)]
    test = [np.array(test_x), np.array(test_y)]

    first_layer = 100
    second_layer = 50
    start_time = time()
    if os.path.exists(model_dir):
        print('DNN model is already trained')
        model = load_model(model_dir)
    else:
        model = get_nn_model(train, first_layer=first_layer, second_layer=second_layer)
        model.save(model_dir)
    train_time = time() - start_time
    print('NN train time: %s' % train_time)

    scores = model.evaluate(test[0], test[1])
    start_time = time()
    classes = model.predict(test[0], batch_size=128)
    passed_time = embedding_time + time() - start_time
    passed_time += hmm_time / len(test[0]) * len(true_embeddings) / 2
    fn = 0.0
    fp = 0.0
    tp = 0.0
    tn = 0.0
    for i in range(len(test[1])):
        majority = int(is_true(classes[i]))# + int(is_true(classes2[i])) + int(is_true(classes3[i]))
        # print(majority)
        if test[1][i][0] == 1:
            if majority >= 1:
                tp += 1
            else:
                fn += 1
        else:
            if majority >= 1:#is_true(classes[i]):
                fp += 1
            else:
                tn += 1
    precision = tp / (tp + fp) if tp > 0 else 0
    recall = tp / (tp + fn)
    accuracy = scores[1]*100

    print('NN: %s' % passed_time)
    print(tp, fp, fn, tn)
    print('Precision:', precision)
    print('Recall:', recall)
    print("\n%s: %.2f%%" % (model.metrics_names[1], accuracy))

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    with open(output_dir + '/dnn.txt', 'w') as outfile:
        outfile.write('%s\n' % train_time)
        outfile.write('%s\n' % passed_time)
        outfile.write('%s\n' % precision)
        outfile.write('%s\n' % recall)
        outfile.write('%s\n' % accuracy)
        outfile.write('%s,%s,%s,%s\n' % (tp, fn, fp, tn))


def main():
    vntr_map = {}
    if hg38:
        reference_vntrs = load_unique_vntrs_data('vntr_data/hg38_selected_VNTRs_Illumina.db')
        vntr_ids = []
        for ref_vntr in reference_vntrs:
            vntr_map[ref_vntr.id] = ref_vntr
            if 100 >= len(ref_vntr.pattern) >= 6:
                vntr_ids.append(ref_vntr.id)
    else:
        reference_vntrs = load_unique_vntrs_data()
        for ref_vntr in reference_vntrs:
            vntr_map[ref_vntr.id] = ref_vntr

        from advntr.advntr_commands import get_tested_vntrs
        vntr_ids = get_tested_vntrs()

    print('len of reference_vntrs:', len(reference_vntrs))
    print('# of vntrs: %s' % len(vntr_ids))

    start, end = int(sys.argv[2]), int(sys.argv[3])

    # run_simulation(vntr_map, 503431)
    # exit(0)

    count = 0
    for vid in vntr_ids:
        count += 1
        if count < start or count > end:
            continue
        run_simulation(vntr_map, vid)

    # best_f, best_s, best_acc = None, None, 0
    # with open('param_training2.txt', 'w') as output:
    #     for f in accuracy_map.keys():
    #         for s in accuracy_map[f].keys():
    #             avg_accuracy = sum(accuracy_map[f][s]) / len(accuracy_map[f][s])
    #             output.write('%s %s %s\n' % (f, s, avg_accuracy))
    #             if avg_accuracy > best_acc:
    #                 best_acc = avg_accuracy
    #                 best_f = f
    #                 best_s = s
    # print(best_f, best_s, best_acc)

if __name__ == '__main__':
    main()
