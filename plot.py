
def plot1():
    stat_files = ['0_size_related_reads.txt', '1_size_sensitivity.txt', '2_size_blast_selected.txt',
                  '3_sim_read_coverage__gc_content.txt']

    x_label = {0: 'Pattern Size', 1: 'Pattern Size', 2: 'Pattern Size', 3: 'Simulated Read Coverage'}
    y_label = {0: 'Reads from the pattern', 1: 'Sensitivity', 2: 'Number of Selected Reads by Blast', 3: 'GC Content'}


    i = 0
    for file_name in stat_files:
        import matplotlib.pyplot as plt

        X = []
        Y = []
        with open(file_name) as input_file:
            lines = input_file.readlines()
            for line in lines:
                line = line.strip()
                if line == '':
                    continue
                nums = [float(n) for n in line.split()]
                X.append(nums[0])
                Y.append(nums[1])
        plt.plot(X, Y, 'o')
        plt.xlabel(x_label[i])
        plt.ylabel(y_label[i])
        plt.savefig('%s.png' % file_name)   # save the figure to file
        plt.close()
        i += 1

def plot2():
    dirs = ['original_case_r1_p1/', 'penalty_3_reward_1/']
    stat_files = ['1_size_sensitivity.txt', '2_size_blast_selected.txt']

    x_label = {1: 'Pattern Size', 2: 'Pattern Size'}
    y_label = {1: 'Sensitivity', 2: 'Number of Selected Reads by Blast'}

    X = []
    Y = [[], []]
    diagram_index = 1
    for file_name in stat_files:
        import matplotlib.pyplot as plt
        for i in range(len(dirs)):
            with open(dirs[i] + file_name) as input_file:
                lines = input_file.readlines()
                for line in lines:
                    line = line.strip()
                    if line == '':
                        continue
                    nums = [float(n) for n in line.split()]
                    if i == 0:
                        X.append(nums[0])
                    Y[i].append(nums[1])
        plt.plot(X, Y[0], 'o', label='same penalty and reward')
        plt.plot(X, Y[1], 'o', label='penalty = 3 * reward')
        plt.xlabel(x_label[diagram_index])
        plt.ylabel(y_label[diagram_index])
        plt.legend(loc=0)
        plt.savefig('compare_%s.png' % file_name)   # save the figure to file
        plt.close()
        diagram_index += 1


def plot_coverage_comparison():
    stat_files = ['10X_ratio.txt', '20X_ratio.txt', '30X_ratio.txt']
    coverages = [10, 20, 30]
    X = [[], [], []]
    Y = [[], [], []]

    import matplotlib.pyplot as plt

    for i in range(len(stat_files)):
        with open(stat_files[i]) as input_file:
            lines = input_file.readlines()
            for line in lines:
                line = line.strip()
                if line == '':
                    continue
                nums = [float(n) for n in line.split()]
                X[i].append(coverages[i])
                Y[i].append(nums[1])
    plt.plot(X[0], Y[0], 'o', label='10X')
    plt.plot(X[1], Y[1], 'o', label='20X')
    plt.plot(X[2], Y[2], 'o', label='30X')
    averages = []
    for i in range(3):
        sum = 0
        for e in Y[i]:
            sum += e
        averages.append(float(sum) / len(Y[i]))
    plt.plot(coverages, averages, label='avg')
    plt.xlabel('Coverage')
    plt.ylabel('Copy Count / True Copy Count')
    plt.legend(loc=0)
    plt.savefig('compare_%s.png' % 'coverages')  # save the figure to file
    plt.close()


def get_x_and_y_from_file(file_name, exclude_x=None):
    points = []
    X = []
    Y = []
    with open(file_name) as input_file:
        lines = input_file.readlines()
        for line in lines:
            line = line.strip()
            if line == '':
                continue
            nums = [float(n) for n in line.split()]
            points.append((nums[0], nums[1]))

    points.sort()
    for x, y in points:
        if exclude_x and x in exclude_x:
            continue
        X.append(x)
        Y.append(y)
    return X, Y


def get_pattern_result_map(file_name):
    res = {}
    min_len = 100
    max_len = 0
    with open(file_name) as input:
        lines = input.readlines()
        for line in lines:
            line = line.strip()
            if line == '':
                continue
            FP, sensitivity, evalue, p_num, len, TP = line.split()
            if float(len) > max_len:
                max_len = float(len)
            if float(len) < min_len:
                min_len = float(len)
            if p_num not in res.keys():
                res[p_num] = []
            res[p_num].append((int(FP), float(sensitivity), float(evalue), int(len), int(TP)))

    return res, min_len, max_len


def plot_sensitivity_over_fallout():
    stat_file = 'FP_and_sensitivity_evalue_min_len50.0.txt'
    # stat_file2 = 'fallout_and_sensitivity_min_len50.0_seq_68.txt'
    X, Y = get_x_and_y_from_file(stat_file)

    import matplotlib.pyplot as plt
    cmap = plt.get_cmap('seismic')

    pattern_results, min_len, max_len = get_pattern_result_map(stat_file)
    for p_num in pattern_results.keys():
        if int(p_num) % 5 != 4:
            continue
        # if p_num != '19':
        #     continue
        X = []
        Y = []
        points = []
        color = None
        for FP, sens, evalue, length, TP in pattern_results[p_num]:
            points.append((FP, sens))
            color = cmap((float(length) - min_len) / (max_len - min_len))
        points.sort()
        if points[len(points)-1][1] < 0.8:
            continue
        if points[0][1] > 0.3:
            continue
        for x, y in points:
            # if x > 1000:
            #     continue
            X.append(x)
            Y.append(y)
        plt.plot(X, Y, label='Pid:%s |Pattern|: %s' % (p_num, length))
    # plt.xscale('log')
    plt.xlabel('False Positives')
    plt.ylabel('Sensitivity')
#    plt.legend(loc=4, prop={'size':9})
    # plt.colorbar()
    plt.savefig('%s.png' % stat_file)  # save the figure to file
    plt.close()


def plot_tandem_copy_number_and_genome_copy_number():
    stat_file = 'copy_number_analysis.txt'
    sens_file = 'original_case_r1_p1/1_size_sensitivity.txt'
    X = []
    Y = []
    Y2 = []

    import matplotlib.pyplot as plt

    with open(stat_file) as input_file:
        lines = input_file.readlines()
        for line in lines:
            line = line.strip()
            if line == '':
                continue
            nums = [float(n) for n in line.split()]
            Y.append(float(nums[1]) / nums[0])

    with open(sens_file) as input_file:
        lines = input_file.readlines()
    for i in range(len(lines)):
        line = lines[i]
        line = line.strip()
        if line == '':
            continue
        nums = [float(n) for n in line.split()]
        X.append(nums[0])
    plt.plot(X, Y, 'o', color='blue', label='CN in 100Kbp region / CN in ~1Kbp region')
    # plt.plot(X, Y2, 'o', color='red', label='CN in 100Kbp region')
    plt.xlabel('Pattern size')
    plt.ylabel('Copies 100k / Copies in 1K')
    plt.legend(loc=0)
    plt.savefig('CN_in_genome_compare.png')  # save the figure to file
    plt.close()


def plot_reference_repeats():
    with open('vntrseek_repeats.txt') as input:
        lines = input.readlines()
        vntrseek_repeats = [int(float(num.strip())) for num in lines]
    with open('pattern_repeat_counts.txt') as input:
        lines = input.readlines()
        our_repeats = [int(num.strip()) for num in lines]
    import matplotlib.pyplot as plt
    X = [i+1 for i, j in enumerate(our_repeats)]
    plt.xlabel('Pattern ID')
    plt.ylabel('Number of Tandem Repeats in Reference Genome')
    plt.plot(X, vntrseek_repeats, '-o',color='blue', label='TRF Repeats')
    plt.plot(X, our_repeats, '--o', color='red', label="Our Method's Repeats")
    plt.legend(loc=0)
    plt.savefig('reference_repeats.png')
    plt.close()


def plot_copy_count_comparison(eliminated_nodes):
    import matplotlib.pyplot as plt
    from math import log
    X, vntr_coverage_ratio = get_x_and_y_from_file('vntr_coverage_ratio.txt', exclude_x=eliminated_nodes)
    _, hmm_y = get_x_and_y_from_file('hmm_repeat_count.txt', exclude_x=eliminated_nodes)
    hmm_y = [log(y, 2) for y in hmm_y]
    plt.xlabel('Pattern ID')
    plt.ylabel('Log of Computed Copy Number divided by True Copy Number')
    plt.plot(X, hmm_y, '--o', color='red', label="Log of Estimared Copy Count Divided by True Copy Count")
    # plt.ylim([0.6, 2])
    # plt.plot(X, vntr_coverage_ratio, color='green', label="VNTR Coverage Ratio in NGS Reads")
    plt.legend(loc=0, fontsize = 'x-small')
    plt.savefig('copy_count_comparison_log.png')
    plt.close()


def plot_FP_for_specific_sensitivity(eliminated_nodes, sensitivity=0.95):
    hmm_fps = {}
    blast_fps = {}
    with open('FP_and_sensitivity_evalue_min_len50.0.txt') as input:
        lines = [line.strip() for line in input.readlines() if line.strip() != '']
        for line in lines:
            FP, sens, _, pattern_id, pattern_len, _ = line.split()
            sens = float(sens)
            FP = int(FP)
            pattern_id = int(pattern_id)
            if pattern_id - 1 in eliminated_nodes:
                continue
            if sens >= sensitivity:
                if pattern_id not in blast_fps.keys():
                    blast_fps[pattern_id] = FP
                blast_fps[pattern_id] = min(blast_fps[pattern_id], FP)
    with open('FP_and_sensitivity_HMM_read_scoring_method.txt') as input:
        lines = [line.strip() for line in input.readlines() if line.strip() != '']
        for line in lines:
            FP, sens, _, pattern_id, pattern_len, _ = line.split()
            sens = float(sens)
            FP = int(FP)
            pattern_id = int(pattern_id)
            if pattern_id - 1 in eliminated_nodes:
                continue
            if sens >= sensitivity:
                if pattern_id not in hmm_fps.keys():
                    hmm_fps[pattern_id] = FP
                hmm_fps[pattern_id] = min(hmm_fps[pattern_id], FP)

    # X = sorted(list(set(hmm_fps.keys()) & set(blast_fps.keys())))
    X = sorted(list(hmm_fps.keys()))
    blast_fps_y = []
    hmm_fps_y = []
    for x in X:
        # blast_fps_y.append(blast_fps[x])
        hmm_fps_y.append(hmm_fps[x])
    import matplotlib.pyplot as plt
    plt.xlabel('Pattern ID')
    plt.ylabel('False Positives for Sensitivity of 0.9')
    # plt.plot(X, blast_fps_y, '-o',color='blue', label='BLAST False Positives')
    plt.plot(X, hmm_fps_y, '-o',color='red', label='False Positive Selected Reads')
    plt.legend(loc=0)
    plt.savefig('false_positives_for_sensitivity_of_09.png')
    plt.close()


def plot_coverage_ratio_histogram():
    m = {}

    from math import log
    with open('vntr_coverage_ratio.txt') as input:
        lines = input.readlines()
        for line in lines:
            index, ratio = line.strip().split()
            ratio = log(float(ratio), 2)
            if ratio not in m.keys():
                m[ratio] = []
            m[ratio].append(index)
    l = list(m.keys())
    xbins = [-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
    import matplotlib.pyplot as plt
    plt.hist(l, xbins, color='blue')
    plt.xlabel('log(coverage ratio)')
    plt.ylabel('# VNTR')
    plt.savefig('coverage_ratio.png')


edges = [(1, 8), (1, 16), (2, 17), (4, 18), (8, 16), (30, 32), (30, 33), (32, 33), (34, 40), (38, 57), (38, 59),
         (38, 67), (57, 59), (57, 67), (59, 67)] + [(5, 53), (47, 19), (71, 3), (31, 9)]
eliminated_nodes = []
for a, b in edges:
    if a-1 not in eliminated_nodes:
        eliminated_nodes.append(a-1)
    if b-1 not in eliminated_nodes:
        eliminated_nodes.append(b-1)

# plot_tandem_copy_number_and_genome_copy_number()
# plot_sensitivity_over_fallout()
# plot_reference_repeats()
# plot_copy_count_comparison(eliminated_nodes)
plot_FP_for_specific_sensitivity(eliminated_nodes)
