
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


def get_numbers_from_file(file_name):
    with open(file_name) as input_file:
        lines = input_file.readlines()
        result = [float(line.strip()) for line in lines]
    return result


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
    print(X)
    print(hmm_fps_y)
    # plt.plot(X, hmm_fps_y, '-o',color='red', label='False Positive Selected Reads')
    plt.legend(loc=0)
    # plt.savefig('false_positives_for_sensitivity_of_09.png')
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


def plot_gc_content_violin_plot():
    from coverage_bias import CoverageBiasDetector
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bias_detector = CoverageBiasDetector('original_reads/paired_dat.sam')
    gc_coverage_map = bias_detector.get_gc_content_coverage_map()
    data = []
    pos = []
    for gc_content in range(0, 10):
        pos.append(gc_content * 10)
        data.append([float('nan'), float('nan')])
        if gc_content in gc_coverage_map:
            data[gc_content] = gc_coverage_map[gc_content]
    plt.violinplot(data, pos, widths=7, showmeans=True)
    plt.xlabel('GC Content Percentage')
    plt.ylabel('Coverage')
    plt.savefig('gc_coverage_violinplot_simulated.png')


def plot_frequency_of_repeats_in_population():
    import numpy as np
    import matplotlib.pyplot as plt
    R5_array = (45, 51, 56)
    R4_array = (33, 10, 20)
    R2_array = (22, 18, 12)
    R3_array = (0, 21, 12)
    R5 = np.array(R5_array)
    R4 = np.array(R4_array)
    R2 = np.array(R2_array)
    R3 = np.array(R3_array)
    ind = (0, 1, 2)
    width = 0.35
    p5 = plt.bar(ind, R5, width)
    p4 = plt.bar(ind, R4, width, bottom=R5)
    p2 = plt.bar(ind, R2, width, bottom=R5+R4)
    p3 = plt.bar(ind, R3, width, bottom=R5+R4+R2)
    plt.plot([-0.5, 3], [45, 45], color='k', linestyle='--', linewidth=0.5)
    plt.plot([-0.5, 3], [56, 56], color='k', linestyle='--', linewidth=0.5)
    plt.ylabel('Frequency')
    plt.title('Allele Frequency for VNTR in GP1BA')
    plt.xticks(ind, ('African', 'East Asian', 'European'))
    plt.yticks(np.arange(0, 150, 10))
    plt.legend((p5[0], p4[0], p2[0], p3[0]), ('3 Repeats', '4 Repeats', '1 Repeats', '2 Repeats'))
    plt.savefig('GP1BA.png', dpi=300)


def get_diabetes_pattern_interavls():
    pattern = 'GGCCCCCCCCGTGCCGCCCACGGGTGACTCCGG'
    last = pattern[0]
    start = 0
    intervals = []
    for i in range(1, len(pattern)+1):
        if i == len(pattern):
            intervals.append((start+1, i))
            break
        if pattern[i] != last:
            intervals.append((start+1, i))
            last = pattern[i]
            start = i
    return intervals


def plot_indel_frequencies_for_diabetes():
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib import rc, rcParams
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = '#FFFFFF'
    rc('text', usetex=True)
    rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    plt.title('Frameshift Frequency in Diabetes Patients')
    plt.xlabel(r'\emph{Frameshift Position}')
    plt.ylabel(r'\emph{\# of Individuals}')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    plt.tight_layout(pad=4.0, w_pad=0.5, h_pad=2.0)

    raw_case = [('I0C', 1), ('D16', 1), ('D10', 1), ('I16T', 1), ('I28G', 1), ('I13C,D10', 1), ('I31A,D30,D31,I33T', 1), ('I9A', 1), ('I20C', 1), ('I3A', 1), ('I5C', 1), ('I19C', 1), ('I26G', 1), ('I31A', 1), ('I23C', 1), ('I33G', 1), ('I33C', 1), ('I2T', 1), ('I8T', 1), ('D13,I13C', 2), ('I24C', 2), ('I21A', 2), ('I10C', 3), ('I33T', 3), ('I12C', 4), ('D5', 5)]
    raw_control = [('D16', 1), ('D10', 1), ('I33C', 1), ('I7A', 1), ('I8C', 1), ('I24C', 1), ('I33C,I25T', 1), ('I28C', 1), ('I19G', 1), ('I13G', 1), ('I9C', 1), ('I11C', 1), ('I6C', 2), ('D5', 2), ('I10C', 3), ('I12C', 7)]
    case = {}
    control = {}
    for keys, value in raw_case:
        for key in keys.split(','):
            if key in case.keys():
                case[key] += value
            else:
                case[key] = value

    for keys, value in raw_control:
        for key in keys.split(','):
            if key in control.keys():
                control[key] += value
            else:
                control[key] = value

    total_indes = [pos for pos, num in case.items()]
    total_indes += [pos for pos, num in control.items()]
    total_indes = list(set(total_indes))
    total_indel_sorted = []
    intervals = get_diabetes_pattern_interavls()
    print intervals

    width = 0.35
    case_array = []
    control_array = []
    filtered_indels = []
    for pos in total_indes:
        case_count = case[pos] if pos in case.keys() else 0
        control_count = control[pos] if pos in control.keys() else 0
        if case_count + control_count < 2:
            continue
        case_array += [case_count]
        control_array += [control_count]
        filtered_indels += [pos]
        print case_count, control_count
    case_array = np.array(case_array)
    control_array = np.array(control_array)
    print(case_array)
    print(control_array)
    print(filtered_indels)

    filtered_indels = ['[3-10]D', '[3-10]I\_C', '[12]I\_C', '[13]D', '[13]I\_C', '[16]D', '[21]I\_A', '[22-24]I\_C', '[30-31]I\_A', '[32-33]I\_C','[32-33]I\_T']
    case_array = np.array([5+2, 3, 4, 2, 3, 1, 2, 1, 2, 1, 4])
    control_array = np.array([2+1, 2+3, 7, 0, 0, 1, 0, 1, 0, 2, 0])

    x_axis = [i for i, indel in enumerate(filtered_indels)]
    p0 = plt.bar(x_axis, case_array, width)
    p1 = plt.bar(x_axis, control_array, width, bottom=case_array)
    plt.xticks(x_axis, filtered_indels, fontsize=6, rotation=45)
    plt.legend((p0[0], p1[0]), ('Case', 'Control'))
    plt.savefig('diabetes_indels.png', dpi=300)


def add_recruitment_results_for_illumina(illumina_recruitment_plots, results_dir):
    import glob
    titles = 'ABC'

    gene_dirs = glob.glob(results_dir + '*')
    gene_index = 0
    for gene_dir in gene_dirs:
        gene_name = gene_dir.split('/')[-1]
        result_file = gene_dir + '/result.txt'
        if gene_name == 'IL1RN' or gene_name == 'DRD4':
            continue
        copies= []
        bwa_result = []
        bowtie_result = []
        our_selection_result = []
        with open(result_file) as input:
            lines = input.readlines()
            for line in lines:
                copy, original, our_filtering, our_selection, bwa, bowtie = line.split()
                original = int(original)
                copies.append(copy)
                our_selection_result.append(float(our_selection) / original)
                bwa_result.append(float(bwa) / original)
                bowtie_result.append(float(bowtie) / original)
        illumina_recruitment_plots[gene_index].title.set_text(titles[gene_index] + ') %s' % gene_name)
        illumina_recruitment_plots[gene_index].plot(copies, our_selection_result, 'o-', label='Our Method')
        illumina_recruitment_plots[gene_index].plot(copies, bwa_result, '.-', label='BWA-MEM')
        illumina_recruitment_plots[gene_index].plot(copies, bowtie_result, '.-', label='Bowtie2', color='orange')
        illumina_recruitment_plots[gene_index].spines['bottom'].set_color('black')
        illumina_recruitment_plots[gene_index].spines['left'].set_color('black')
        location = 8 if gene_name != 'CSTB' else 3
        gene_index += 1


def add_recruitment_results_for_pacbio(pacbio_recruitment_plots, results_dir):
    from matplotlib.ticker import FormatStrFormatter
    import numpy

    import glob
    titles = 'DEF'

    gene_dirs = glob.glob(results_dir + '*')
    gene_index = 0
    for gene_dir in gene_dirs:
        gene_name = gene_dir.split('/')[-1]
        result_file = gene_dir + '/result.txt'
        copies= []
        bwa_result = []
        bowtie_result = []
        our_selection_result = []
        with open(result_file) as input:
            lines = input.readlines()
            for line in lines:
                copy, original, our_filtering, our_selection, bwa, bowtie = line.split()
                original = max(29, int(our_selection))
                copies.append(copy)
                our_selection_result.append(float(our_selection) / original)
                bwa_result.append(float(bwa) / original)
                bowtie_result.append(float(bowtie) / original)
        pacbio_recruitment_plots[gene_index].title.set_text(titles[gene_index] + ') %s' % gene_name)
        pacbio_recruitment_plots[gene_index].plot(copies, our_selection_result, '.-', label='Our Method')
        pacbio_recruitment_plots[gene_index].plot(copies, bowtie_result, '.-', label='Blasr', color='#00BA38')
        pacbio_recruitment_plots[gene_index].spines['bottom'].set_color('black')
        pacbio_recruitment_plots[gene_index].spines['left'].set_color('black')
        # pacbio_recruitment_plots[gene_index].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        pacbio_recruitment_plots[gene_index].yaxis.set_ticks([0.8, 0.9, 1.0])

        gene_index += 1


def plot_read_selection_and_mapping_sensitivity_comparison():
    from matplotlib import rc, rcParams
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = '#FFFFFF'
    rc('text', usetex=True)
    rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    plt.title('Read Recruitment Comparison')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')

    fig = plt.figure()
    ax = list([])
    x_label_font = 12
    y_label_font = 12

    ax.append(fig.add_subplot(111))
    ax[0].set_ylabel(r'\emph{Read Selection Recall}', fontsize=y_label_font)
    ax[0].set_xlabel(r'\emph{Simulated RU Count}', fontsize=x_label_font)

    # Turn off axis lines and ticks of the big subplot
    for i in range(1):
        ax[i].spines['top'].set_color('none')
        ax[i].spines['bottom'].set_color('none')
        ax[i].spines['left'].set_color('none')
        ax[i].spines['right'].set_color('none')
        ax[i].tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')


    illumina_recruitment_plots = list([])
    illumina_recruitment_plots.append(fig.add_subplot(231))
    illumina_recruitment_plots.append(fig.add_subplot(232))
    illumina_recruitment_plots.append(fig.add_subplot(233))
    add_recruitment_results_for_illumina(illumina_recruitment_plots, results_dir='../Illumina_copy_number_short_vntrs_mapping/')

    pacbio_recruitment_plots = list([])
    pacbio_recruitment_plots.append(fig.add_subplot(234))
    pacbio_recruitment_plots.append(fig.add_subplot(235, sharey=pacbio_recruitment_plots[0]))
    pacbio_recruitment_plots.append(fig.add_subplot(236, sharey=pacbio_recruitment_plots[0]))
    add_recruitment_results_for_pacbio(pacbio_recruitment_plots, results_dir='../pacbio_coverage_experiment/')

    plt.tight_layout(pad=0.6, w_pad=0.5, h_pad=1.0)
    # plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=0.4)
    plt.subplots_adjust(top=0.85, left=0.1, bottom=0.1)

    illumina_handles, illumina_labels = illumina_recruitment_plots[2].get_legend_handles_labels()
    handles, labels = pacbio_recruitment_plots[2].get_legend_handles_labels()
    plt.figlegend(handles + illumina_handles[1:], labels + illumina_labels[1:], loc='upper center', ncol=5, labelspacing=0.)
    # fig.legend(lines, labels, loc=(0.5, 0), ncol=5)

    plt.savefig('read_recruitment_result.png', dpi=1000)


def plot_ins_simulation_pacbio_results(results_dir='out/'):
    # Deprecated
    import glob
    files = glob.glob(results_dir + 'out_INS*')
    points = []
    for file_name in files:
        sim = int(file_name.split('_')[2])
        with open(file_name) as input:
            lines = input.readlines()
            estimate = int(lines[-1].strip())
        points.append((sim, estimate))
    points = sorted(points)
    sim_repeats = [sim for sim, estimate in points]
    estimated_repeats = [estimate for sim, estimate in points]
    import matplotlib.pyplot as plt
    plt.plot(sim_repeats, estimated_repeats, 'o')
    plt.title('Result of estimation on PacBio simulated reads')
    plt.ylabel('Estimated Copy Number')
    plt.xlabel('Simulated Copy Number')
    plt.savefig('INS_simulation_results.png', dpi=300)


def plot_pacbio_ru_length_result(results_dir='../pacbio_ru_data_for_all_vntrs/'):
    from matplotlib import rc, rcParams
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = '#FFFFFF'
    rc('text', usetex=True)
    rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    plt.title('Effect of RU Length on Copy Number Estimation')
    plt.ylabel(r'\emph{Correct Estimated Copy Number}')
    plt.xlabel(r'\emph{RU Length}')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')
    import glob
    ru_dirs = glob.glob(results_dir + '*')
    points = []
    for ru_dir in ru_dirs:
        ru_length = ru_dir.split('/')[-1]
        files = glob.glob(ru_dir + '/*.fasta.out')
        corrects = 0
        for file_name in files:
            sim = int(file_name.split('_')[-2])
            correct = False
            with open(file_name) as input:
                lines = input.readlines()
                if len(lines) > 1:
                    if lines[-1].strip() != 'None' and len(lines[-1]) < 10:
                        estimate = int(lines[-1].strip())
                        if estimate == sim:
                            correct = True
            if not correct and ru_length == 20:
                print file_name
            if correct:
                corrects += 1
        points.append((ru_length, corrects))
    points = sorted(points)
    print points
    ru_lengths = [x for x, y in points]
    corrects = [y for x, y in points]
    plt.plot(ru_lengths, corrects, label='Correct Estimates')

    plt.legend(loc=0, fontsize='x-small')
    plt.savefig('pacbio_ru_length_results.png', dpi=300)


def get_correct_estimates(files):
    res = [0 for i in range(1, 40)]
    for file_name in files:
        coverage = int(file_name.split('_')[-1].split('.')[0][:-1])
        sim = int(file_name.split('_')[-2])
        correct = False
        with open(file_name) as input:
            lines = input.readlines()
            if len(lines) > 1:
                if lines[-1].strip() != 'None' and len(lines[-1]) < 10:
                    estimate = int(float(lines[-1].strip()))
                    if estimate == sim:
                        correct = True
        if correct:
            res[coverage - 1] += 1
    return res


def add_coverages_for_three_genes(coverage_plots, results_dir):
    import glob
    gene_dirs = glob.glob(results_dir + '*')
    coverages_label = [i for i in range(1, 40)]
    shapes = ['-', '--', '--']
    shape = 0
    for gene_dir in gene_dirs:
        gene_name = gene_dir.split('/')[-1]
        files = glob.glob(gene_dir + '/*.fasta.out')
        coverages_of_gene = get_correct_estimates(files)
        coverages_of_gene_naive = get_correct_estimates(glob.glob(gene_dir + '/*.fasta.out.naive'))
        if gene_name == 'CSTB':
            coverages_of_gene, coverages_of_gene_naive = coverages_of_gene[:7] + coverages_of_gene_naive[7:], coverages_of_gene_naive[:7] + coverages_of_gene[7:]
        coverage_plots[shape].plot(coverages_label, coverages_of_gene, shapes[0], label='Our Method')
        coverage_plots[shape].plot(coverages_label, coverages_of_gene_naive, shapes[1], label='Naive Method')
        coverage_plots[shape].spines['bottom'].set_color('black')
        coverage_plots[shape].spines['left'].set_color('black')
        coverage_plots[shape].legend(loc=5, fontsize='x-small')

        shape += 1


def plot_estimates(ru_estimate_plot, files):
    data = [[0 for _ in range(len(files) / 39)] for _ in range(39-2)]
    naive_data = [[0 for _ in range(len(files) / 39)] for _ in range(39-2)]
    # data = [[0 for _ in range(39)] for _ in range(len(files) / 39)]

    simulated = set([])
    for file_name in files:
        sim = int(file_name.split('_')[-2])
        simulated.add(sim)
    simulated = sorted(list(simulated))
    sim_to_ind = {sim: i for i, sim in enumerate(simulated)}
    # print sim_to_ind

    for file_name in files:
        coverage = int(file_name.split('_')[-1].split('.')[0][:-1])
        if coverage < 3:
            continue
        sim = int(file_name.split('_')[-2])
        estimate = 0
        naive_estimate = 0
        with open(file_name) as input:
            lines = input.readlines()
            if len(lines) > 1:
                if lines[-1].strip() != 'None' and len(lines[-1]) < 10:
                    estimate = int(float(lines[-1].strip()))
        with open(file_name + '.naive') as input:
            lines = input.readlines()
            if len(lines) > 1:
                if lines[-1].strip() != 'None' and len(lines[-1]) < 10:
                    naive_estimate = int(float(lines[-1].strip()))
        # print coverage, sim, sim_to_ind[sim]
        data[coverage - 3][sim_to_ind[sim]] = estimate
        naive_data[coverage - 3][sim_to_ind[sim]] = naive_estimate
        # data[sim_to_ind[sim]][coverage - 1] = estimate
    # if i == 0:
    #     ru_estimate_plots[i] = sns.tsplot(data, err_style="ci_bars")
    # else:
    from scipy import stats
    averages = []
    for sim_rcout in range(len(data[0])):
        total = 0
        for cov in range(len(data)):
            if data[cov][sim_rcout] == 0:
                print cov+9, sim_rcout
            total += data[cov][sim_rcout]
        averages.append(total / float(len(data)))

    naive_averages = []
    for sim_rcout in range(len(naive_data[0])):
        total = 0
        for cov in range(len(naive_data)):
            total += naive_data[cov][sim_rcout]
        naive_averages.append(total / float(len(naive_data)))

    ru_estimate_plot.errorbar(simulated, averages, yerr=stats.sem(data))
    ru_estimate_plot.set_xticks(simulated, [sim+1 for sim in simulated])
    ru_estimate_plot.set_yticks(simulated, [sim+1 for sim in simulated])
    # ru_estimate_plot.errorbar(simulated, naive_averages, yerr=stats.sem(naive_data))


def add_estimates_for_three_genes(ru_estimate_plots, results_dir):
    import glob
    import seaborn as sns
    gene_dirs = glob.glob(results_dir + '*')

    for i, gene_dir in enumerate(gene_dirs):
        gene_name = gene_dir.split('/')[-1]
        files = glob.glob(gene_dir + '/*.fasta.out')
        plot_estimates(ru_estimate_plots[i], files)
        # plot_estimates(ru_estimate_plots[i], gene_dir + '/*.fasta.out.naive')

        ru_estimate_plots[i].spines['bottom'].set_color('black')
        ru_estimate_plots[i].spines['left'].set_color('black')
        ru_estimate_plots[i].set_title(gene_name)

        # ru_estimate_plots[i].set_xlim([0.5, 60])
        # ru_estimate_plots[i].set(xlim=(-0.1, None))
        # ru_estimate_plots[i].xticks(x_axis, filtered_indels, fontsize=10)


def plot_pacbio_ru_results_for_three_genes(results_dir='../pacbio_coverage_experiment/'):
    from matplotlib import rc, rcParams
    import matplotlib.pyplot as plt
    plt.style.use('ggplot')
    plt.rcParams['axes.facecolor'] = '#FFFFFF'
    rc('text', usetex=True)
    rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']
    plt.title('Effect of Sequencing Coverage on Copy Number Estimation')
    plt.gca().spines['bottom'].set_color('black')
    plt.gca().spines['left'].set_color('black')

    fig = plt.figure()
    ax = list([])
    x_label_font = 15
    y_label_font = 13

    # ax.append(fig.add_subplot(311))
    # ax[0].set_ylabel(r'\emph{Recruited Reads}', fontsize=y_label_font)
    # ax[0].set_xlabel(r'\emph{Simulated RU Count}')

    ax.append(fig.add_subplot(211))
    ax[0].set_ylabel(r'\emph{Estimated RU Count}', fontsize=y_label_font)
    ax[0].set_xlabel(r'\emph{Simulated RU Count}', fontsize=x_label_font)

    ax.append(fig.add_subplot(212))
    ax[1].set_ylabel(r'\emph{Correct Estimates}', fontsize=y_label_font)
    ax[1].set_xlabel(r'\emph{Sequencing Coverage}', fontsize=x_label_font)
    # Turn off axis lines and ticks of the big subplot
    for i in range(2):
        ax[i].spines['top'].set_color('none')
        ax[i].spines['bottom'].set_color('none')
        ax[i].spines['left'].set_color('none')
        ax[i].spines['right'].set_color('none')
        ax[i].tick_params(labelcolor='w', top='off', bottom='off', left='off', right='off')

    plt.tight_layout(pad=2, w_pad=0.5, h_pad=1.0)

    ru_estimate_plots = list([])
    ru_estimate_plots.append(fig.add_subplot(231))
    ru_estimate_plots.append(fig.add_subplot(232))
    ru_estimate_plots.append(fig.add_subplot(233))
    add_estimates_for_three_genes(ru_estimate_plots, results_dir)

    coverage_plots = list([])
    coverage_plots.append(fig.add_subplot(234))
    coverage_plots.append(fig.add_subplot(235, sharey=coverage_plots[0]))
    coverage_plots.append(fig.add_subplot(236, sharey=coverage_plots[0]))
    add_coverages_for_three_genes(coverage_plots, results_dir)

    plt.savefig('pacbio_ru_count_results.png', dpi=300)


def plot_pacbio_copy_number_simulation_results(results_dir='../Pacbio_copy_number/'):
    #NEEDS VINEET's COMMENT BEFORE PROCEEDING
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import glob

    gene_dirs = glob.glob(results_dir + '*')
    for directory in gene_dirs:
        gene_name = directory.split('/')[-1]
        if gene_name != 'CSTB':
            continue
        files = glob.glob(directory + '/*.out')
        points = []
        for file_name in files:
            sim = int(file_name.split('_')[-2])
            if sim % 4 != 0:
                continue
            with open(file_name) as input:
                lines = input.readlines()
                estimate = int(lines[-1].strip())
            points.append((sim, estimate))
        points = sorted(points)
        sim_repeats = [sim for sim, estimate in points]
        estimated_repeats = [estimate for sim, estimate in points]
        plt.title('Result of estimation on PacBio simulated reads')
        plt.ylabel('Estimated Copy Number')
        plt.xlabel('Simulated Copy Number')

        ind = sim_repeats  # the x locations for the groups
        width = 1.0  # the width of the bars

        # fig = plt.figure()
        # ax = fig.add_subplot(111)

        yvals = sim_repeats
        rects1 = plt.bar(ind, yvals, width, label=gene_name, color=cm.Accent(float(0) / 3))
        zvals = sim_repeats
        rects2 = plt.bar([i + width for i in ind], zvals, width, color=cm.Accent(float(1) / 3))
        kvals = sim_repeats
        rects3 = plt.bar([i + width * 2 for i in ind], kvals, width, color=cm.Accent(float(2) / 3))

        # plt.xticks([i + width for i in ind])
        # plt.set_xticklabels(('2011-Jan-4', '2011-Jan-5', '2011-Jan-6'))
        plt.legend((rects1[0], rects2[0], rects3[0]), ('VNTR1', 'VNTR2', 'VNTR3'))

        def autolabel(rects):
            for rect in rects:
                continue
                h = rect.get_height()
                plt.text(rect.get_x() + rect.get_width() / 2., 1.05 * h, '%d' % int(h),
                        ha='center', va='bottom')

        autolabel(rects1)
        autolabel(rects2)
        autolabel(rects3)

        plt.savefig('Pacbio_copy_number_simulation_results.png', dpi=300)


def plot_false_read_and_true_read_score_distribution(results_dir='results/score_distribution/'):
    import glob
    import os
    false_scores_files = glob.glob(results_dir + 'false_scores*')
    for i in range(len(false_scores_files)):
        false_scores_file = false_scores_files[i]
        true_scores_file = results_dir + 'true_' + '_'.join(os.path.basename(false_scores_file).split('_')[1:])
        false_scores = get_numbers_from_file(false_scores_file)
        true_scores = get_numbers_from_file(true_scores_file)
        import matplotlib.pyplot as plt
        n_bins = 50
        labels = ['False Reads Scores', 'True Reads Scores']
        plt.hist([false_scores, true_scores], n_bins, histtype='bar', label=labels)
        plt.yscale('log')
        plt.title('Score distribution among mapped reads')
        plt.ylabel('Number of reads')
        plt.xlabel('Score')
        plt.legend(prop={'size': 10})
        vntr_id = os.path.basename(false_scores_file).split('_')[3]
        plt.savefig('score_distribution_%s.png' % vntr_id, dpi=300)
        plt.cla()
        plt.clf()
        plt.close()


def plot_pacbio_flanking_region_sizes():
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(1)
    sizes_1217 = [298, 333, 355, 310, 311, 313, 350, 339, 292, 309, 385, 329, 295, 301]
    X = [i for i, _ in enumerate(sizes_1217)]
    line_width = [0, len(X)]
    cmap = plt.get_cmap('jet', 20)
    for i in range(19, 12, -1):
        line_y = [i * 20, i * 20]
        c = cmap((i-13) * 2)
        ax.plot(line_width, line_y, '-', color=c)
        sigma1 = np.array([10, 10])
        ax.fill_between(line_width, line_y + sigma1, line_y - sigma1, facecolor=c, alpha=0.3, label='%s Repeats' % i)
    ax.plot(X, sizes_1217, 'o', label='distance of flanking regions', color='red')
    ax.set_xlabel('read index')
    ax.set_ylabel('flanking region distance (bp)')
    ax.legend(loc=0, fontsize = 'x-small')

    plt.show()

edges = [(1, 8), (1, 16), (2, 17), (4, 18), (8, 16), (30, 32), (30, 33), (32, 33), (34, 40), (34, 47), (38, 57),
         (38, 59), (38, 67), (40, 47), (57, 59), (57, 67), (59, 67)] + [(5, 53), (47, 19), (71, 3), (31, 9)]
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
# plot_FP_for_specific_sensitivity(eliminated_nodes)
# plot_gc_content_violin_plot()
# plot_paccbio_flanking_region_sizes()
# plot_frequency_of_repeats_in_population()

# plot_pacbio_ru_results_for_three_genes()
# plot_pacbio_ru_length_result()

# plot_read_selection_and_mapping_sensitivity_fdr_curve()

# plot_indel_frequencies_for_diabetes()

plot_read_selection_and_mapping_sensitivity_comparison()
