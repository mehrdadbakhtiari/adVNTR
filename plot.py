
stat_files = ['0_size_related_reads.txt', '1_size_sensibility.txt', '2_size_blast_selected.txt',
              '3_sim_read_coverage__gc_content.txt']

x_label = {0: 'Pattern Size', 1: 'Pattern Size', 2: 'Pattern Size', 3: 'Simulated Read Coverage'}
y_label = {0: 'Reads from the pattern', 1: 'Sensibility', 2: 'Number of Selected Reads by Blast', 3: 'GC Content'}


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
