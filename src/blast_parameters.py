

class ValidParameter():
    def __init__(self, reward, penalty, gapopen, gapextend):
        self.reward = reward
        self.penalty = penalty
        self.gapopen = gapopen
        self.gapextend = gapextend


def read_blast_params_list(file_name='blast_parameters.txt'):
    result = []
    with open(file_name) as input:
        lines = input.readlines()
        for line in lines:
            reward_penalty, gaps_list = line.strip().split('\t')
            reward, penalty = reward_penalty.split('/')
            gaps_list = gaps_list.split(',')
            gaps_list = [gap.split('/') for gap in gaps_list]
            for gapopen, gapextend in gaps_list:
                result.append(ValidParameter(reward, penalty, gapopen, gapextend))

    return result
