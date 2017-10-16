import ast
import glob


def find_vntrs():
    files = glob.glob('extend*')
    vntr_map = {}
    for input_file in files:
        with open(input_file) as infile:
            lines = infile.readlines()
            for i in range(len(lines)):
                if i % 2 == 1:
                    continue
                vntr_id = int(lines[i].strip())
                if lines[i+1][0] == '[':
                    copy_numbers = ast.literal_eval(lines[i+1].strip())
                elif lines[i+1].strip() != 'None':
                    copy_numbers = [int(lines[i+1].strip())] * 2
                else:
                    continue
                if vntr_id not in vntr_map.keys():
                    vntr_map[vntr_id] = set([])
                vntr_map[vntr_id] |= set(copy_numbers)

    result = []
    for vntr_id in vntr_map.keys():
        if len(vntr_map[vntr_id]) > 1:
            result.append(vntr_id)
            # print('%s: %s' % (vntr_id, vntr_map[vntr_id]))
    return result


def get_genotypes(file_name):
    genotypes = {} # genotypes[1214] = (31, 30)
    with open(file_name) as input:
        lines = input.readlines()
    for i in range(len(lines)):
        if i % 2 == 1:
            continue
        vntr_id = int(lines[i].strip())
        if lines[i + 1][0] == '[':
            copy_numbers = ast.literal_eval(lines[i + 1].strip())
        elif lines[i + 1].strip() != 'None':
            copy_numbers = [int(lines[i + 1].strip())]
        else:
            continue
        if len(copy_numbers) < 2:
            copy_numbers *= 2
        genotypes[vntr_id] = copy_numbers
    return genotypes


def inherited(allele, parent_genotype):
    min_diff = 1e9
    for parent_allele in parent_genotype:
        min_diff = min(min_diff, abs(parent_allele - allele) / float(allele))
    return min_diff < 1.0 / 15


def is_consistent(vntr_id, father_genotypes, mother_genotypes, child_genotypes):
    if vntr_id not in child_genotypes.keys():
        return True
    if vntr_id not in father_genotypes.keys() and vntr_id not in mother_genotypes.keys():
        return True

    child = child_genotypes[vntr_id]
    if vntr_id not in father_genotypes.keys():
        for allele in child:
            if allele in mother_genotypes[vntr_id]:
                return True
        return False

    if vntr_id not in mother_genotypes.keys():
        for allele in child:
            if allele in father_genotypes[vntr_id]:
                return True
        return False

    father = father_genotypes[vntr_id]
    mother = mother_genotypes[vntr_id]
    if inherited(child[0], father) and inherited(child[1], mother):
        return True
    if inherited(child[1], father) and inherited(child[0], mother):
        return True

    return False


def check_trio_consistency(father_file, mother_file, child_file):
    father_genotypes = get_genotypes(father_file)
    mother_genotypes = get_genotypes(mother_file)
    child_genotypes = get_genotypes(child_file)

    vntr_ids = set(father_genotypes.keys() + mother_genotypes.keys() + child_genotypes.keys())
    print('Total vntrs: %s' % len(vntr_ids))
    inconsistents = []
    for vid in vntr_ids:
        if not is_consistent(vid, father_genotypes, mother_genotypes, child_genotypes):
            inconsistents.append(vid)
            print('%s: %s %s %s' % (vid, father_genotypes[vid], mother_genotypes[vid], child_genotypes[vid]))

    print('Total inconsistencies: %s' % len(inconsistents))
