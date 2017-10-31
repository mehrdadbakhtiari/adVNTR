import os

from Bio import SearchIO
from Bio import pairwise2
from Bio import Seq, SeqRecord, SeqIO

from hmm_utils import build_reference_repeat_finder_hmm, get_repeat_segments_from_visited_states_and_region
from utils import *
import settings


class ReferenceVNTR:
    def __init__(self, vntr_id, pattern, start_point, chromosome, estimated_repeats=None):
        self.non_overlapping = True
        self.has_homologous = False
        self.id = vntr_id
        self.pattern = pattern
        self.start_point = start_point
        self.chromosome = chromosome
        self.estimated_repeats = estimated_repeats
        self.repeat_segments = None
        self.left_flanking_region = None
        self.right_flanking_region = None

    def init_from_vntrseek_data(self):
        corresponding_region_in_ref = self.get_corresponding_region_in_ref()
        repeat_segments = self.find_repeat_segments(corresponding_region_in_ref)
        self.repeat_segments = repeat_segments
        flanking_region_size = 500
        self.left_flanking_region, self.right_flanking_region = self.get_flanking_regions(flanking_region_size)

    def init_from_xml(self, repeat_segments, left_flanking_region, right_flanking_region):
        self.repeat_segments = repeat_segments
        self.left_flanking_region = left_flanking_region
        self.right_flanking_region = right_flanking_region

    def is_non_overlapping(self):
        return self.non_overlapping

    def has_homologous_vntr(self):
        return self.has_homologous

    def get_length(self):
        return sum([len(e) for e in self.repeat_segments])

    def get_repeat_segments(self):
        return self.repeat_segments

    def is_homologous_vntr(self, another):
        structure1 = self.left_flanking_region[-20:] + self.pattern + self.right_flanking_region[:20]
        structure2 = another.left_flanking_region[-20:] + another.pattern + another.right_flanking_region[:20]
        alignment_score = pairwise2.align.localms(structure1, structure2, 1, -1, -1, -1, score_only=True)
        if float(alignment_score) / len(structure1) > 0.66 or float(alignment_score) / len(structure2) > 0.66:
            return True
        return False

    def find_repeat_segments(self, region_in_ref):
        patterns = [self.pattern]
        model = build_reference_repeat_finder_hmm(patterns, copies=self.estimated_repeats)
        logp, path = model.viterbi(region_in_ref)
        visited_states = [state.name for idx, state in path[1:-1]]
        repeat_segments = get_repeat_segments_from_visited_states_and_region(visited_states, region_in_ref)

        return repeat_segments

    def __get_chromosome_reference_sequence(self):
        return get_chromosome_reference_sequence(self.chromosome)

    def get_corresponding_region_in_ref(self):
        ref_sequence = self.__get_chromosome_reference_sequence()
        estimated_length = len(self.pattern) * self.estimated_repeats
        corresponding_region_in_ref = ref_sequence[self.start_point:self.start_point + estimated_length].upper()
        return corresponding_region_in_ref

    def get_flanking_regions(self, flanking_region_size=140):
        ref_sequence = self.__get_chromosome_reference_sequence()
        left_flanking = ref_sequence[self.start_point - flanking_region_size:self.start_point].upper()
        end_of_repeats = self.start_point + self.get_length()
        right_flanking = ref_sequence[end_of_repeats:end_of_repeats + flanking_region_size].upper()
        return left_flanking, right_flanking


def load_unprocessed_vntrseek_data(vntrseek_output, chromosome=None):
    vntrs = []
    with open(vntrseek_output) as input_file:
        input_lines = [line.strip() for line in input_file.readlines() if line.strip() != '']
        for vntr_id, line in enumerate(input_lines):
            vntrseek_repeat, _, pattern, chromosome_number, start = line.split()
            estimated_repeats = int(float(vntrseek_repeat) + 5)
            if chromosome and chromosome_number != chromosome:
                continue
            vntrs.append(ReferenceVNTR(vntr_id, pattern, int(start)-1, chromosome_number, estimated_repeats))
    return vntrs


def find_non_overlapping_vntrs(vntrseek_output='repeats_length_patterns_chromosomes_starts.txt'):
    vntrs = load_unprocessed_vntrseek_data(vntrseek_output)
    skipped_vntrs = []
    for i in range(len(vntrs)):
        print(i)
        estimated_end = len(vntrs[i].pattern) * vntrs[i].estimated_repeats + vntrs[i].start_point
        if i < len(vntrs)-1 and vntrs[i].chromosome == vntrs[i+1].chromosome and estimated_end > vntrs[i+1].start_point:
            vntrs[i].estimated_repeats += vntrs[i+1].estimated_repeats
        vntrs[i].init_from_vntrseek_data()
        repeat_segments = vntrs[i].get_repeat_segments()
        if i in skipped_vntrs:
            vntrs[i].non_overlapping = False
        else:
            j = i + 1
            end_point = len(vntrs[i].pattern) * len(repeat_segments) + vntrs[i].start_point
            while j < len(vntrs) and vntrs[i].chromosome == vntrs[j].chromosome and end_point > vntrs[j].start_point:
                skipped_vntrs.append(j)
                j += 1
    return vntrs


def process_vntrseek_data(processed_vntrs_file='repeats_and_segments.txt'):
    vntrs = find_non_overlapping_vntrs()
    for vntr in vntrs:
        comma_separated_segments = ','.join(vntr.get_repeat_segments())
        with open(processed_vntrs_file, 'a') as out:
            out.write('%s %s %s %s %s\n' % (vntr.id, vntr.is_non_overlapping(), vntr.left_flanking_region,
                                            vntr.right_flanking_region, comma_separated_segments))


def identify_homologous_vntrs(vntrs, chromosome=None):
    for i in range(len(vntrs)):
        for j in range(i + 1, len(vntrs)):
            if chromosome and (chromosome != vntrs[i].chromosome and chromosome != vntrs[j].chromosome):
                continue
            if vntrs[i].is_homologous_vntr(vntrs[j]):
                vntrs[i].has_homologous = True
                vntrs[j].has_homologous = True
    return vntrs


def load_unique_vntrs_data(vntrseek_output='repeats_length_patterns_chromosomes_starts.txt'):
    vntrs = []
    with open(vntrseek_output) as input_file:
        vntrseek_data = [line.strip() for line in input_file.readlines() if line.strip() != '']
    with open('repeats_and_segments.txt') as input_file:
        segments_lines = input_file.readlines()
    for vntr_id in range(len(vntrseek_data)):
        _id, non_overlapping, left_flanking_region, right_flanking_region, segments = segments_lines[vntr_id].split()
        repeat_segments = segments.split(',')
        vntrseek_repeat, _, pattern, chromosome, start = vntrseek_data[vntr_id].split()
        vntr = ReferenceVNTR(vntr_id, pattern, int(start)-1, chromosome, vntrseek_repeat)
        vntr.init_from_xml(repeat_segments, left_flanking_region, right_flanking_region)
        vntr.non_overlapping = True if non_overlapping == 'True' else False
        vntrs.append(vntr)
    return vntrs


def is_false_vntr_hit(qresult, ref_vntr, position, threshold):
    for hit in qresult:
        for hsp in hit:
            score = hsp.match_num - hsp.mismatch_num - hsp.hit_gapopen_num
            if score > 65:
                if ref_vntr.chromosome == hit.id and abs(position - hsp.hit_start) <= threshold:
                    continue
                else:
                    print('found in ', hit.id, hsp.hit_start)
                    return True
    return False


def find_similar_region_for_vntr(sema, reference_vntr, vntr_id, result_list):
    vntr_len = reference_vntr.get_length()
    searches = list([])
    searches.append((reference_vntr.pattern * 2, reference_vntr.start_point, 2 * vntr_len + 30))
    searches.append((reference_vntr.left_flanking_region[-100:], reference_vntr.start_point, vntr_len + 20))
    searches.append((reference_vntr.right_flanking_region[:100], reference_vntr.start_point + vntr_len, vntr_len + 20))
    ref_file = 'hg19_chromosomes/CombinedHG19_Reference.fa'
    found = 0
    for search in searches:
        qfile = settings.BLAST_TMP_DIR + str(i) + '_query.fasta'
        with open(qfile, "w") as output_handle:
            my_rec = SeqRecord.SeqRecord(seq=Seq.Seq(search[0]), id='query', description='')
            SeqIO.write([my_rec], output_handle, 'fasta')
        output = 'output.psl'
        command = 'blat -q=dna -oneOff=1 -tileSize=8 -stepSize=6 -minIdentity=70 %s %s %s' % (ref_file, qfile, output)
        os.system(command)
        os.system('rm %s' % qfile)
        qresult = SearchIO.read('output.psl', 'blat-psl')
        if is_false_vntr_hit(qresult, reference_vntr, search[1], search[2]):
            found += 1
    if found > 1:
        print('there is similar sequence for %s' % vntr_id)
        result_list.append(vntr_id)
    sema.release()


def identify_similar_regions_for_vntrs_using_blat():
    from multiprocessing import Process, Semaphore, Manager

    reference_vntrs = load_unique_vntrs_data()
    sema = Semaphore(24)
    manager = Manager()
    result_list = manager.list()
    process_list = []
    for i in range(len(reference_vntrs)):
        if not reference_vntrs[i].is_non_overlapping() or reference_vntrs[i].has_homologous_vntr():
            continue
        sema.acquire()
        p = Process(target=find_similar_region_for_vntr, args=(sema, reference_vntrs[i], i, result_list))
        process_list.append(p)
        p.start()

    for p in process_list:
        p.join()
    result_list = list(result_list)
    with open('similar_vntrs.txt', 'w') as out:
        for vntr_id in result_list:
            out.write('%s\n' % vntr_id)


def extend_flanking_regions_in_processed_vntrs(flanking_size=500, output_file='repeats_and_segments2.txt'):
    vntrs = load_unique_vntrs_data()
    reference_genomes = {}
    for vntr in vntrs:
        comma_separated_segments = ','.join(vntr.get_repeat_segments())
        if vntr.chromosome not in reference_genomes.keys():
            reference_genomes[vntr.chromosome] = get_chromosome_reference_sequence(vntr.chromosome)
        start = vntr.start_point
        left_flanking_region = reference_genomes[vntr.chromosome][start-flanking_size:start].upper()
        end = vntr.start_point + vntr.get_length()
        right_flanking_region = reference_genomes[vntr.chromosome][end:end+flanking_size].upper()
        with open(output_file, 'a') as out:
            out.write('%s %s %s %s %s\n' % (vntr.id, vntr.is_non_overlapping(), left_flanking_region,
                                            right_flanking_region, comma_separated_segments))

# process_vntrseek_data()
identify_similar_regions_for_vntrs_using_blat()
