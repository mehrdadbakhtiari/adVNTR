from settings import *
from utils import get_chromosome_reference_sequence, get_gc_content
import pysam


class CoverageBiasDetector:
    """Find the coverage distribution based on GC content."""

    def __init__(self, alignment_file='original_reads/paired_dat.sam', chromosome=None, alignment_reference='HG19'):
        """
        :param alignment_file: The alignment file. It must be SAM/BAM file.
        :param chromosome: It specifies the chromosome in which the bias should be detected. ex: chr10.
            Defaults to None which corresponds to all chromosomes.
        :param alignment_reference: The alignment reference used in alignment file. Defaults to HG19.
        """
        self.alignment_file = alignment_file
        self.chromosome = chromosome
        self.alignment_reference = alignment_reference

    def get_gc_contents_of_reference_windows(self):
        gc_map = {}
        for chromosome in self.chromosome or CHROMOSOMES:
            gc_map[chromosome] = {}
            chromosome_reference = get_chromosome_reference_sequence(chromosome)
            for window_number in range(len(chromosome_reference) / GC_CONTENT_WINDOW_SIZE):
                start = window_number * GC_CONTENT_WINDOW_SIZE
                end = start + GC_CONTENT_WINDOW_SIZE
                gc_content = get_gc_content(chromosome_reference[start:end])
                gc_map[chromosome][window_number] = gc_content
        return gc_map

    def __add_bp_to_coverage_map(self, covered_bps, chromosome, window_number, read_start, read_end):
        start = max(window_number * GC_CONTENT_WINDOW_SIZE, read_start)
        end = min(window_number * GC_CONTENT_WINDOW_SIZE + GC_CONTENT_WINDOW_SIZE, read_end)
        if window_number not in covered_bps[chromosome]:
            covered_bps[chromosome][window_number] = 0
        covered_bps[chromosome][window_number] += end - start
        if read_end > window_number * GC_CONTENT_WINDOW_SIZE + GC_CONTENT_WINDOW_SIZE:
            self.__add_bp_to_coverage_map(covered_bps, chromosome, window_number+1, end, read_end)

    def get_covered_base_pairs_of_reference_windows(self):
        read_mode = 'r' if self.alignment_file.endswith('sam') else 'rb'
        samfile = pysam.AlignmentFile(self.alignment_file, read_mode)
        covered_bps = {}
        for chromosome in self.chromosome or CHROMOSOMES:
            covered_bps[chromosome] = {}

        reference = self.chromosome
        if self.chromosome and self.alignment_reference != 'HG19':
            reference = self.chromosome[3:]

        for read in samfile.fetch(reference, until_eof=True):
            window_number = read.reference_start / GC_CONTENT_WINDOW_SIZE
            read_start = read.reference_start
            read_end = read.reference_end if read.reference_end else read_start + len(read.seq)
            reference_name = read.reference_name
            if not reference_name.startswith('chr'):
                reference_name = 'chr' + reference_name
            if read.reference_name in covered_bps.keys():
                self.__add_bp_to_coverage_map(covered_bps, reference_name, window_number, read_start, read_end)
        return covered_bps

    def get_gc_content_coverage_map(self):
        reference_gc_map = self.get_gc_contents_of_reference_windows()
        covered_bps = self.get_covered_base_pairs_of_reference_windows()

        gc_coverage_map = {}
        for chromosome in covered_bps.keys():
            for window_number in covered_bps[chromosome]:
                windows_gc = reference_gc_map[chromosome][window_number]
                windows_gc = int(windows_gc * 10)
                windows_coverage = covered_bps[chromosome][window_number] / GC_CONTENT_WINDOW_SIZE
                if windows_coverage > OUTLIER_COVERAGE:
                    continue
                if windows_gc not in gc_coverage_map:
                    gc_coverage_map[windows_gc] = []
                gc_coverage_map[windows_gc].append(windows_coverage)

        return gc_coverage_map


class CoverageCorrector:
    """Normalize the coverage based on the coverage distribution."""

    def __init__(self, gc_coverage_map):
        """
        :param gc_coverage_map: The map of coverage distribution for each GC content,
            which is found by CoverageBiasDetector
        """
        self.gc_coverage_map = gc_coverage_map

    def get_sequencing_mean_coverage(self):
        windows_coverages = []
        for gc_content, coverages in self.gc_coverage_map.items():
            windows_coverages.append(coverages)
        return sum(windows_coverages) / float(len(windows_coverages))

    def get_mean_coverage_of_gc_content(self, gc_content):
        return sum(self.gc_coverage_map[gc_content]) / float(len(self.gc_coverage_map[gc_content]))

    def get_scaled_coverage(self, reference_vntr, observed_coverage):
        gc_content = get_gc_content(''.join(reference_vntr.get_repeat_segments()))
        scale_ratio = self.get_sequencing_mean_coverage() / self.get_mean_coverage_of_gc_content(gc_content)
        return observed_coverage * scale_ratio
