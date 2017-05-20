from settings import *
from utils import get_chromosome_reference_sequence, get_gc_content


class CoverageBiasDetector:
    """Find the coverage distribution based on GC content."""

    def __init__(self, sam_file='original_reads/paired_dat.sam'):
        pass

    def get_gc_contents_of_reference_windows(self):
        gc_map = {}
        for chromosome in CHROMOSOMES:
            gc_map[chromosome] = {}
            chromosome_reference = get_chromosome_reference_sequence(chromosome)
            for window_number in range(len(chromosome_reference) / GC_CONTENT_WINDOW_SIZE):
                start = window_number * GC_CONTENT_WINDOW_SIZE
                end = start + GC_CONTENT_WINDOW_SIZE
                gc_content = get_gc_content(chromosome_reference[start:end])
                gc_map[chromosome][window_number] = gc_content
        return gc_map

    def get_gc_content_coverage_map(self):
        reference_gc_map = self.get_gc_contents_of_reference_windows()


class CoverageCorrector:
    """Normalize the coverage based on the coverage distribution."""
    pass
