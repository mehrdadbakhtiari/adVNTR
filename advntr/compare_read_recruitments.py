import os
import glob
import pysam

from Bio import SeqIO

from reference_vntr import load_unique_vntrs_data
from sam_utils import get_id_of_reads_mapped_to_vntr_in_bamfile, make_bam_and_index
from vntr_finder import VNTRFinder


def count_reads(bam_file):
    print bam_file
    alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    return alignment_file.unmapped + alignment_file.mapped


def count_reads_mapped_to_vntr_in_bamfile(bam_file, reference_vntr):
    return len(get_id_of_reads_mapped_to_vntr_in_bamfile(bam_file, reference_vntr))


def get_pacbio_true_read_ids(bam_file, reference_vntr, ref_length):
    alignment_file = pysam.AlignmentFile(bam_file, 'rb')
    true_read_ids = []
    for read in alignment_file.fetch():
        read_end = read.reference_start + len(read.seq)
        if read.reference_start < ref_length - 1000 and read_end > 1000:
            true_read_ids.append(read.qname)
            continue
    return true_read_ids


def bowtie_alignment(fq_file):
    bowtie_alignment_file = fq_file[:-3] + '_bowtie_aln.sam'
    os.system('bowtie2 -x hg19_chromosomes/hg19_bt2_idx -U %s -S %s --threads 20' % (fq_file, bowtie_alignment_file))
    make_bam_and_index(bowtie_alignment_file)
    return bowtie_alignment_file[:-4] + '.bam'


def bwamem_alignment(fq_file):
    bwa_alignment_file = fq_file[:-3] + '_bwa_aln.sam'
    os.system('bwa mem -t 20 hg19_chromosomes/CombinedHG19_Reference.fa %s > %s' % (fq_file, bwa_alignment_file))
    make_bam_and_index(bwa_alignment_file)
    return bwa_alignment_file[:-4] + '.bam'


def bwasw_alignment(fq_file):
    bwa_alignment_file = fq_file[:-3] + '_bwasw_aln.sam'
    if not os.path.exists(bwa_alignment_file[:-4] + '.bam'):
        os.system('bwa bwasw -t 20 hg19_chromosomes/CombinedHG19_Reference.fa %s > %s' % (fq_file, bwa_alignment_file))
        make_bam_and_index(bwa_alignment_file)
    return bwa_alignment_file[:-4] + '.bam'


def blasr_alignment(fq_file):
    blasr_alignment_file = fq_file[:-3] + '_blasr_aln.sam'
    if not os.path.exists(blasr_alignment_file[:-4] + '.bam'):
        print('Running Blasr for %s' % blasr_alignment_file[:-4] + '.bam')
        os.system('blasr %s hg19_chromosomes/CombinedHG19_Reference.fa -sam -noSplitSubreads -out %s -nproc 8' % (fq_file, blasr_alignment_file))
        make_bam_and_index(blasr_alignment_file)
    return blasr_alignment_file[:-4] + '.bam'


def get_our_selected_reads_count(fq_file, vntr_finder):
    hmm = vntr_finder.get_vntr_matcher_hmm(read_length=150)
    min_score_to_count_read = vntr_finder.get_min_score_to_select_a_read(hmm, None, 150)
    selected_reads = []
    fasta_sequences = SeqIO.parse(open(fq_file), 'fastq')
    for fasta in fasta_sequences:
        logp, vpath = hmm.viterbi(fasta.seq)
        rev_logp, rev_vpath = hmm.viterbi(str(fasta.seq.reverse_complement()))
        logp = max(rev_logp, logp)
        if logp > min_score_to_count_read:
            selected_reads.append(fasta)
    return len(selected_reads)


def get_our_filtered_reads_count(fq_file, vntr_finder):
    unmapped_reads_file = fq_file[:-3] + 'fa'
    os.system("cat %s | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > %s" % (fq_file, unmapped_reads_file))
    filtered_ids = vntr_finder.filter_reads_with_keyword_matching('working_directory/', unmapped_reads_file)
    return len(filtered_ids)


def get_out_pacbio_filtered_counts(fq_file, vntr_finder):
    unmapped_reads_file = fq_file[:-3] + 'fa'
    # os.system("cat %s | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > %s" % (fq_file, unmapped_reads_file))
    os.system('../seqtk/seqtk seq -a %s > %s' % (fq_file, unmapped_reads_file))
    filtered_ids = vntr_finder.filter_reads_with_keyword_matching('working_directory/', unmapped_reads_file, False)
    return len(filtered_ids)

    unmapped_reads = SeqIO.parse(unmapped_reads_file, 'fasta')
    filtered_reads = []
    for read in unmapped_reads:
        for read.id in filtered_ids:
            filtered_reads.append(read)
    return len(vntr_finder.get_spanning_reads_of_unaligned_pacbio_reads(filtered_reads))


def get_pacbio_comparison_result():
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {1221: 'CSTB', 1216: 'HIC1', 1215: 'INS'}
    genes = glob.glob('../Pacbio_copy_number/*')
    for gene_dir in genes:
        print(gene_dir)
        files = glob.glob(gene_dir + '/*30x.fastq.sam')
        gene_name = gene_dir.split('/')[-1]
        mapped_reads = {}
        for file_name in files:
            copies = int(file_name.split('_')[-2])
            make_bam_and_index(file_name)
            base_name = file_name[:-4]
            original_bam = base_name + '.bam'
            bwasw_alignment(base_name)
            blasr_alignment(base_name)
            bwasw_alignment_file = base_name[:-3] + '_bwasw_aln.bam'
            blasr_alignment_file = base_name[:-3] + '_blasr_aln.bam'

            vntr_id = None
            for vid, gname in id_to_gene.items():
                if gname == gene_name:
                    vntr_id = vid

            ref_length = copies * len(reference_vntrs[vntr_id].pattern) + 2000
            true_ids = get_pacbio_true_read_ids(original_bam, reference_vntrs[vntr_id], ref_length)
            blasr_ids = get_id_of_reads_mapped_to_vntr_in_bamfile(blasr_alignment_file, reference_vntrs[vntr_id])
            bwasw_ids = get_id_of_reads_mapped_to_vntr_in_bamfile(bwasw_alignment_file, reference_vntrs[vntr_id])
            blasr_tp = [read_id for read_id in blasr_ids if read_id in true_ids]
            bwasw_tp = [read_id for read_id in bwasw_ids if read_id in true_ids]
            vntr_finder = VNTRFinder(reference_vntrs[vntr_id])
            our_filtering = get_out_pacbio_filtered_counts(base_name, vntr_finder)
            our_selection = our_filtering
            mapped_reads[copies] = [len(true_ids), our_filtering, our_selection, len(bwasw_tp), len(blasr_tp)]
        with open(gene_dir + '/result.txt', 'w') as out:
            for copies in sorted(mapped_reads.iterkeys()):
                original, our_filtering, our_selection, bwasw, blasr = mapped_reads[copies]
                out.write('%s %s %s %s %s %s\n' % (copies, original, our_filtering, our_selection, bwasw, blasr))


def get_illumina_comparison_result():
    reference_vntrs = load_unique_vntrs_data()
    id_to_gene = {119: 'DRD4', 1220: 'GP1BA', 1221: 'CSTB', 1214: 'MAOA', 1219: 'IL1RN'}
    genes = glob.glob('../Illumina_copy_number/*')
    for gene_dir in genes:
        print(gene_dir)
        files = glob.glob(gene_dir + '/*30x.sam')
        gene_name = gene_dir.split('/')[-1]
        print(len(files))
        mapped_reads = {}
        for file_name in files:
            copies = file_name.split('_')[-2]
            make_bam_and_index(file_name)
            base_name = file_name[:-4]
            bowtie_bam = bowtie_alignment(base_name + '.fq')
            bwa_bam = bwamem_alignment(base_name + '.fq')
            original_bam = file_name[:-4] + '.bam'

            vntr_id = None
            for vid, gname in id_to_gene.items():
                if gname == gene_name:
                    vntr_id = vid

            vntr_finder = VNTRFinder(reference_vntrs[vntr_id])
            original = count_reads(original_bam)
            our_selection = get_our_selected_reads_count(base_name + '.fq', vntr_finder)
            our_filtering = get_our_filtered_reads_count(base_name + '.fq', vntr_finder)
            bwa = count_reads_mapped_to_vntr_in_bamfile(bwa_bam, reference_vntrs[vntr_id])
            bowtie = count_reads_mapped_to_vntr_in_bamfile(bowtie_bam, reference_vntrs[vntr_id])
            mapped_reads[int(copies)] = [original, our_filtering, our_selection, bwa, bowtie]
        with open(gene_dir + '/result.txt', 'w') as out:
            for copies in sorted(mapped_reads.iterkeys()):
                original, our_filtering, our_selection, bwa, bowtie = mapped_reads[copies]
                out.write('%s %s %s %s %s %s\n' % (copies, original, our_filtering, our_selection, bwa, bowtie))

get_pacbio_comparison_result()
