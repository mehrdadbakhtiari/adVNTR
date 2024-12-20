[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/advntr/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/advntr/badges/downloads.svg)](https://anaconda.org/bioconda/advntr)
[![Documentation Status](https://readthedocs.org/projects/pip/badge/?version=stable)](http://pip.pypa.io/en/stable/?badge=stable)

adVNTR - A tool for genotyping VNTRs
------------------------------------
[adVNTR](https://github.com/mehrdadbakhtiari/adVNTR/) is a tool for genotyping Variable Number Tandem Repeats (VNTR)
from sequence data. It works with both NGS short reads (Illumina HiSeq) and SMRT reads (PacBio) and finds
diploid repeating counts for VNTRs and identifies possible mutations in the VNTR sequences.

[code-adVNTR](https://github.com/mehrdadbakhtiari/adVNTR/tree/enhanced_hmm), a tool specialized in detecting small indel variants within motifs using short reads is now available. 
This tool employs multiple motif Hidden Markov Models (HMMs) to identify small variants within motifs and estimate diploid repeat counts for VNTRs specifically in coding regions. For more details, please refer to this [readme](https://github.com/mehrdadbakhtiari/adVNTR/tree/enhanced_hmm).

Installation
------------
If you are using the conda packaging manager (e.g. [miniconda](https://docs.conda.io/en/latest/miniconda.html) or anaconda),
you can install adVNTR from the [bioconda  channel](https://bioconda.github.io/):

    conda config --add channels bioconda
    conda install -c conda-forge -c bioconda advntr

adVNTR could be invoked from command line with ``advntr``

Alternatively, you can install dependencies and [install the adVNTR from source](http://advntr.readthedocs.io/en/latest/installation.html#install-from-source-not-recommended).


Data Requirements and Pre-trained Models (Databases)
-----------------
In order to genotype VNTRs, you need to either train models for loci of interest or use pre-trained models (recommended):
* To run adVNTR on trained VNTR models:
    - Download [vntr_data_recommended_loci.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_recommended_loci.zip)
    and extract it inside the project directory. This includes a set of pre-trained VNTR models in hg19 for Illumina (6719 loci)
    and Pacbio (8960 loci) sequencing data. You can use [vntr_data_recommended_loci_hg38.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_recommended_loci_hg38.zip) for VNTRs in GRCh38.
    - You can also download and use [vntr_data_genic_loci.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip)
    for 158522 VNTRs in hg19 that results in having much longer running time.
    - You can download the database for the gene-proximal and phenotype associated VNTRs [G-VNTRs and P-VNTRs database](https://drive.google.com/file/d/1rF1CIliwzFcJmrCU2ibMVVZncsU6uNcs/view?usp=share_link) based on the GRCh38 reference.
    - Alternatively you can download the BED file for gene-proximal VNTRs [gene_proximal_vntrs.bed](https://drive.google.com/file/d/1DetpBQySPNe2YAJa4FsjHn9qiRNS3wEV/view?usp=share_link) for VNTRs in GRCh38.
    - The VNTR coordinates (GRCh38) utilized for GIAB TR benchmarking can be downloaded from [this link](https://drive.google.com/file/d/1FCk4HeXQMBiwnXeLucmIgSrz2qcmEcts/view?usp=drive_link).

Alternatively, you can add model for custom VNTR. See [Add Custom VNTR](http://advntr.readthedocs.io/en/latest/tutorial.html#add-custom-vntr-label) for more information about training models for custom VNTRs.

[Optional] For faster genotyping with adVNTR-NN, pretrained neural network models can be downloaded from [here](https://drive.google.com/drive/folders/1xeIoaE_iX4JojfKjlUkqXQ0iONPR5Zax?usp=sharing).

Execution:
----------
Use following command to see the help for running the tool.

    advntr --help

The program outputs the RU count genotypes of trained VNTRs. To specify a single VNTR by its ID use ``--vntr_id <id>`` option.
The list of some known VNTRs and their ID is available at [Disease-linked-VNTRs page](https://github.com/mehrdadbakhtiari/adVNTR/wiki/Disease-linked-VNTRs) in wiki.

See the demo below or [Quickstart](http://advntr.readthedocs.io/en/latest/quickstart.html) page to see an example
data set with step-by-step genotyping commands.

Demo input in BAM format
------------------------
* ``--alignment_file`` specifies the alignment file containing mapped and unmapped reads:

```sh
    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/
```

* With ``--pacbio``, adVNTR assumes the alignment file contains PacBio sequencing data:

```sh
    advntr genotype --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio
```

* Use ``--frameshift`` to find the possible frameshifts in VNTR:

```sh
    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift
```

Documentation:
--------------
Documentation is available at [advntr.readthedocs.io](http://advntr.readthedocs.io).

See [Quickstart](http://advntr.readthedocs.io/en/latest/quickstart.html) page to see an example data set with step-by-step genotyping commands.

Citation:
---------
- code-adVNTR:

    Jonghun Park, Mehrdad Bakhtiari, Bernt Popp, Michael Wiesener, Vineet Bafna.
    <b>[Detecting tandem repeat variants in coding regions using code-adVNTR](https://doi.org/10.1016/j.isci.2022.104785) </b> iScience vol. 25,8 104785. (2022),

- adVNTR-NN (v1.4.0):

    Mehrdad Bakhtiari, Jonghun Park, Yuan-Chun Ding, Sharona Shleizer-Burko, Susan L. Neuhausen, Bjarni V. Halldorsson, 
    Kari Stefansson, Melissa Gymrek, Vineet Bafna. <b>[Variable Number Tandem Repeats mediate the expression of proximal genes](https://doi.org/10.1038/s41467-021-22206-z) </b> 
    Nature Communications 12, 2075 (2021).

- Original publication (adVNTR):

    Bakhtiari, M., Shleizer-Burko, S., Gymrek, M., Bansal, V. and Bafna, V., 2018. 
    <b>[Targeted genotyping of variable number tandem repeats with adVNTR](https://genome.cshlp.org/content/28/11/1709). </b> Genome research vol. 28,11 (2018)
