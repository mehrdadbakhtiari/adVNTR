<!-- [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/advntr/README.html) -->
<!-- [![Anaconda-Server Badge](https://anaconda.org/bioconda/advntr/badges/downloads.svg)](https://anaconda.org/bioconda/advntr) -->
[![Documentation Status](https://readthedocs.org/projects/pip/badge/?version=stable)](http://pip.pypa.io/en/stable/?badge=stable)

code-adVNTR - A tool for genotyping coding VNTRs
------------------------------------
[code-adVNTR](https://github.com/mehrdadbakhtiari/adVNTR/tree/enhanced_hmm) is a tool for genotyping Variable Number Tandem Repeats (VNTR)
from sequence data. 

code-adVNTR utilizes multiple motif HMMs to identify small variants within motifs and estimate diploid repeat counts for VNTRs. 
It takes short reads, and a pre-trained HMM model for reference VNTRs as input and outputs either estimated diploid repeat count or small variants within the target VNTRs. 

Installation
------------
Currently, you can install code-adVNTR from source [install the adVNTR from source](http://advntr.readthedocs.io/en/latest/installation.html#install-from-source-not-recommended) with dependencies.

code-adVNTR could be invoked from command line with ``advntr``. In the future, code-adVNTR will be merged into the master branch of adVNTR for easy access.

Here are the instructions to install from source:

0) Prepare a directory that you want to download code-adVNTR and go to the directory, `cd "directory_name"`
1) Download the source code using the following command, `git clone https://github.com/mehrdadbakhtiari/adVNTR.git --branch enhanced_hmm`
2) Go to the directory that you downloaded the source code, `cd adVNTR`
3) Install adVNTR using the following command, `python setup install`
4) Download the reference VNTRs from this link below.
5) Run adVNTR for a VNTR using the following command, `advntr genotype -fs --vntr_id [id] --alignment_file [bam_file] -m [reference_vntr] --working_directory [working_dir]`
-fs is a parameter for variant detection, --vntr_id option is for specifying the target VNTR. For example **25561** for MUC1 VNTR with [hg19 reference VNTRs](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip).

Data Requirements
-----------------
In order to genotype VNTRs, you need to either train models for loci of interest or use pre-trained models (recommended):
* To run adVNTR on trained VNTR models:
    - Download [vntr_data_recommended_loci.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_recommended_loci.zip)
    and extract it inside the project directory. This includes a set of pre-trained VNTR models for Illumina (6719 loci)
    and Pacbio (8960 loci) sequencing data.
    - You can also download and use [vntr_data_genic_loci.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip)
    for 158522 VNTRs that results in having much longer running time.

Alternatively, you can add model for custom VNTR. See [Add Custom VNTR](http://advntr.readthedocs.io/en/latest/tutorial.html#add-custom-vntr-label) for more information about training models for custom VNTRs.

Execution:
----------
Use following command to see the help for running the tool.

    advntr --help

The program outputs either the RU count genotypes or small variants of trained VNTRs. To specify a single VNTR by its ID use ``--vntr_id <id>`` option.
The list of some known VNTRs and their ID is available at [Disease-linked-VNTRs page](https://github.com/mehrdadbakhtiari/adVNTR/wiki/Disease-linked-VNTRs) in wiki.

See the demo below or [Quickstart](http://advntr.readthedocs.io/en/latest/quickstart.html) page to see an example
data set with step-by-step genotyping commands.

Demo input in BAM format
------------------------
* ``--alignment_file`` specifies the alignment file containing mapped and unmapped reads:

```sh
    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/
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
Jonghun Park, Mehrdad Bakhtiari, Bernt Popp, Michael Wiesener, Vineet Bafna.
    <b>[Detecting tandem repeat variants in coding regions using code-adVNTR](https://doi.org/10.1016/j.isci.2022.104785) </b> iScience vol. 25,8 104785. (2022)
