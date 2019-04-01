[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/advntr/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/advntr/badges/downloads.svg)](https://anaconda.org/bioconda/advntr)
[![Documentation Status](https://readthedocs.org/projects/pip/badge/?version=stable)](http://pip.pypa.io/en/stable/?badge=stable)

adVNTR - A tool for genotyping VNTRs
------------------------------------
[adVNTR](https://github.com/mehrdadbakhtiari/adVNTR/) is a tool for genotyping Variable Number Tandem Repeats (VNTR)
from sequence data. It works with both NGS short reads (Illumina HiSeq) and SMRT reads (PacBio) and finds
diploid repeating counts for VNTRs and identifies possible mutations in the VNTR sequences.

Installation
------------
If you are using the conda packaging manager (e.g. [miniconda](https://docs.conda.io/en/latest/miniconda.html) or anaconda),
you can install adVNTR from the [bioconda  channel](https://bioconda.github.io/):

    conda config --add channels bioconda
    conda install advntr

adVNTR could be invoked from command line with ``advntr``

Alternatively, you can install dependencies and [install the adVNTR from source](http://advntr.readthedocs.io/en/latest/installation.html#install-from-source-not-recommended).


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
Bakhtiari, M., Shleizer-Burko, S., Gymrek, M., Bansal, V. and Bafna, V., 2018. [Targeted genotyping of variable number tandem repeats with adVNTR](https://genome.cshlp.org/content/28/11/1709). Genome Research, 28(11), pp.1709-1719.
