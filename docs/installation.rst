Installation
============
In order to use adVNTR, it is recommended to (1) install adVNTR using conda packaging manager and (2) download the
predefined models for human genome from :ref:`data-requirements` section. However, you can install it from the source
and/or use custom models.

.. _install-with-conda:

Install adVNTR with conda
-------------------------
If you are using the conda packaging manager (e.g. `miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ or anaconda),
you can install adVNTR from the `bioconda  channel <https://bioconda.github.io/>`_:

::

    conda config --add channels bioconda
    conda install advntr

adVNTR could be invoked from command line with ``advntr``


.. _data-requirements:

Data Requirements
-----------------
In order to genotype VNTRs, you need to either train models for loci of interest or use pre-trained models (recommended):
* To run adVNTR on trained VNTR models:
    - Download `vntr_data_recommended_loci.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_recommended_loci.zip>`_ and extract it inside the project directory. This includes a set of pre-trained VNTR models in hg19 for Illumina (6719 loci) and Pacbio (8960 loci) sequencing data. You can use `vntr_data_recommended_loci_hg38.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_recommended_loci_hg38.zip>` for VNTRs in GRCh38.
    - You can also download and use `vntr_data_genic_loci.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data_genic_loci.zip>`_ for 158522 VNTRs in hg19 that results in having much longer running time.
Alternatively, you can add model for custom VNTR. See :ref:`add-custom-vntr-label` for more information about training models for custom VNTRs.

Execution:
----------
Use following command to see the help for running the tool.

::

    advntr --help

The program outputs the RU count genotypes of VNTRs. To specify a single VNTR by its ID use ``--vntr_id <id>`` option.
The list of some known VNTRs and their ID is available at `Disease-linked-VNTRs page <https://github.com/mehrdadbakhtiari/adVNTR/wiki/Disease-linked-VNTRs>`_ in wiki.

See the demo execution here or :ref:`quickstart` page to see an example data set with step-by-step genotyping commands.

Demo: input in BAM format
-------------------------
* ``--alignment_file`` specifies the alignment file containing mapped and unmapped reads:

::

    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/

* With ``--pacbio``, adVNTR assumes the alignment file contains PacBio sequencing data:

::

    advntr genotype --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio

* Use ``--frameshift`` to find the possible frameshifts in VNTR:

::

    advntr genotype --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift


Install from source (Not recommended)
-------------------------------------
You can also install adVNTR from the source instead of using conda. First, you need to install the following packages:

Dependencies
^^^^^^^^^^^^
1. Following libraries are required to be installed on the system:
    -   ``python2.7``
    -   ``python-pip``
    -   ``python-tk``
    -   ``libz-dev``
    -   ``samtools``
    -   ``muscle``

You can install these requirement in Ubuntu Linux by running ``sudo apt-get install python2.7 python-pip python-tk libz-dev samtools muscle``
To install these packages on Mac OS, you can use `Homebrew <https://brew.sh/>`_.

2. Following python2.7 packages are required:
    -   ``biopython``
    -   ``pysam`` version 0.9.1.4 or above
    -   ``cython``
    -   ``networkx`` version 1.11
    -   ``scipy``
    -   ``joblib``
    -   ``scikit-learn``

You can install required python libraries by running ``pip install -r requirements.txt``

To Install
^^^^^^^^^^
To get the latest version and install it locally, run:

::

    git clone https://github.com/mehrdadbakhtiari/adVNTR
    cd adVNTR
    make; make install
    python setup.py install

adVNTR could be invoked from command line with ``advntr``
