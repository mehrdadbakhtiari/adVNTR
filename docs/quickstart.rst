Quick Start
===========
To help verify the installation and showing the workflow, we include a small data set and commands to genotype
this simulated dataset. If you have already installed adVNTR, jump to :ref:`genotype-simulated-dataset`.

Install
-------
The easiest way to get started is to :ref:`Install adVNTR with conda <install-with-conda>`.
To install adVNTR, run these commands:

::

    conda config --add channels bioconda
    conda install advntr


.. _genotype-simulated-dataset:

Genotype Predefined VNTR in Simulated Data
------------------------------------------
To genotype a VNTR in the simulated dataset, one option is to use predefined models.
Download `vntr_data.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data.zip>`_ and extract it inside the
project directory to use these models from human genome. Here, we genotype a VNTR with id 301645 that corresponds to a
disease-linked VNTR. The list of some known VNTRs and their ID is available at
`Disease-linked-VNTRs page <https://github.com/mehrdadbakhtiari/adVNTR/wiki/Disease-linked-VNTRs>`_ in wiki.

Then, download `simulated sequencing data of a human sample <https://conda.io/miniconda.html>`_.
It only includes reads around a VNTR in CSTB gene which is known to be linked to progressive myoclonus epilepsies.
Run this command to get 2/5 genotype for this VNTR.

::

    advntr genotype --vntr_id 301645 --alignment_file CSTB_2_5_testdata.bam --working_directory working_dir

Genotype Custom VNTR
--------------------
You can train a new model for a VNTR that doesn't exist in predefined models. Instead of downloading
`vntr_data.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data.zip>`_, you need the organism (here, human) reference
genome to train a model for a specific VNTR.
`Download chromosome 21 of hg19 <http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz>`_ and extract it.
It is recommended to have full reference genome of the organism to add the model, however, we use a single chromosome
in quickstart since it is easier to download and runs faster.
Run this command to add the VNTR in CSTB gene and train VNTR-specific scores:

::

    advntr addmodel -r chr21.fa -p CGCGGGGCGGGG -s 45196324 -e 45196360 -c chr21

If you run the above command without using predefined models, this VNTR gets the first id.
Run ``genotype`` command to get 2/5 genotype:

::

    advntr genotype --vntr_id 1 --alignment_file CSTB_2_5_testdata.bam --working_directory working_dir

