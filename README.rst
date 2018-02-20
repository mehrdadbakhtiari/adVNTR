adVNTR - A tool for genotyping VNTRs
------------------------------------
.. include:: docs/description.rst


Software Requirements
---------------------
1. Following libraries are required
    -   ``python2.7``
    -   ``python-pip``
    -   ``python-tk``
    -   ``libz-dev``
    -   ``samtools``

You can install these requirement in Ubuntu Linux by running ``sudo apt-get install python2.7 python-pip python-tk libz-dev samtools``

2. Following python2.7 packages are required:
    -   ``biopython``
    -   ``pysam`` version 0.9.1.4 or above
    -   ``cython``
    -   ``networkx`` version 1.11
    -   ``scipy``
    -   ``joblib``

You can install required python libraries by running ``pip install -r requirements.txt``

3. In addition, ``ncbi-blast`` version 2.2.29 or above is required


Data Requirements
-----------------
* To run adVNTR on trained VNTR models:
    - Download `vntr_data.zip <https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data.zip/>`_ and extract it inside the project directory.
Alternatively, you can add model for custom VNTR. See :ref:`add-custom-vntr-label` for more information.

Execution:
----------
Use following command to see the help for running the tool.

.. code:: bash
    
    python advntr.py --help

The program outputs the RU count genotypes for all VNTRs in ``vntr_data`` directory. To specify a single VNTR by its ID use ``--vntr_id <id>`` option. 

Demo 1: input in BAM format
---------------------------
* ``--alignment_file`` specifies the alignment file containing mapped and unmapped reads:

.. code:: bash
    
    python advntr.py --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/

* With ``--pacbio``, adVNTR assumes the alignment file contains PacBio sequencing data:

.. code:: bash
    
    python advntr.py --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio

* Use ``--frameshift`` to find the possible frameshifts in VNTR:

.. code:: bash
    
    python advntr.py --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift

Demo 2: input in fasta format
-----------------------------
* Use the following command to genotype the RU count using fasta file:

.. code:: bash
    
    python advntr.py --fasta unaligned_illumina_reads.fasta --working_directory ./log_dir/

Citation:
---------
.. include:: docs/publication.rst
