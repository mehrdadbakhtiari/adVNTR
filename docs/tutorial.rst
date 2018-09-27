Tutorial
========

Inputs
------
- NGS short reads (Illumina HiSeq)
- SMRT reads (PacBio)

Outputs
-------
- Text
- VCF (Under construction)

Usage
-----
adVNTR runs as follows:
::

    usage: advntr <command> [options]

There are four commands:

genotype
  Determine repeat unit count and sequence variation within VNTRs.

viewmodel
  Show the structure and information about the VNTRs in the database.

addmodel
  Add a custom VNTR to the database.

delmodel
  Delete a VNTR from the database.

Each of these commands and their options is described below.

Genotype
^^^^^^^^
Use :code:`advntr genotype [options]` to genotype a VNTR using sequencing data.
Alignment file and working directory are required.

Summary of options:

:code:`--frameshift`: Use this option to identify frameshift instead of finding copy number of a VNTR.

:code:`--pacbio`: Use this flag to genotype VNTRs using PacBio sequencing data.

:code:`--update`: Use this option to iteratively update the model using real data before finding the genotype.

Input/output options:

+---------------------------+--------------------------------------------------------------------------------+
| -a/--alignment_file <file>| Alignment file in BAM format or SAM format                                     |
+---------------------------+--------------------------------------------------------------------------------+
| -f/--fasta <file>         | Fasta file containing raw reads                                                |
+---------------------------+--------------------------------------------------------------------------------+
| -p/--pacbio               | set this flag if input file contains PacBio reads instead of Illumina reads    |
+---------------------------+--------------------------------------------------------------------------------+
| -n/--nanopore             | set this flag if input file contains Nanopore MinION reads instead of Illumina |
+---------------------------+--------------------------------------------------------------------------------+

Algorithm options:

+---------------------------+--------------------------------------------------------------------------------+
| -fs/--frameshift          | set this flag to search for frameshifts in VNTR instead of copy number.        |
+---------------------------+--------------------------------------------------------------------------------+
| -e/--expansion            | set this flag to determine long expansion from PCR-free data                   |
+---------------------------+--------------------------------------------------------------------------------+
| -c/--coverage <float>     | average sequencing coverage in PCR-free sequencing                             |
+---------------------------+--------------------------------------------------------------------------------+
| --haploid                 | set this flag if the organism is haploid                                       |
+---------------------------+--------------------------------------------------------------------------------+
| -naive/--naive            | use naive approach for PacBio reads                                            |
+---------------------------+--------------------------------------------------------------------------------+

Other options:

+---------------------------+--------------------------------------------------------------------------------+
| -h/--help                 | show this help message and exit                                                |
+---------------------------+--------------------------------------------------------------------------------+
| --working_directory <path>| working directory for creating temporary files needed for computation          |
+---------------------------+--------------------------------------------------------------------------------+
| -m/--models <file>        | file containing VNTRs information [vntr_data/hg19_VNTRs.db]                    |
+---------------------------+--------------------------------------------------------------------------------+
| -t/--threads <int>        | number of threads [4]                                                          |
+---------------------------+--------------------------------------------------------------------------------+
| -u/--update               | set this flag to iteratively update the model                                  |
+---------------------------+--------------------------------------------------------------------------------+
| -vid/--vntr_id <text>     | comma-separated list of VNTR IDs                                               |
+---------------------------+--------------------------------------------------------------------------------+


View VNTRs
^^^^^^^^^^
Under construction ...

.. _add-custom-vntr-label:

Add Custom VNTR
^^^^^^^^^^^^^^^
Use :code:`advntr addmodel [options]` to add a VNTR to the database.
The structure of VNTR and its genomic coordinate are required.

Required arguments:

+-----------------------+----------------------------------------------------------------+
| -r/--reference <text> | Reference genome                                               |
+-----------------------+----------------------------------------------------------------+
| -c/--chromosome <text>| Chromosome (e.g. chr1)                                         |
+-----------------------+----------------------------------------------------------------+
| -p/--pattern <text>   | First repeating pattern of VNTR in forward (5' to 3') direction|
+-----------------------+----------------------------------------------------------------+
| -s/--start <int>      | Start coordinate of VNTR in forward (5' to 3') direction       |
+-----------------------+----------------------------------------------------------------+
| -e/--end <int>        |  End coordinate of VNTR in forward (5' to 3') direction        |
+-----------------------+----------------------------------------------------------------+

Other options:

+-------------------------+--------------------------------+
| -g/--gene <text>        |Gene name                       |
+-------------------------+--------------------------------+
| -a/--annotation <text>  |Annotation of VNTR region       |
+-------------------------+--------------------------------+
| -h/--help               |show this help message and exit |
+-------------------------+--------------------------------+

You can use :code:`--update` in genotyping step to iteratively update the model using real data.


Delete a VNTR
^^^^^^^^^^^^^
Use :code:`advntr delmodel --vntr_id <ID>` to remove a VNTR from database.
