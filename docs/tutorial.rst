Tutorial
========

Inputs
------
- NGS short reads (Illumina HiSeq)
- SMRT reads (PacBio)

Outputs
-------
Currently there are two possible formats to get the genotyping output:

- Text
    Writes two lines in the output for each VNTR. The first contains the VNTR ID and the second line contains R1/R2
    as the repeating unit counts.
    Below is an example output in text format for one VNTR:

  | 301645
  | 2/3
  |

- BED
    BED format contains one line per locus and it is a tab-delimited output comprised of 9 columns: 1. The name of the
    chromosome, 2. Start position of the VNTR, 3. End position of the VNTR, 4. VNTR ID, 5. Name of the gene that
    contains the VNTR, 6. Repeating motif, 7. Number of repeats in reference genome, 8 and 9. Number of repeats in the
    sample.
    Below is an example output in BED format for one VNTR:

+--------+---------+---------+---------+------+-------------+---------+---+---+
| #CHROM | Start   | End     | VNTR_ID | Gene | Motif       | RefCopy | R1| R2|
+--------+---------+---------+---------+------+-------------+---------+---+---+
| chr21  |45196324 |45196360 | 301645  | CSTB |CGCGGGGCGGGG | 3       |2  | 3 |
+--------+---------+---------+---------+------+-------------+---------+---+---+


- VCF
    (Under construction)

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

+-------------------------------+---------------------------------------------------------------------------------+
| -a/--alignment_file <file>    | alignment file in SAM/BAM/CRAM format                                           |
+-------------------------------+---------------------------------------------------------------------------------+
| -r/--reference_filename <file>| path to a FASTA-formatted reference file for CRAM files.                        |
+-------------------------------+---------------------------------------------------------------------------------+
| -f/--fasta <file>             | Fasta file containing raw reads                                                 |
+-------------------------------+---------------------------------------------------------------------------------+
| -p/--pacbio                   | set this flag if input file contains PacBio reads instead of Illumina reads     |
+-------------------------------+---------------------------------------------------------------------------------+
| -n/--nanopore                 | set this flag if input file contains Nanopore MinION reads instead of Illumina  |
+-------------------------------+---------------------------------------------------------------------------------+
| -o/--outfile <file>           | file to write results. adVNTR writes output to stdout if oufile is not specified|
+-------------------------------+---------------------------------------------------------------------------------+
| -of/--outfmt <format>         | output format. Allowed values are {text, bed} [text]                            |
+-------------------------------+---------------------------------------------------------------------------------+

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
