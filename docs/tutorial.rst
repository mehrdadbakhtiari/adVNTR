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
Under construction. In the meantime, please refer to the examples in the Installation page.


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


Delete a VNTR
^^^^^^^^^^^^^
Use :code:`advntr delmodel --vntr_id <ID>` to remove a VNTR from database.
