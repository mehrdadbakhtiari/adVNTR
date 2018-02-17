.. _faq:

FAQ
===


.. contents::
  :local:


How do I cite adVNTR?
---------------------
    If you found adVNTR useful, we would appreciate it if you could cite our manuscript describing adVNTR:

    .. include:: publication.rst

Can adVNTR work with repeating units that are shorter than 6bp?
---------------------------------------------------------------
    Tandem repeats with period below 6bp are classified as Short Tandem Repeats (STRs). Although adVNTR can detect STRs
    expansions, we do not recommend to use it on STRs.

Can I run adVNTR to study expansion in other organisms?
-------------------------------------------------------
    You can run adVNTR for other organisms if you add custom VNTR to its database. However, it always returns diploid
    RU counts for the number of repeats and it is expected to get homozygous RU counts on haploid organisms.

What sequencing platforms does adVNTR support?
----------------------------------------------
    adVNTR is designed to analyze **Illumina** or **PacBio** sequencing data. We generally do not recommend to use it on
    sequencing data from other technologies as their error model is different.
