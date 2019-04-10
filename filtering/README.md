adVNTR-Filtering
----------------
This tool is part of [adVNTR](https://github.com/mehrdadbakhtiari/adVNTR) tool to speedup its read recruitment process.
It searches a set of keywords from multiple repetitive loci in DNA in the pool of sequencing reads and select the reads
that could be mapped to each locus.

Installation
------------
This package is part of adVNTR tool. Please refer to [adVNTR installation](https://github.com/mehrdadbakhtiari/adVNTR#installation) page to install the tool.

Usage
-----
adVNTR-Filtering could be invoked from command line with ``adVNTR-Filtering``

It takes the file containing fasta reads as the input, reads the keywords from stdin, and writes the output to stdout.
An example of command to run adVNTR-Filtering is as follows:
```
./adVNTR-Filtering sequences.fa --min_matches 4 <keywords.txt > output.txt
```

The `keywords.txt` is a text file that contains one line for each locus. In each line, there are multiple entries that
are separated by whitespace. The first item corresponds to the ID of the locus and all following items are keywords
assigned to this locus. An example of `keywords.txt` would look like this:
```
1 ACCC CACC CCAC CCCA
2 TGGT TTGG GTTG GGTT
...
n TGGG GTGG GGTG GGGT
```
In this example, we are looking to recruit reads corresponding to (ACCC)<sub>N</sub> and (TGGT)<sub>N</sub> loci.
For the command above, a reads must contain at least 4 exact matches to any of the keywords of one locus in order to be recruited for that locus.

For this hypothetical sequences.fa, the output for each command is shown below. 
```
>one
ACCCNNNNNNNNNNNN
>two
ACCCACCCNNNNNNNN
>three
ACCCACCCNNNNCCCT
>four
ACCCACCCACCCACCC
>one_ACCC_one_TGGG
ACCCTGGGNNNNNNNN
```

```
[user@linux:~]$ ./adVNTR-Filtering sequences.fa --min_matches 4 <keywords.txt 
1 1 four_ACCC
2 0
four_ACCC ACCCACCCACCCACCC

[user@linux:~]$ ./adVNTR-Filtering sequences.fa --min_matches 1 <keywords.txt 
1 5 four_ACCC two_ACCC three_ACCC one_ACCC_one_TTGG one_ACCC
2 1 one_ACCC_one_TTGG
four_ACCC ACCCACCCACCCACCC
one_ACCC ACCCNNNNNNNNNNNN
one_ACCC_one_TTGG ACCCTTGGNNNNNNNN
three_ACCC ACCCACCCNNNNCCCT
two_ACCC ACCCACCCNNNNNNNN
```

For each locus, there is one line starting with its ID followed by an integer `T`, number of reads that match to it.
The next `T` items in the line are ID of reads that match to the locus.

At the end, pair of fasta headers and fasta sequences of all reads that match to at least one locus will be written. (Each fasta record will appear at most once)
