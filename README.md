Software Requirements
==========
1. Following libraries are reguired
    -   ```python-pip```
    -   ```python-tk```
    -   ```libz-dev```

You can install these requirement in linux by running ```sudo apt-get install python-pip python-tk libz-dev```

2. ```Python2.7``` and the following packages are required:
    -   ```biopython```
    -   ```pysam``` version 0.9.1.4 or above
    -   ```cython```
    -   ```networkx``` version 1.11
    -   ```scipy```
    -   ```joblib```

You can install required python libraries by running ```pip install -r requirements.txt```

3. ```ncbi-blast``` version 2.2.29 or above is required


Data Requirements
===========
* To run adVNTR on trained VNTR models:
    - Download [vntr_data.zip](https://cseweb.ucsd.edu/~mbakhtia/adVNTR/vntr_data.zip) and extract it inside the project directory.

Execution:
===========
Use following command to see the help for running the tool.
```sh
python advntr.py --help
```
The program outputs the RU count genotypes for all VNTRs in ```vntr_data``` directory. To specify a single VNTR by its ID use ```--vntr_id <id>``` option. 

Demo 1: input in [BAM] format
===========
* ```--alignment_file``` specifies the alignment file containing mapped and unmapped reads:
```sh
python advntr.py --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/
```
* With ```--pacbio```, adVNTR assumes the alignment file contains PacBio sequencing data:
```sh
python advntr.py --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio
```
* Use ```--frameshift``` to find the possible frameshifts in VNTR:
```sh
python advntr.py --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift
```

Demo 2: input in [fasta] format
===========
* Use the following command to genotype the RU count using fasta file:
```sh
python advntr.py --fasta unaligned_illumina_reads.fasta --working_directory ./log_dir/
```
