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

You can install python required packages by running ```pip install -r requirements.txt```

Execution:
===========
Use following command to see the help for running the tool.
```sh
python main.py --help
```

Demo 1: input in [BAM] format (aligned reads), genotyping RU counts
===========
```--alignment_file``` specifies the alignment file containing mapped and unmapped reads:
```sh
python main.py --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir
```
With ```--pacbio```, adVNTR assumes the alignment file contains PacBio sequencing data:
```sh
python main.py --alignment_file aligned_pacbio_reads.bam --working_directory ./log_dir/ --pacbio
```sh

Demo 2: input in [BAM] format (aligned reads), finding frameshift in VNTR
===========
With ```--frameshift```, adVNTR searches for frameshift in VNTR instead of finding RU count:
```sh
python main.py --alignment_file aligned_illumina_reads.bam --working_directory ./log_dir/ --frameshift
```sh
