# ToxinFinderPy
## ToxinFinderPy.py
THis script discribes a highthrough method for toxin gene searching in BT's genome or Metagenome.
## Dependency
### 1.Python2
### 2.Profigal
### 3.Blast+

usage: ToxinFinderPy.py [-h] -i INPUT_FILES [INPUT_FILES ...] -db LOCAL_DB -o
                        OUTPUT_FILES [-t THREADS_NUM]
a script to find the toxin gene in BT's genome or metagenome
optional arguments:

  -h, --help                          show this help message and exit
  
  -i INPUT_FILES                      input genome file in fasta(required)
  
  -db LOCAL_DB                        local database_files for blast in fasta(required)
  
  -o OUTPUT_FILES                     the output file (required)
  
  -t THREADS_NUM                      the num of threads used for blast

## NCBIBLAST_and_parse.py
This script enbale to run NCBI Blast automatically

usage: NCBIBLAST_and_parse.py [-h] [-i INPUT_FASTA [INPUT_FASTA ...]] -o OUTPUT_DIR

a script to running NCBI blast automated

optional arguments:

  -h, --help                          show this help message and exit
  
  -i INPUT_FASTA                      input AA fasta file
  
  -o OUTPUT_DIR                       the output dir (required)
