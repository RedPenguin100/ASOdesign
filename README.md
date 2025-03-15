# ASO Designer

This is the 2025 TAU iGEM project, dedicated to fighting cancer with ASO - antisense oligonucleotides.

The project is designed to run primarily on Linux systems.

# Setup

## Python
Steps:
1. Install Linux, preferably Ubuntu 24.04 
2. Install python3.11
3. Install miniconda
4. Create miniconda environment
5. Install dependencies using the miniconda environment `conda env create -f environment.yml`
6. "Install" the `asodesigner` package as `conda develop .`
7. Run `pytest` from the project root. If it works, then the setup is done!


## Data
The data is not saved directly in git because it is too large.

### Human
The human transcriptome is available in this link  https://www.gencodegenes.org/human/release_34.html <br/>
and the exact download link used was this:
https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz

You should download and unzip inside the data/human_dna folder.

### Yeast
For S. Cerevisiae I did not find a transcriptome as the only available data was ORFs that did not contain the 3'/5' UTR.

The S. Cerevisiae genome is the RefSeq Assembly downloaded via ncbi:
https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/ <br>
Select the Genome assembly, RefSeq only
Unzip the result to create a GCF folder inside data/yeast/yeast_data folder


For the 3/5 UTR information I used 
http://sgd-archive.yeastgenome.org/sequence/S288C_reference/ <br>

Download the files named "SGD_all_ORFs_3prime_UTRs.fsa.zip", "SGD_all_ORFs_5_prime_UTRs.fsa"

There is no account for alternative splicing in yeast, but it is extremely rare.

## External dependencies

### RISearch

When cloning, either
1. Clone recursively with `git clone --recursive`
2. After cloning, run git submodule update --init --recursive 

This will initialize the RIsearch repository that is linked this under this repository

To compile you need to install first gcc, make, cmake.

To compile RISearch1, simply enter the folder external/risearch/RIsearch1 and run `make`. The executable will be created.

To compile RIsearch2 you need to install first `sudo apt-get install libpcre3-dev` and then run `./rebuild.sh`
