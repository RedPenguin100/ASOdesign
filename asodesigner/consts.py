import os
from pathlib import Path

PROJECT_PATH = Path(os.path.dirname(os.path.dirname(__file__)))

TEST_PATH = PROJECT_PATH

DATA_PATH = PROJECT_PATH / 'data'

GFP1_PATH = DATA_PATH / 'gfp1_seq.txt'
GFP_FIRST_EXP_FASTA = DATA_PATH / 'GFP_first_exp.fasta'

# Yeast
YEAST_DATA = DATA_PATH / 'yeast' / 'yeast_data'
YEAST_FASTA_PATH = YEAST_DATA / 'GCF_000146045.2' / 'GCF_000146045.2_R64_genomic.fna'
YEAST_GFF_PATH = YEAST_DATA / 'GCF_000146045.2' / 'genomic.gff'
YEAST_GFF_DB_PATH = YEAST_DATA / 'dbs' / 'yeast_gff.db'
YEAST_FIVE_PRIME_UTR = YEAST_DATA / 'SGD_all_ORFs_5prime_UTRs.fsa'
YEAST_THREE_PRIME_UTR = YEAST_DATA / 'SGD_all_ORFs_3prime_UTRs.fsa'

# Human
HUMAN_DATA = DATA_PATH / 'human_dna'
HUMAN_GTF_BASIC = HUMAN_DATA / 'gencode.v34.basic.annotation.gtf'
HUMAN_GFF = HUMAN_DATA / 'gencode.v34.annotation.gff3'
HUMAN_GFF_GZ = HUMAN_DATA / 'gencode.v34.annotation.gff3.gz'
HUMAN_TRANSCRIPTS_FASTA_GZ = HUMAN_DATA / 'gencode.v34.transcripts.fa.gz'
HUMAN_TRANSCRIPTS_FASTA = HUMAN_DATA / 'gencode.v34.transcripts.fa'
HUMAN_GENOME_FASTA_GZ = HUMAN_DATA / 'GRCh38.p13.genome.fa.gz'
HUMAN_GENOME_FASTA = HUMAN_DATA / 'GRCh38.p13.genome.fa'
HUMAN_DB_PATH = HUMAN_DATA / 'dbs'
HUMAN_DB = HUMAN_DB_PATH / 'human_gff.db'
HUMAN_DB_BASIC = HUMAN_DB_PATH / 'human_gff_basic.db'
HUMAN_DB_BASIC_INTRONS = HUMAN_DB_PATH / 'human_gff_basic_introns.db'
