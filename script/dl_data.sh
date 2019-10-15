#!/bin/bash

mkdir -p data

# Streptococcus pneumoniae PRJNA521678
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR855/006/SRR8556426/SRR8556426_1.fastq.gz | seqtk seq -A - > data/SRR8556426.fasta

# Candida vartiovaarae
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/006/ERR1877966/ERR1877966.fastq.gz | seqtk seq -A - > data/ERR1877966.fasta
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/007/ERR1877967/ERR1877967.fastq.gz | seqtk seq -A - > data/ERR1877967.fasta
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/008/ERR1877968/ERR1877968.fastq.gz | seqtk seq -A - > data/ERR1877968.fasta
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/009/ERR1877969/ERR1877969.fastq.gz | seqtk seq -A - > data/ERR1877969.fasta
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR187/000/ERR1877970/ERR1877970.fastq.gz | seqtk seq -A - > data/ERR1877970.fasta

# Escherichia coli CFT073
curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR849/000/SRR8494940/SRR8494940_1.fastq.gz | seqtk seq -A - > data/SRR8494940.fasta
