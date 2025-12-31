#!/bin/bash

# 1. Create the 'refs' subdirectory if it doesn't exist
mkdir -p refs

echo "Starting downloads..."

# 2. Download Human GRCh38
echo "Downloading Human GRCh38 to refs/human_GRCh38.fa.gz..."
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz \
  -O refs/human_GRCh38.fa.gz

# 3. Download Rat mRatBN7.2
echo "Downloading Rat mRatBN7.2 to refs/rat_mRatBN7.2.fa.gz..."
wget http://ftp.ensembl.org/pub/release-110/fasta/rattus_norvegicus/dna/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz \
  -O refs/rat_mRatBN7.2.fa.gz

echo "All downloads complete."
