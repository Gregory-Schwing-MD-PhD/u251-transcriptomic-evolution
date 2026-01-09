#!/bin/bash
#SBATCH -q secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=3-00:00:00
#SBATCH -o xengsort_prod_%j.out
#SBATCH -e xengsort_prod_%j.err

# 1. Environment Setup
source "${HOME}/mambaforge/etc/profile.d/conda.sh"
source activate nextflow

# 2. Step 1: Species Sorting (Custom main.nf)
nextflow run main.nf -profile singularity -resume

# 3. Step 2: RNA-seq Alignment & Quantification
nextflow run nf-core/rnaseq \
    -r 3.14.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/human_evolution_paired.csv" \
    --outdir "$(pwd)/ANALYSIS/results_human_final" \
    --fasta "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --bbsplit_fasta_list "$(pwd)/ANALYSIS/bbsplit.csv" \
    --skip_bbsplit false \
    --save_bbsplit_reads \
    --save_unaligned \
    --extra_star_align_args '--outFilterMultimapNmax 1' \
    --aligner star_salmon \
    --max_cpus 16 \
    --max_memory '64.GB' \
    -resume

# 4. Step 3: Differential Abundance & Clustering
nextflow run nf-core/differentialabundance \
    -r 1.5.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/metadata_therapy.csv" \
    --contrasts "$(pwd)/ANALYSIS/contrasts_therapy.csv" \
    --matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_counts.tsv" \
    --transcript_length_matrix "$(pwd)/ANALYSIS/results_human_final/star_salmon/salmon.merged.gene_lengths.tsv" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --exploratory_main_variable Classification \
    --shinyngs_build_app \
    -params-file therapy_params.yaml \
    -resume
