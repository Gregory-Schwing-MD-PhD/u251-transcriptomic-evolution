#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH -o xengsort_prod_%j.out
#SBATCH -e xengsort_prod_%j.err

# Trace every command to the .err file
set -x

# 1. Environment Setup
echo "LOG: Sourcing conda.sh"
source "${HOME}/mambaforge/etc/profile.d/conda.sh"

echo "LOG: Bypassing activate and setting PATH directly for Nextflow and Java"
# Manually pointing to the environment's bin folder to avoid Conda deadlock
export CONDA_PREFIX="${HOME}/mambaforge/envs/nextflow"
export PATH="${CONDA_PREFIX}/bin:$PATH"
export CONDA_DEFAULT_ENV="nextflow"

# Crucial: Unset JAVA_HOME so it defaults to the java in the PATH we just set
unset JAVA_HOME

# --- Sanity Checks ---
echo "LOG: Verifying Java Path..."
which java
java -version

echo "LOG: Verifying Nextflow Path..."
which nextflow
nextflow -v

# 2. Path Setup
echo "LOG: Setting XDG_RUNTIME_DIR"
export XDG_RUNTIME_DIR="${HOME}/xdr"
mkdir -p $XDG_RUNTIME_DIR

# 3. Singularity Image Prep
echo "LOG: Setting NXF_SINGULARITY_CACHEDIR"
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
mkdir -p $NXF_SINGULARITY_CACHEDIR

echo "Pre-pulling Singularity image from Docker Hub..."
singularity pull --name ${NXF_SINGULARITY_CACHEDIR}/go2432-xengsort-latest.img docker://go2432/xengsort:latest

# --- End of Setup / Turn off Trace ---
set +x

# 4. Clean up stale Nextflow locks
echo "LOG: Cleaning up locks"
find .nextflow/cache -name "LOCK" -delete 2>/dev/null

# 5. Production Run
echo "RUNNING STEP 1: XENGSORT"
# Passing parameters directly here since they were removed from nextflow.config
nextflow run main.nf -profile singularity \
    --input "ANALYSIS/samplesheet.csv" \
    --host_fasta "/wsu/home/go/go24/go2432/u251-transcriptomic-evolution/ANALYSIS/refs/rat/Rattus_norvegicus.mRatBN7.2.dna.toplevel.fa.gz" \
    --graft_fasta "/wsu/home/go/go24/go2432/u251-transcriptomic-evolution/ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --outdir "ANALYSIS/xengsort_out" \
    -resume -ansi-log false

# Clear parameters between pipelines to avoid "invalid input" warnings
unset NXF_PARAMS

echo "RUNNING STEP 2: RNA-SEQ"
nextflow run nf-core/rnaseq \
    -r 3.14.0 \
    -profile singularity \
    --input "$(pwd)/ANALYSIS/human_evolution_paired.csv" \
    --outdir "$(pwd)/ANALYSIS/results_human_final" \
    --fasta "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.genome.fa.gz" \
    --gtf "$(pwd)/ANALYSIS/refs/human/GRCh38.primary_assembly.annotation.gtf.gz" \
    --skip_bbsplit true \
    --save_unaligned \
    --extra_star_align_args '--outFilterMultimapNmax 1' \
    --aligner star_salmon \
    --max_cpus 16 \
    --max_memory '128.GB' \
    -resume \
    -ansi-log false

unset NXF_PARAMS

echo "RUNNING STEP 3: DIFFERENTIAL ABUNDANCE"
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
    -resume \
    -ansi-log false
