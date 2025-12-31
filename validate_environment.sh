#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH -o verify_%j.out
#SBATCH -e verify_%j.err

# 1. Load the environment
source "${HOME}/mambaforge/etc/profile.d/mamba.sh"
mamba activate nextflow

# 2. Setup Singularity Cache (important to test if this path is writable)
export NXF_SINGULARITY_CACHEDIR=${HOME}/singularity_cache
mkdir -p $NXF_SINGULARITY_CACHEDIR

# 3. Verify tool versions
echo "--- Checking Software Versions ---"
nextflow -v
singularity --version

# 4. Pull the pipeline (tests GitHub API connection)
echo "--- Pulling Pipeline ---"
nextflow pull nf-core/rnaseq -r 3.22.2

echo "--- Check Complete ---"
