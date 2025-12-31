#!/bin/bash
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH -o verify_%j.out
#SBATCH -e verify_%j.err

# --- THE FIX STARTS HERE ---
# 1. Force the shell to recognize Mamba functions
export MAMBA_ROOT_PREFIX="${HOME}/mambaforge"
eval "$(${MAMBA_ROOT_PREFIX}/bin/mamba shell hook -s bash)"

# 2. Activate the environment using the full path to be safe
mamba activate nextflow

# 3. Purge system Java and force the Mamba Java to the front
unset JAVA_HOME
unset _JAVA_OPTIONS
export PATH="${HOME}/mambaforge/envs/nextflow/bin:$PATH"
# --- THE FIX ENDS HERE ---

# 4. Setup Nextflow/Singularity paths
export NXF_SINGULARITY_CACHEDIR="${HOME}/singularity_cache"
export XDG_RUNTIME_DIR="${HOME}/xdr"
mkdir -p $NXF_SINGULARITY_CACHEDIR

# 5. Validation Check
echo "DEBUG: Which Java am I using?"
which java
java -version

echo "DEBUG: Which Nextflow am I using?"
which nextflow
nextflow -v

# 6. Run the Pull
nextflow pull nf-core/rnaseq -r 3.22.2
