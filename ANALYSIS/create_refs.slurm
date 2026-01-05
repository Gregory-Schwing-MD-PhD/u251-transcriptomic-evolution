#!/bin/bash
#SBATCH --job-name=create_refs      # Job name
#SBATCH --output=create_refs_%j.log # Standard output and error log
#SBATCH -q primary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00               # Time limit hrs:min:sec

# 1. Ensure the script is executable
chmod +x create_refs.sh

module load r/4.5.0

# 2. Run the script
./create_refs.sh
