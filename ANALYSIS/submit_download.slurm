#!/bin/bash
#SBATCH --job-name=download_refs      # Job name
#SBATCH --output=download_refs_%j.log # Standard output and error log
#SBATCH --nodes=1                     # Run all processes on a single node
#SBATCH --ntasks=1                    # Run a single task
#SBATCH --cpus-per-task=1             # Number of CPU cores per task
#SBATCH --mem=1G                      # Memory limit (1GB is plenty for downloads)
#SBATCH --time=01:00:00               # Time limit hrs:min:sec

# 1. Ensure the script is executable
chmod +x download_phase2_refs.sh

# 2. Run the script
./download_phase2_refs.sh
