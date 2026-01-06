#!/bin/bash
# EXECUTE THIS SCRIPT INSIDE THE 'ANALYSIS' DIRECTORY
# Usage: bash create_rat_samplesheet.sh

OUT_FILE="samplesheet_host_filtered.csv"
BBSPLIT_DIR="results/bbsplit"

# 1. Safety Check
if [ ! -d "$BBSPLIT_DIR" ]; then
    echo "Error: Directory '$BBSPLIT_DIR' not found."
    echo "Make sure you are inside the ANALYSIS folder and Phase 1 has finished the 'bbsplit' stage."
    exit 1
fi

echo "Found bbsplit directory. Generating samplesheet..."
echo "sample,fastq_1,fastq_2,strandedness" > "$OUT_FILE"

# 2. Loop through the Rat-specific reads
# Note: Matches files ending in 'rat_ref_1.fastq.gz' (based on your bbsplit.csv key)
COUNT=0
for R1 in "$BBSPLIT_DIR"/*rat_ref_1.fastq.gz; do
    
    # Check if file exists (handles empty directories)
    [ -e "$R1" ] || continue

    # Get Full Absolute Path (Critical for Nextflow)
    ABS_R1="$(pwd)/$R1"
    
    # Determine R2 Path (Replace _1 with _2)
    ABS_R2="${ABS_R1/_1.fastq.gz/_2.fastq.gz}"

    # Extract Sample Name
    # Handles both "Sample.rat_ref..." and "Sample_rat_ref..." patterns
    FILENAME=$(basename "$R1")
    SAMPLE=$(echo "$FILENAME" | sed -E 's/[._]rat_ref_1.fastq.gz//')

    # 3. Filter out the Culture sample (C2B) - No rat cells there
    if [[ "$SAMPLE" == *"C2B"* ]]; then
        echo "  [SKIP] Excluding Culture sample: $SAMPLE"
        continue
    fi

    # 4. Append to CSV
    echo "${SAMPLE},${ABS_R1},${ABS_R2},auto" >> "$OUT_FILE"
    ((COUNT++))
done

echo "------------------------------------------------"
echo "Success! Added $COUNT samples to $OUT_FILE"
echo "------------------------------------------------"
