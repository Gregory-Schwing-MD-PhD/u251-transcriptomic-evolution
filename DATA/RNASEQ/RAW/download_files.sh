#!/bin/bash

# 1. Check if gdown is installed
if ! command -v gdown &> /dev/null; then
    echo "Error: 'gdown' is not installed or not in your PATH."
    echo "Please install it first by running: pip install gdown"
    exit 1
fi

# 2. Define the Folder URL
FOLDER_URL="https://drive.google.com/drive/u/0/folders/1bOFogbLTm_i-JQRfNKiLcOCNHOyBLUT7"

# 3. Create a directory for the downloads
# Note: gdown --folder will usually create a subfolder with the Drive folder's name inside this directory.
mkdir -p downloaded_files
cd downloaded_files

# 4. Download the folder with retry logic
echo "Starting download of the entire folder..."
echo "---------------------------------------------"

while true; do
    # --folder: tells gdown this is a folder link
    # --continue: resumes partially downloaded files if they fail
    # --remaining: (implied in newer gdown versions) skips files that are already fully downloaded
    gdown "$FOLDER_URL" --folder --continue

    # Check if the download was successful
    if [ $? -eq 0 ]; then
        echo "✅ Folder download completed successfully."
        echo "---------------------------------------------"
        break # Exit the loop
    else
        echo "❌ Download interrupted or failed."
        echo "⚠️  Network issue or Timeout detected. Sleeping for 60 seconds before retrying..."
        echo "   (It will verify existing files and resume where it left off)"
        sleep 60
    fi
done

echo "All tasks finished."
