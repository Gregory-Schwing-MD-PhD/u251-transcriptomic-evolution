#!/bin/bash

# 1. Define the Folder URL
FOLDER_URL="https://drive.google.com/drive/u/0/folders/1bOFogbLTm_i-JQRfNKiLcOCNHOyBLUT7"

# 2. Setup Directory
mkdir -p downloaded_files
cd downloaded_files

# 3. Check for cookies.txt (Crucial for old gdown versions)
if [ -f "cookies.txt" ]; then
    echo "🍪 Found 'cookies.txt' in current directory."
    echo "   Your version of gdown should auto-detect this file."
else
    echo "⚠️  WARNING: 'cookies.txt' NOT found in $(pwd)"
    echo "   Please upload 'cookies.txt' to this folder to bypass the quota error."
fi

# 4. Download Loop
echo "Starting download..."
echo "---------------------------------------------"

while true; do
    # REMOVED "--cookies cookies.txt" because your version doesn't support it.
    # It will attempt to auto-load cookies.txt from the current folder.
    gdown "$FOLDER_URL" --folder --continue

    if [ $? -eq 0 ]; then
        echo "✅ Download successful."
        break
    else
        echo "❌ Download failed or interrupted."
        echo "   If you see 'Too many users', ensure cookies.txt is valid and in this folder."
        echo "⚠️  Retrying in 60 seconds..."
        sleep 60
    fi
done

echo "Done."
