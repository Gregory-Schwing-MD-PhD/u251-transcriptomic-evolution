#!/bin/bash

# 1. Check if gdown is installed
if ! command -v gdown &> /dev/null; then
    echo "Error: 'gdown' is not installed."
    echo "Please run: pip install gdown"
    exit 1
fi

# 2. Define the Folder URL
FOLDER_URL="https://drive.google.com/drive/u/0/folders/1bOFogbLTm_i-JQRfNKiLcOCNHOyBLUT7"

# 3. Check for cookies.txt
# This is crucial for bypassing the "Too many users" error
if [ -f "cookies.txt" ]; then
    echo "🍪 Found cookies.txt! Using authenticated download."
    COOKIE_FLAG="--cookies cookies.txt"
else
    echo "⚠️  cookies.txt NOT found."
    echo "   You are attempting an anonymous download."
    echo "   If you get 'Too many users' errors, please upload cookies.txt to this folder."
    COOKIE_FLAG=""
fi

# 4. Setup Directory
mkdir -p downloaded_files
cd downloaded_files

# 5. Download Loop
echo "Starting download..."
echo "---------------------------------------------"

while true; do
    # We pass the $COOKIE_FLAG variable into the command
    # --folder: Download whole folder
    # --continue: Resume partial downloads
    gdown "$FOLDER_URL" --folder --continue $COOKIE_FLAG

    if [ $? -eq 0 ]; then
        echo "✅ Download successful."
        break
    else
        echo "❌ Download failed or interrupted."
        echo "   (If the error is 'Permission denied' with cookies, your session may have expired."
        echo "    Export a fresh cookies.txt and try again.)"
        echo "⚠️  Retrying in 60 seconds..."
        sleep 60
    fi
done

echo "Done."
