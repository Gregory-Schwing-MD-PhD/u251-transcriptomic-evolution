#!/bin/bash

# 1. Define the Folder URL
FOLDER_URL="https://drive.google.com/drive/u/0/folders/1bOFogbLTm_i-JQRfNKiLcOCNHOyBLUT7"

# 2. Prepare the directory
mkdir -p downloaded_files

# 3. Handle cookies.txt
# If cookies.txt exists in the current folder, copy it into the download folder
if [ -f "cookies.txt" ]; then
    echo "🍪 Found cookies.txt. Copying it to 'downloaded_files/'. "
    cp cookies.txt downloaded_files/
else
    echo "⚠️  WARNING: cookies.txt not found in current directory."
    echo "   The download might fail with 'Too many users' errors."
fi

# 4. Move into the directory
cd downloaded_files

# 5. Download Loop
echo "Starting download..."
echo "---------------------------------------------"

while true; do
    # gdown (older versions) will now find cookies.txt in this folder automatically
    gdown "$FOLDER_URL" --folder --continue

    if [ $? -eq 0 ]; then
        echo "✅ Download successful."
        break
    else
        echo "❌ Download failed or interrupted."
        echo "⚠️  Retrying in 60 seconds..."
        sleep 60
    fi
done

echo "Done."
