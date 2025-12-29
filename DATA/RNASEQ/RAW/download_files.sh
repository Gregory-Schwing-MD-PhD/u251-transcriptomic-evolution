#!/bin/bash

# 1. Check if gdown is installed
if ! command -v gdown &> /dev/null; then
    echo "Error: 'gdown' is not installed or not in your PATH."
    echo "Please install it first by running: pip install gdown"
    exit 1
fi

# 2. Define the list of File IDs
FILE_IDS=(
    "1hHWpUuy7kEL6kPJwB8S_sYWko4mrxWi9"
    "1tMMAsvmjpGkgqZUUJjkX1mrTFSbcMYQa"
    "1wYbPQ0GrXSvPAwqwXAKy7UtdcUY0vMy9"
    "1yD5KS48ET4wLBfwIJqL4Aq2cFg8uC3lh"
    "1xkJzJWVWsxwbrNLIBVt4DwqWbAMZ6-Y6"
    "1ObR6hmAoZjsPIAZ3fQQuzQPoCQpXtjXa"
    "17Y7meKjwzZrs5gZwYoIT0mj_kN8ffXHp"
    "1REOAwjMkyN4Zvmcwa1LDQGojlbWNsLR5"
    "12L_fh_0O3PIj3WxfSbqw7IetAzh8JyM_"
    "1i6d4fmJrf0MSRVsTZ5EyoF_S90dx0tUM"
    "1iJ2nT7j_J5PUHbZcK3F0TaDaHdYnCEm0"
    "1YthLCcXW_QfeZn6AAxz3d8dgsxu4D1ZM"
    "1YIJVV2ztNec5qS00tabBAaYc99ZB1uXF"
    "14oiZ8Ko_eoMMf1zzq0cjsQiaVb5RExg2"
    "15Crkgjkypn-ASqpPjk0pofDTknrJ2GEG"
    "1TfWlimRyYU4ZHyvTGy1FjU-F9POiE_-c"
    "14MP9k0KXoCgjrAR7VhxTi4Yf6shiuOok"
    "12Xbm-2gcY1sNMV3RHZI3CdZ03Rss1g6I"
    "1BvPGioQA-LK8hOPveSjo6339RvPgF2Yw"
    "1SNbeJjJ1cIMBCkmSx6MRq9PV4JXp3PuT"
    "1_NZ4uNFYdi3mbxMOzwU94ZF5Y_7PrJCz"
)

# 3. Create a directory for the downloads (Optional, keeps things clean)
mkdir -p downloaded_files
cd downloaded_files

# 4. Loop through the array and download each file
echo "Starting download of ${#FILE_IDS[@]} files..."
echo "---------------------------------------------"

for id in "${FILE_IDS[@]}"; do
    echo "Downloading ID: $id"
    gdown "https://drive.google.com/uc?id=${id}"
    
    if [ $? -eq 0 ]; then
        echo "✅ Success"
    else
        echo "❌ Failed to download ID: $id"
    fi
    echo "---------------------------------------------"
done

echo "All tasks finished."
