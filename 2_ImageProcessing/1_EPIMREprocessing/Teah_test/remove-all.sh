#!/bin/bash

#rm -f ~/Desktop/* 
#rm -rf ~/Desktop/* 
echo "The following files were deleted on $(date):">> ~/Desktop/removed-files.txt
ls ~/Desktop -1  >> ~/Desktop/removed-files.txt
ls ~/Download -1  >> ~/Desktop/removed-files.txt

echo ""

find ~/Desktop -type f ! -name 'removed-files.txt' ! -name 'remove-all.sh' -exec rm -f {} +
find ~/Desktop -type d ! -name 'McIlvainDrive2' -exec rm -rf {} +

find ~/Downloads -type f ! -name 'removed-files.txt' ! -name 'remove-all.sh' -exec rm -f {} +
find ~/Downloads -type d ! -name 'McIlvainDrive2' -exec rm -rf {} +
