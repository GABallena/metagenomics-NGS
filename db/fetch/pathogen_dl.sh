#!/bin/bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
set -euo pipefail

mkdir -p reference_genomes
clear

echo -e "\nInitializing environment..." | lolcat
sleep 1
echo -e "Syncing remote indices..." | lolcat
sleep 1
echo -e "Parsing organism identifiers..." | lolcat
sleep 1

total=$(($(wc -l < gcpathogen_complete.tsv)-1))
count=0

tail -n +2 gcpathogen_complete.tsv | while IFS=$'\t' read -r pathogen _; do
    ((count++))
    echo -e "\n[$count/$total] Query: ${pathogen}" | lolcat
    safe_name=$(echo "$pathogen" | sed 's/ /_/g')

    echo " > Attempting reference genome fetch..." | lolcat
    datasets download genome taxon "$pathogen" \
        --reference --dehydrated \
        --filename "${safe_name}.zip" 2>/dev/null

    if [ ! -f "${safe_name}.zip" ]; then
        echo " > Reference unavailable. Attempting fallback strategy..." | lolcat
        datasets download genome taxon "$pathogen" \
            --representative --dehydrated \
            --filename "${safe_name}.zip" 2>/dev/null
    fi

    if [ -f "${safe_name}.zip" ]; then
        echo " > Extraction in progress..." | lolcat
        unzip -q "${safe_name}.zip" -d "reference_genomes/${safe_name}"
        rm "${safe_name}.zip"
        echo " > Archive processed: ${safe_name}" | lolcat
    else
        echo " > No available assembly: ${safe_name}" | lolcat
    fi

    sleep 1
done

echo -e "\nAll targets processed. Session complete." | lolcat
