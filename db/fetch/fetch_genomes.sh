#!/bin/bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
set -euo pipefail

INPUT="clean_pathogens.txt"
OUTDIR="reference_genomes"
mkdir -p "$OUTDIR"

# Pre-count total lines
total=$(grep -cve '^\s*$' "$INPUT")
count=0

while IFS= read -r pathogen; do
    ((count++))
    # Trim whitespace
    pathogen=$(echo "$pathogen" | sed 's/^[ \t]*//;s/[ \t]*$//')
    [ -z "$pathogen" ] && continue  # Skip empty lines

    echo -e "\n[$count/$total] Query: ${pathogen}" | lolcat 2>/dev/null || echo -e "\n[$count/$total] Query: ${pathogen}"

    safe_name=$(echo "$pathogen" | sed 's/ /_/g')
    target_dir="$OUTDIR/$safe_name"

    # Skip if already exists
    if [ -d "$target_dir" ]; then
        echo "✅ Already exists: $target_dir, skipping" | lolcat 2>/dev/null || echo "✅ Already exists: $target_dir, skipping"
        continue
    fi

    echo " > Attempting reference genome fetch..." | lolcat 2>/dev/null || echo " > Attempting reference genome fetch..."
    datasets download genome taxon "$pathogen" \
        --reference \
        --filename "${safe_name}.zip" 2>/dev/null

    # Fallback if reference doesn't exist
    if [ ! -f "${safe_name}.zip" ]; then
        echo " > Reference unavailable. Trying representative genome..." | lolcat 2>/dev/null || echo " > Reference unavailable. Trying representative genome..."
        datasets download genome taxon "$pathogen" \
            --representative \
            --filename "${safe_name}.zip" 2>/dev/null
    fi

    # Extraction
    if [ -f "${safe_name}.zip" ]; then
        echo " > Extracting..." | lolcat 2>/dev/null || echo " > Extracting..."
        unzip -q "${safe_name}.zip" -d "$target_dir"
        rm -f "${safe_name}.zip"
        echo " > ✅ Downloaded and extracted: $target_dir" | lolcat 2>/dev/null || echo " > ✅ Downloaded and extracted: $target_dir"
    else
        echo " ❌ No genome available for: $pathogen" | lolcat 2>/dev/null || echo " ❌ No genome available for: $pathogen"
    fi

    sleep 1
done < "$INPUT"

echo -e "\n🎉 All targets processed. Session complete." | lolcat 2>/dev/null || echo -e "\n🎉 All targets processed. Session complete."
