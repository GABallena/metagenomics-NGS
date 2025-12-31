#!/bin/bash
#
# Portfolio-safe copy: paths/identifiers generalized; databases not included.
#
set -euo pipefail

# 1. Extract pathogenwatch names (clean spaces)
cut -d'_' -f1,2 --output-delimiter=' ' pathogenwatch.list | sed 's/__.*//' | sed 's/[[:space:]]*$//' | sort -u > tmp_watch.txt

# 2. Extract gcpathogen names (skip header)
tail -n +2 gcpathogen.tsv | cut -f1 | sed 's/[[:space:]]*$//' | sort -u > tmp_gc.txt

# 3. Combine and deduplicate
cat tmp_watch.txt tmp_gc.txt | sed 's/[[:space:]]\+/ /g' | sort -u > combined_pathogens.txt

# 4. Cleanup
rm tmp_watch.txt tmp_gc.txt

echo "✅ Cleaned and deduplicated list saved to: combined_pathogens.txt"
