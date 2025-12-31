# Portfolio-safe copy (shell script)
# Identifiers generalized; secrets removed; paths parameterized.

#!/bin/bash

# Step 0: Create env folder with yamls
echo "Step 0: Exporting Conda environments to YAML files..."
[ ! -d env ] && mkdir env

for env in $(conda env list | awk '{print $1}' | tail -n +4); do
    echo "Exporting $env..."
    conda env export --name $env --no-builds > env/${env}.yaml
done
echo "YAML export completed."

# Step 1: Check Dependencies
echo "Step 1: Running check_dependencies.py..."
if python check_dependencies.py; then
    echo "check_dependencies.py completed successfully."
else
    echo "Error: check_dependencies.py encountered an issue." >&2
    exit 1
fi

# Step 2: Filter Dependencies with Bioconda Tools
echo "Step 2: Running bioconda_tools.py..."
if python bioconda_tools.py; then
    echo "bioconda_tools.py completed successfully."
else
    echo "Error: bioconda_tools.py encountered an issue." >&2
    exit 1
fi

# Step 3: Check for new versions and append updates to the YAMLs
echo "Step 3: Checking for new versions and appending updates..."
if python append_new.py; then
    echo "append_new.py completed successfully."
else
    echo "Error: append_new.py encountered an issue." >&2
    exit 1
fi

# Step 4: Recreate updated environments from the updated YAML files
echo "Step 4: Recreating environments from updated YAML files..."
if bash update_yamls.bash; then
    echo "Environment recreation completed successfully."
else
    echo "Error: update_yamls.bash encountered an issue." >&2
    exit 1
fi

echo "Automation workflow completed successfully!"
