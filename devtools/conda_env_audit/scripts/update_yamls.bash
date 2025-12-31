# Portfolio-safe copy (shell script)
# Identifiers generalized; secrets removed; paths parameterized.

#!/bin/bash

yaml_dir="env/"

for yaml_file in "$yaml_dir"*.yaml; do
    env_name=$(basename "$yaml_file" .yaml)
    echo "Recreating environment: $env_name"
    conda env create -f "$yaml_file" --force
done
