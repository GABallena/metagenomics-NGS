# Portfolio-safe copy (shell script)
# Identifiers generalized; secrets removed; paths parameterized.

#!/bin/bash

# Array of packages and their corresponding environments
declare -A PACKAGES_ENVS=(
    [emboss]="emboss_env"
    [kma]="kma_env"
    [pandaseq]="pandaseq_env"
    [spades]="spades_env"
    [manta]="manta_env"
    [blast]="blast_env"
    [samtools]="samtools_env"
    [bowtie2]="bowtie2_env"
    [minimap2]="minimap2_env"
    [mash]="mash_env"
    [hmmer]="hmmer_env"
    [recycler]="recycler_env"
    [plasmidfinder]="plasmidfinder_env"
    [mob_suite]="mob_suite_env"
    [gapfiller]="gapfiller_env"
)

# Configure Bioconda channels globally
echo "Configuring Bioconda channels..."
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Loop through the packages and create individual environments
for package in "${!PACKAGES_ENVS[@]}"; do
    env_name="${PACKAGES_ENVS[$package]}"
    echo "Creating and setting up environment: $env_name for package: $package"

    # Create the environment and install the package
    conda create -n "$env_name" -y
    conda activate "$env_name"
    conda install -n "$env_name" -y "bioconda::$package"

    # Deactivate the environment
    conda deactivate
done

echo "All environments have been created and packages installed."
echo "To activate an environment, run: conda activate <environment_name>"
