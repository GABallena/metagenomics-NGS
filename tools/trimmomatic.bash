#!/bin/bash

# Activate the Conda environment for Trimmomatic
conda activate trimmomatic_env

# Trimmomatic settings
TRIMMOMATIC_ADAPTERS="NexteraPE-PE.fa"
TRIMMOMATIC_SETTINGS="LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:60"

# Input and output files passed as arguments
R1_INPUT=$1
R2_INPUT=$2
R1_PAIRED_OUTPUT=$3
R1_UNPAIRED_OUTPUT=${R1_PAIRED_OUTPUT/_paired/_unpaired}
R2_PAIRED_OUTPUT=$4
R2_UNPAIRED_OUTPUT=${R2_PAIRED_OUTPUT/_paired/_unpaired}

# Run Trimmomatic
trimmomatic PE -phred33 \
  $R1_INPUT $R2_INPUT \
  $R1_PAIRED_OUTPUT $R1_UNPAIRED_OUTPUT \
  $R2_PAIRED_OUTPUT $R2_UNPAIRED_OUTPUT \
  ILLUMINACLIP:$TRIMMOMATIC_ADAPTERS:2:30:10:1:true $TRIMMOMATIC_SETTINGS

echo "Processing complete for $R1_INPUT and $R2_INPUT."
