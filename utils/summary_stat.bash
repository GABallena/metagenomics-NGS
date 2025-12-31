#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# Optional: activate FastQC environment if conda is available
if command -v conda >/dev/null 2>&1; then
  # In non-interactive shells, conda may require initialization (source conda.sh) in your environment.
  conda activate fastqc_env >/dev/null 2>&1 || true
fi

# Directories
RAW_READS_DIR="raw_reads"
TRIMMED_READS_DIR="trimmed_reads"
SUMMARY_DIR="summary_statistics"
FASTQC_RAW_DIR="$SUMMARY_DIR/fastqc_raw"
FASTQC_TRIMMED_DIR="$SUMMARY_DIR/fastqc_trimmed"
EXTRACTION_RESULTS_DIR="$SUMMARY_DIR/extracted_metrics"

# Create output directories
echo "Creating necessary directories..."
mkdir -p "$FASTQC_RAW_DIR" "$FASTQC_TRIMMED_DIR" "$EXTRACTION_RESULTS_DIR"
echo "Directories created."

count_reads() {
  local in_dir="$1"
  local out_file="$2"
  echo "Counting reads in ${in_dir}..."
  : > "$out_file"

  local files=( "${in_dir}"/*.fastq.gz )
  if [ ${#files[@]} -eq 0 ]; then
    echo "No .fastq.gz files found in ${in_dir}"
    return 0
  fi

  for file in "${files[@]}"; do
    echo "Processing file: ${file}"
    # FASTQ has 4 lines per read
    local line_count
    line_count=$(zcat "$file" | wc -l)
    local read_count=$(( line_count / 4 ))
    echo "${file}: ${read_count} reads" >> "$out_file"
  done

  echo "Read counting completed for ${in_dir}."
}

extract_fastqc_module() {
  # Args: module_name, fastqc_data.txt, output_file
  local module="$1"
  local data_file="$2"
  local out_file="$3"

  awk -v name="$module" '
    $0 ~ "^>>"name {flag=1; next}
    flag && $0 ~ "^>>END_MODULE" {exit}
    flag {print}
  ' "$data_file" > "$out_file"
}

extract_quality_metrics() {
  local fastqc_dir="$1"
  local out_dir="$2"
  echo "Extracting quality metrics from FastQC results in ${fastqc_dir}..."
  mkdir -p "$out_dir"

  local zips=( "${fastqc_dir}"/*_fastqc.zip )
  if [ ${#zips[@]} -eq 0 ]; then
    echo "No FastQC zip files found in ${fastqc_dir}"
    return 0
  fi

  for file in "${zips[@]}"; do
    echo "Unzipping ${file}..."
    unzip -q -o "$file" -d "$fastqc_dir/"

    local summary_dir="${file%.zip}"
    local summary_file="${summary_dir}/fastqc_data.txt"
    local base_name
    base_name="$(basename "${file%_fastqc.zip}")"

    local base_output_file="${out_dir}/${base_name}_per_base_quality.txt"
    local seq_output_file="${out_dir}/${base_name}_per_sequence_quality.txt"

    if [ ! -f "$summary_file" ]; then
      echo "WARNING: Missing ${summary_file}; skipping."
      continue
    fi

    extract_fastqc_module "Per base sequence quality" "$summary_file" "$base_output_file"
    extract_fastqc_module "Per sequence quality scores" "$summary_file" "$seq_output_file"
  done

  echo "Quality metrics extraction completed for ${fastqc_dir}."
}

# Read counts
count_reads "$RAW_READS_DIR" "$SUMMARY_DIR/raw_read_counts.txt"
count_reads "$TRIMMED_READS_DIR" "$SUMMARY_DIR/trimmed_read_counts.txt"

# FastQC on raw reads
echo "Running FastQC on raw reads..."
raw_files=( "$RAW_READS_DIR"/*.fastq.gz )
if [ ${#raw_files[@]} -gt 0 ]; then
  fastqc "${raw_files[@]}" -o "$FASTQC_RAW_DIR"
else
  echo "No raw reads found for FastQC."
fi
echo "FastQC completed for raw reads."

# FastQC on trimmed reads
echo "Running FastQC on trimmed reads..."
trim_files=( "$TRIMMED_READS_DIR"/*.fastq.gz )
if [ ${#trim_files[@]} -gt 0 ]; then
  fastqc "${trim_files[@]}" -o "$FASTQC_TRIMMED_DIR"
else
  echo "No trimmed reads found for FastQC."
fi
echo "FastQC completed for trimmed reads."

# Extract metrics
extract_quality_metrics "$FASTQC_RAW_DIR" "$EXTRACTION_RESULTS_DIR"
extract_quality_metrics "$FASTQC_TRIMMED_DIR" "$EXTRACTION_RESULTS_DIR"

echo "Read counting, FastQC analysis, and metric extraction complete."
echo "All output files have been saved in the ${SUMMARY_DIR} directory."
