# Portfolio-safe copy (python script)
# Identifiers generalized; secrets removed; paths parameterized.

import pandas as pd
from datetime import datetime
import sys
from pathlib import Path

# Add keyword lists as constants
INCLUSION_KEYWORDS = [
    "phylo", "kmer", "populat", "metagen", "antimicrob", "antibio", "resistance",
    "ortholog", "paralog", "bayesian", "trim", "read", "mapp", "card", "amr", "illumina",
    "assembl", "cluster", "mag", "genom", "pipeline", "alignment", "graph", "tensor",
    "informat", "algorithm", "geograph", "phage", "bootstrap", "statist", "pathogen",
    "quantif", "gpu", "loop", "gatk", "flanking", "slurm"
]

EXCLUSION_KEYWORDS = [
    "rnaseq", "nanopore", "pacbio", "longread", "cancer", "16s", "its", "singlecell",
    "pcr", "minion", "server", "windows", "plant", "mouse", "neural", "epigen", "soma",
    "chrom", "crispr", "expression", "cellculture", "library", "api", "tissue", "minion",
    "proteome", "lcms", "massspec", "microarray", "rnai", "rnaseq", "chipseq",
    "wholegenome", "msms", "hic", "hiseq", "hichip", "cytometry", "methylat",  
    "10x", "10xgenomics", "hg19", "hg38", "mm10", "mm9", "grch37", "grch38",
    "singlemolecule", "singlecell", "radseq", "rad", "wgs", "wes", "wholeexome"
]

def load_package_data(file_path):
    """
    Load and validate package data from TSV file.
    
    Args:
        file_path (str): Path to TSV file
        
    Returns:
        pandas.DataFrame: Loaded and processed dataframe
    """
    try:
        if not Path(file_path).exists():
            raise FileNotFoundError(f"Input file not found: {file_path}")
            
        # Read TSV with first row as header
        df = pd.read_csv(file_path, sep='\t')
        
        # Validate required columns
        required_cols = ['Package_Name', 'Description', 'Updated_Date']
        if not all(col in df.columns for col in required_cols):
            raise ValueError("Missing required columns in input file")
        
        # Rename columns to match expected names
        df = df.rename(columns={'Package_Name': 'Package'})
        
        # Apply package filtering with keywords
        df['Package'] = df['Package'].apply(filter_package_name)
        df = df.dropna(subset=['Package'])  # Remove filtered out packages
            
        # Convert date and validate 
        df['Updated_Date'] = pd.to_datetime(df['Updated_Date'], format='%Y-%m-%d')
        
        return df
        
    except Exception as e:
        print(f"Error loading package data: {str(e)}")
        sys.exit(1)

def filter_package_name(name):
    """
    Filter and standardize package names using inclusion/exclusion keywords.
    Returns None if package should be filtered out.
    """
    if not name or len(name.strip()) < 2:
        return None
        
    # Normalize name
    name = str(name).lower()
    name = name.replace('.', '-').replace('_', '-').replace(' ', '-')
    
    # Clean up the name first
    clean_name = ''.join(c for c in name if c.isalnum() or c == '-')
    clean_name = clean_name.strip('-')
    while '--' in clean_name:
        clean_name = clean_name.replace('--', '-')
        
    # Check inclusion keywords (at least one must match)
    has_inclusion = any(keyword in clean_name for keyword in INCLUSION_KEYWORDS)
    if not has_inclusion:
        return None
        
    # Check exclusion keywords (none should match)
    has_exclusion = any(keyword in clean_name for keyword in EXCLUSION_KEYWORDS)
    if has_exclusion:
        return None
        
    return clean_name

def clean_package_name(name):
    """
    Clean and standardize package names.
    Handles special cases like dots, underscores, etc.
    """
    name = str(name).lower()
    # Replace dots and underscores with hyphens
    name = name.replace('.', '-').replace('_', '-')
    # Remove any non-alphanumeric chars except hyphens
    name = ''.join(c for c in name if c.isalnum() or c == '-')
    # Remove duplicate hyphens
    while '--' in name:
        name = name.replace('--', '-')
    # Remove leading/trailing hyphens
    name = name.strip('-')
    return name

def clean_text(text):
    """Helper function to clean and standardize text fields"""
    # Convert to lowercase 
    text = str(text).lower()
    # Remove special characters but keep hyphens
    text = ''.join(c for c in text if c.isalnum() or c == '-' or c.isspace())
    # Remove extra spaces
    text = ' '.join(text.split())
    return text

def analyze_bioconda_packages(file_path):
    """
    Analyze Bioconda package data and print statistics.
    
    Args:
        file_path (str): Path to TSV file
    """
    df = load_package_data(file_path)
    
    # Clean package names
    df['Package'] = df['Package'].apply(clean_text)
    
    # Sort by update date, most recent first
    df_sorted = df.sort_values('Updated_Date', ascending=False)
    
    # Calculate statistics
    stats = {
        'Total packages': len(df),
        'Unique update dates': df['Updated_Date'].nunique(),
        'Bioconductor packages': df['Package'].str.startswith('bioconductor-').sum(),
        'R packages': df['Package'].str.startswith('r-').sum(),
        'Python packages': df['Package'].str.startswith('python-').sum(),
        'Perl packages': df['Package'].str.startswith('perl-').sum(),
        'Most recent update': df['Updated_Date'].max().strftime('%Y-%m-%d'),
        'Oldest update': df['Updated_Date'].min().strftime('%Y-%m-%d'),
        'Average desc length': df['Description'].str.len().mean()
    }
    
    # Print formatted statistics
    print("\n=== Package Statistics ===")
    for key, value in stats.items():
        if isinstance(value, float):
            print(f"{key}: {value:.1f}")
        else:
            print(f"{key}: {value}")
            
    print("\n=== Most Recently Updated Packages ===")
    print(df_sorted[['Package', 'Updated_Date']].head(10).to_string(index=False))
    
    print("\n=== Package Type Distribution ===")
    print("(Top 10 package prefixes)")
    prefixes = df['Package'].str.split('-').str[0].value_counts().head(10)
    for prefix, count in prefixes.items():
        print(f"{prefix}: {count}")

def clean_descriptions(input_file, output_file):
    """
    Clean package names and descriptions by standardizing case,
    removing special characters, and normalizing text.
    
    Args:
        input_file (str): Input TSV file path
        output_file (str): Output TSV file path  
    """
    try:
        df = load_package_data(input_file)
        
        # Clean descriptions only - package names already filtered
        df['Description'] = df['Description'].apply(clean_text)
        
        # Remove duplicates after filtering
        original_count = len(df)
        df = df.drop_duplicates(subset=['Package'])
        
        # Save cleaned data
        df.to_csv(output_file, sep='\t', index=False)
        
        # Print summary
        print(f"\nCleaning summary:")
        print(f"Original packages: {original_count}")
        print(f"After cleaning: {len(df)}")
        print(f"Removed duplicates: {original_count - len(df)}")
        print(f"Output saved to: {output_file}")
        
    except Exception as e:
        print(f"Error cleaning data: {str(e)}")
        sys.exit(1)

if __name__ == '__main__':
    input_file = "bioconda_filtered_packages.tsv"
    output_file = "bioconda_filtered_packages_clean.tsv"
    
    try:
        analyze_bioconda_packages(input_file)
        clean_descriptions(input_file, output_file)
    except KeyboardInterrupt:
        print("\nProcess interrupted by user")
        sys.exit(1)
