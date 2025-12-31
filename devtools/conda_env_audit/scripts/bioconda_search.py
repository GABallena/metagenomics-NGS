# Portfolio-safe copy (python script)
# Identifiers generalized; secrets removed; paths parameterized.

import requests
from bs4 import BeautifulSoup
import time
from typing import List
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

BIOCONDA_REPO_URL = "https://anaconda.org/bioconda/repo"

def setup_session():
    session = requests.Session()
    retries = Retry(
        total=5, 
        backoff_factor=2,  # Gradual backoff
        status_forcelist=[524, 502, 503, 504], 
        allowed_methods=["GET"]
    )
    session.mount('http://', HTTPAdapter(max_retries=retries))
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

def filter_package_name(name: str) -> bool:
    """
    Specific filtering for package names, handling common prefixes.
    Returns True if package should be kept.
    """
    name = name.lower()
    
    # Explicitly exclude certain package prefixes/patterns
    exclude_prefixes = [
        'bioconductor-', 'bioconductor:', 'r-',
        'ucsc-', 'chipseq-', 'rnaseq-', 'scrnaseq-'
    ]
    
    if any(name.startswith(prefix) for prefix in exclude_prefixes):
        return False
        
    return True

def filter_text(text: str, keywords: List[str], exclusion_keywords: List[str]) -> bool:
    """
    Check if text matches inclusion/exclusion criteria.
    Returns True if text should be kept, False if it should be filtered out.
    """
    text = text.lower()
    
    # First check package name specific exclusions
    if not filter_package_name(text):
        return False
    
    # Then do regular keyword filtering
    has_inclusion = any(keyword.lower() in text for keyword in keywords)
    has_exclusion = any(ex_keyword.lower() in text for ex_keyword in exclusion_keywords)
    
    return has_inclusion and not has_exclusion

def fetch_package_details(session, keywords, exclusion_keywords, output_file):
    page = 1
    found_packages = False  # To detect when we've fetched at least one package

    while True:
        try:
            print(f"Searching page {page}...")
            response = session.get(f"{BIOCONDA_REPO_URL}?page={page}", timeout=30)
            response.raise_for_status()

            soup = BeautifulSoup(response.text, 'html.parser')
            package_table = soup.find('table')

            # Stop if no more packages are found on the page
            if not package_table:
                if found_packages:
                    print("No more packages found. Ending search.")
                else:
                    print("No packages found at all. Check if the URL or repository structure has changed.")
                break

            with open(output_file, "a") as f:
                for row in package_table.find_all('tr')[1:]:  # Skip header row
                    columns = row.find_all('td')
                    if len(columns) < 4:
                        continue  # Skip rows that do not have enough columns

                    package_name = columns[0].find('a').text.strip()
                    description = columns[2].text.strip() if len(columns) > 2 else "No description available"
                    updated_date = columns[3].text.strip() if len(columns) > 3 else "No date available"

                    print(f"Processing package: {package_name}")

                    # First check package name
                    if not filter_package_name(package_name):
                        print(f"Filtered out {package_name} - excluded prefix")
                        continue

                    # Then check content
                    name_ok = filter_text(package_name, keywords, exclusion_keywords)
                    desc_ok = filter_text(description, keywords, exclusion_keywords)
                    
                    if not (name_ok or desc_ok):
                        print(f"Filtered out {package_name} - no matching keywords")
                        continue

                    f.write(f"{package_name}\t{description}\t{updated_date}\n")
                    print(f"Added {package_name} to TSV file.")
                    found_packages = True  # Mark that we found a valid package
                    
            page += 1
            time.sleep(5)  # Avoid rate-limiting

        except requests.exceptions.Timeout:
            print(f"Timeout occurred for page {page}. Retrying...")
        except requests.exceptions.HTTPError as e:
            print(f"HTTP Error for page {page}: {e}")
            break
        except requests.exceptions.RequestException as e:
            print(f"Request failed: {e}")
            break


def search_and_write_package_details(keywords, exclusion_keywords, output_file):
    # Open and write the header of the TSV file
    with open(output_file, "w") as f:
        f.write("Package_Name\tDescription\tUpdated_Date\n")

    session = setup_session()
    fetch_package_details(session, keywords, exclusion_keywords, output_file)

# Example usage
keywords = [
    "phylo", "k-mer", "populat", "metagen", "antimicrob", "antibio", "resistance",
    "ortholog", "paralog", "Bayesian", "trim", "read", "mapp", "CARD", "AMR", "Illumina",
    "assembl", "cluster", "MAG", "genom", "pipeline", "alignment", "graph", "tensor",
    "informat", "algorithm", "geograph", "phage", "bootstrap", "statist", "pathogen",
    "quantif", "GPU", "loop", "GATK", "flanking", "slurm"
]

exclusion_keywords = [
    "RNA-seq", "Nanopore", "PacBio", "long-read", "cancer", "16S", "ITS", "single-cell",
    "PCR", "MinION", "server", "Windows", "plant", "mouse", "neural", "epigen", "soma",
    "chrom", "CRISPR", "expression", "cell-culture", "library", "API", "tissue", "MinION",
    "proteome", "LC-MS", "mass-spec", "microarray", "RNAi", "RNA-Seq", "RNAseq", "ChIP-Seq",
    "Whole-genome", "MS-MS", "Hi-C", "HiC", "HiSeq", "Hi-Seq", "HiChIP", "Hi-C", "HiC",
    "bioconductor-*", "bioconductor::", "bioconductor-", "bioconductor/", "bioconductor.",
    "UCSC", "CHIP-Seq", "ChIPSeq", "ChIP-seq", "ChIP", "ChIP-Seq", "ChIP-seq", "ChIPSeq",
    "cytometry", "methylat", "methylation", "methy", "methyseq", "methy-seq", "methy-Seq",
    "Genome wide", "10X", "10-X", "10x", "10-x", "10XGenomics", "10X-Genomics", "10X-Genomics",
    "hg19", "hg38", "mm10", "mm9", "GRCh37", "GRCh38", "GRCh37.p13", "GRCh38.p12", "GRCh38.p13",
    "GRCh37.p13", "GRCh38.p12", "GRCh38.p13", "GRCh37.p13", "GRCh38.p12", "GRCh38.p13",
    "GRCh37.p13", "GRCh38.p12", "GRCh38.p13", "GRCh37.p13", "GRCh38.p12", "GRCh38.p13",
    "GRCh37.p13", "GRCh38.p12", "GRCh38.p13", "GRCh37.p13", "GRCh38.p12", "GRCh38.p13",
    "single molecule", "single-molecule", "singlemolecule", "single-cell", "single cell",
    "singlecell", "single-cell", "single-cell", "singlecell", "single-cell", "single-cell",
    "RAD-Seq", "RADSeq", "RAD-seq", "RADseq", "RAD", "RAD-Seq", "RAD-seq", "RADSeq", "RAD",
    "methylation", "methy", "methyseq", "methy-seq", "methy-Seq", "methy-seq", "methyseq",
    "RNA-Seq", "RNAseq", "RNA-seq", "RNAseq", "RNA-Seq", "RNA-seq", "RNAseq", "RNA-seq",
    "genomic sequencing", "genomic-sequencing", "genomicsequencing", "genomic-sequencing",
    "WGS", "WES", "whole exome", "whole-exome", "wholeexome", "whole-exome", "whole-genome",
    "ONT", "nanopore", "nanopore sequencing", "nanopore-sequencing", "nanoporesequencing",
    "long read", "long-read", "longread", "long-read", "long-read", "longread", "long-read", 
    "long reads", "long-reads", "longreads", "long-reads", "long-reads", "longreads", "long-reads",
    "transcriptomic", "transcriptome", "transcriptomics", "HiFi", "Hi-Fi", "HiFi", "Hi-Fi",
    "HiFi sequencing", "HiFi-sequencing", "HiFisequencing", "HiFi-sequencing", "HiFi-sequencing",
    "mitochondria", "mitochondrial", "Bisulfite", "bisulfite", "bisulfite sequencing"
 

]

output_file = "bioconda_filtered_packages.tsv"

search_and_write_package_details(keywords, exclusion_keywords, output_file)
