# Portfolio-safe copy: paths/identifiers generalized; databases not included.
import requests
from bs4 import BeautifulSoup
import urllib3
import time

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

INPUT_FILE = "unique_isfinder_names.txt"
OUTPUT_FASTA = "ISfinder_sequences.fasta"
FAILED_LOG = "failed_is_names.txt"

with open(INPUT_FILE, "r") as infile:
    is_names = [line.strip() for line in infile if line.strip()]

with open(OUTPUT_FASTA, "w") as fasta_out, open(FAILED_LOG, "w") as failed_out:
    for i, is_name in enumerate(is_names, 1):
        url = f"https://isfinder.biotoul.fr/scripts/ficheIS.php?name={is_name}"
        try:
            r = requests.get(url, verify=False, timeout=15)
            soup = BeautifulSoup(r.text, 'html.parser')

            tag = soup.find("p", string="DNA sequence ")
            if not tag:
                print(f"[{i}/{len(is_names)}] ❌ Tag not found for {is_name}")
                failed_out.write(f"{is_name}\n")
                continue

            dna_div = tag.find_next("div")
            if not dna_div:
                print(f"[{i}/{len(is_names)}] ❌ DNA div not found for {is_name}")
                failed_out.write(f"{is_name}\n")
                continue

            dna_seq = dna_div.get_text().replace('\n', '')
            fasta_out.write(f">{is_name}\n{dna_seq}\n")
            print(f"[{i}/{len(is_names)}] ✅ {is_name}")
            time.sleep(0.5)

        except Exception as e:
            print(f"[{i}/{len(is_names)}] ❌ Error for {is_name} → {e}")
            failed_out.write(f"{is_name}\n")
