# Portfolio-safe copy: paths/identifiers generalized; databases not included.
import os
from time import sleep
import requests
from bs4 import BeautifulSoup

from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

# === CONFIG ===
snapshot_ts = "20211127085601"
base_url = f"https://web.archive.org/web/{snapshot_ts}/http://integrall.bio.ua.pt/?list"
n_pages = 240
debug = True

output_file = "unique_integrases.txt"

# === Setup retry-capable session ===
session = requests.Session()
retry = Retry(
    total=5,
    backoff_factor=1,
    status_forcelist=[500, 502, 503, 504],
    allowed_methods=["GET"],
)
adapter = HTTPAdapter(max_retries=retry)
session.mount("https://", adapter)
session.mount("http://", adapter)


def fetch_page(offset):
    url = f"{base_url}&s={offset}&ob=org"
    if debug:
        print(f"\n[DEBUG] Fetching: {url}")
    try:
        resp = session.get(url, timeout=10)
        resp.raise_for_status()
        return resp.text
    except Exception as e:
        print(f"[ERROR] Failed to fetch page at offset {offset}: {e}")
        return None


def parse_integrases(html):
    soup = BeautifulSoup(html, "html.parser")
    tables = soup.select('table')
    if len(tables) < 3:
        return []

    rows = tables[2].find_all("tr")[3:-1]
    integrases = []
    for i, row in enumerate(rows):
        cols = row.find_all("td")
        if len(cols) < 3:
            continue
        integrase = cols[2].text.strip()
        if debug:
            print(f"[DEBUG] Row {i}: Integrase = '{integrase}'")
        if integrase:
            integrases.append(integrase)
    return integrases


def main():
    seen = set()

    for page in range(n_pages):
        offset = page * 50
        html = fetch_page(offset)
        if not html:
            continue

        integrases = parse_integrases(html)
        seen.update(integrases)

    print(f"\n✅ Total unique integrase genes found: {len(seen)}")

    with open(output_file, "w") as out:
        for integrase in sorted(seen):
            out.write(f"{integrase}\n")

    print(f"📄 Written to: {output_file}")

    # Equivalent to: wc -l output.txt
    print(f"\n📊 Line count (wc -l): {len(seen)}")


if __name__ == "__main__":
    main()
