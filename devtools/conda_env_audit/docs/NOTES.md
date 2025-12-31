## Practical notes

- The Bioconda repo HTML layout may change. If `bioconda_search.py` stops finding packages, update the parser logic.
- `check_dependencies.py` uses HTTP requests. Cache results or throttle to avoid rate limits.
- Prefer pinning versions conservatively in YAMLs unless reproducibility is critical.
