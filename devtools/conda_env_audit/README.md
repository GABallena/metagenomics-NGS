# Conda environment audit & update toolkit (portfolio-safe)

A small collection of scripts to:
- export conda environments to YAML (`env/*.yaml`)
- scan YAML dependencies and identify bioinformatics tools likely coming from Bioconda
- optionally check for newer versions (via Anaconda API / Bioconda pages)
- recreate environments from YAMLs

## Layout
- `scripts/`
  - `versioncheck.bash`: orchestration driver
  - `check_dependencies.py`: parse YAML envs and compare versions (via API)
  - `bioconda_search.py`: crawl Bioconda repo pages and filter packages by keywords
  - `analyze_packages.py`: post-process / summarize filtered Bioconda package TSV
  - `append_new.py`: update YAMLs with newer versions when detected
  - `update_yamls.bash`: rebuild envs from YAMLs
  - `install_individual_envs.sh`: create “one tool per env” installs (optional strategy)
  - `create_aliases.py`: generate a `generated_aliases.sh` file for convenience (portfolio version does **not** auto-modify `~/.bashrc`)
  - `git_status.bash`: bulk git status/add/commit/push helper (paths parameterized)
  - `checkingAPI.sh`: WakaTime heartbeat fetcher (**does not print API keys**; supports `WAKATIME_API_KEY` env var)

## Safety notes
- Secrets are never printed in the portfolio version.
- Repo/personal paths are generalized; you must adapt `DESKTOP_DIR`, `env_dir`, etc. for your machine.

## Minimal usage
```bash
cd devtools/conda_env_audit
bash scripts/versioncheck.bash
```
