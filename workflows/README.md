## Modular Workflows

- **Preprocessing (`preprocessing.smk`)**  
  Read-level trimming/QC scaffold intended for inclusion in larger DAGs. Validates input layout, handles FastQC + Sickle trimming, and wires in downstream contaminant filtering and assembly stubs. Best used as a template to integrate into a project-specific config.

- **Phylogenetics (`phylogenetics.smk`)**  
  Minimal alignment/tree pipeline (MUSCLE → FastTree/MEGA-CC) expecting upstream clustering/BLAST artifacts. Provided as a lightweight example; extend with your own configs/scripts before production use.

- **Mobile Genetic Elements (`MGE.smk`)**  
  Sketch of an MGE discovery/annotation pipeline with configurable envs and directories. It demonstrates rule layout for plasmid prediction, MOB-suite, and GC skew analysis; wire it to your own inputs and environments before running.

Supporting configs live under `workflows/config/`. These files are anonymized and meant as starting points, not strict specifications.
