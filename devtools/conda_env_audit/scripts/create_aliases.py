# Portfolio-safe copy (python script)
# Identifiers generalized; secrets removed; paths parameterized.

import os

def create_aliases_bash_file(dependencies, env_dir="/home/<user>/miniconda3/envs"):
    bash_file_path = "generated_aliases.sh"
    
    with open(bash_file_path, "w") as bash_file:
        bash_file.write("echo "
")
        for dep in dependencies:
            alias_command = f"alias {dep}='{env_dir}/{dep}_env/bin/{dep}'\n"
            bash_file.write(alias_command)
        bash_file.write("\" >> ~/.bashrc && source ~/.bashrc\n")
    
    return bash_file_path

def execute_bash_file(bash_file_path):
    # os.system(f"bash {bash_file_path}")  # disabled for safety in portfolio version

def main():
    # Sample list of dependencies - in reality, you would extract these from the YAML files or bioconda list
    dependencies = [
        "cutadapt", "fastani", "graphviz", "scikit-learn", "hmmer",
        "iqtree", "diamond", "quast", "trnascan-se", "busco", "mafft", 
        "pplacer", "samtools", "kneaddata", "trim-galore", "megahit",
        "metabat2", "cd-hit", "anvio"
    ]
    
    # Step 1: Create the bash file with the aliases
    bash_file_path = create_aliases_bash_file(dependencies)
    
    # Step 2: Execute the bash file to apply the aliases
    execute_bash_file(bash_file_path)

if __name__ == "__main__":
    main()
