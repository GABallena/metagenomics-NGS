# Portfolio-safe copy (python script)
# Identifiers generalized; secrets removed; paths parameterized.

import requests
import os
import yaml

BIOCONDA_BASE_URL = "https://anaconda.org/bioconda/"

bioinformatics_keywords = ["bioconda"]

def is_bioconda_package(package_name):
    # Check if the package is listed in bioconda
    url = f"{BIOCONDA_BASE_URL}{package_name}"
    response = requests.get(url)
    
    if response.status_code == 200 and any(keyword in response.text for keyword in bioinformatics_keywords):
        return True  # Package is a bioinformatics tool in bioconda
    else:
        return False  # Package is not found in bioconda

def extract_dependencies_from_yaml(env_folder):
    dependencies = set()

    # List all YAML files in the environment directory
    yaml_files = [f for f in os.listdir(env_folder) if f.endswith(".yaml")]

    for yaml_file in yaml_files:
        with open(os.path.join(env_folder, yaml_file), "r") as file:
            env_data = yaml.safe_load(file)
            # Check for the dependencies section
            if "dependencies" in env_data:
                for dep in env_data["dependencies"]:
                    if isinstance(dep, str):
                        # If it's a string, extract package name
                        dep_name = dep.split("=")[0]
                        dependencies.add(dep_name)
                    elif isinstance(dep, dict):
                        # If it's a dictionary (like pip packages), extract pip deps
                        if "pip" in dep:
                            for pip_dep in dep["pip"]:
                                dep_name = pip_dep.split("=")[0]
                                dependencies.add(dep_name)
    return dependencies

def filter_bioconda_dependencies(dependencies):
    bioinformatics_tools = []
    for dep in dependencies:
        if is_bioconda_package(dep):
            bioinformatics_tools.append(dep)
    return bioinformatics_tools

if __name__ == "__main__":
    env_folder = "env/"  # Adjust this to point to the correct directory
    dependencies = extract_dependencies_from_yaml(env_folder)
    bioinformatics_tools = filter_bioconda_dependencies(dependencies)
    
    # Print or write the results to a file
    with open("bioinformatics_tools.txt", "w") as f:
        for tool in bioinformatics_tools:
            f.write(f"{tool}\n")
    print(f"Bioinformatics tools saved to bioinformatics_tools.txt")
