# Portfolio-safe copy (python script)
# Identifiers generalized; secrets removed; paths parameterized.

import requests
import yaml
import os

# Directory containing the YAML files
yaml_dir = "env/"

# Function to check if the package exists on Bioconda and fetch its latest version
def check_bioconda_package(dependency, version=None):
    package_url = f"https://bioconda.github.io/recipes/{dependency}/README.html#package-{dependency}"
    print(f"Checking {package_url}")
    
    response = requests.get(package_url)
    
    # If the page exists (status code 200), consider it a Bioconda package
    if response.status_code == 200:
        print(f"{dependency} is a Bioconda package.")
        
        # Check version if provided
        latest_version = fetch_bioconda_version(dependency)
        if version:
            print(f"Installed version: {version}, Latest Bioconda version: {latest_version}")
            if version == latest_version:
                print(f"{dependency} is up-to-date.")
                return False  # No update needed
            else:
                print(f"{dependency} has a newer version available: {latest_version}")
                return latest_version  # Return the new version
        return None
    else:
        print(f"{dependency} is NOT a Bioconda package.")
        return None

# Function to fetch the latest version of a Bioconda package using the Anaconda API
def fetch_bioconda_version(package_name):
    api_url = f"https://api.anaconda.org/package/bioconda/{package_name}"
    response = requests.get(api_url)
    
    if response.status_code == 200:
        package_info = response.json()
        latest_version = package_info['latest_version']
        return latest_version
    else:
        print(f"Failed to fetch version info for {package_name}.")
        return None

# Function to collect unique dependencies from the YAML files
def collect_unique_dependencies(yaml_dir):
    unique_dependencies = set()
    
    for yaml_file in os.listdir(yaml_dir):
        if yaml_file.endswith(".yaml"):
            print(f"Processing {yaml_file}...")
            dependencies = parse_yaml_for_dependencies(os.path.join(yaml_dir, yaml_file))
            unique_dependencies.update(dependencies)
    
    return unique_dependencies

# Function to parse the dependencies from a YAML file, extracting both package names and versions
def parse_yaml_for_dependencies(yaml_file):
    with open(yaml_file) as f:
        environment = yaml.safe_load(f)
    
    dependencies = environment.get('dependencies', [])
    parsed_dependencies = []
    
    for dep in dependencies:
        if isinstance(dep, str):
            # Handle simple string dependencies like "package=1.0"
            dep_parts = dep.split('=')
            package_name = dep_parts[0]
            package_version = dep_parts[1] if len(dep_parts) > 1 else None
            parsed_dependencies.append((package_name, package_version))
        elif isinstance(dep, dict):
            # Handle dictionaries (sometimes conda has complex dependencies)
            for key in dep.keys():
                parsed_dependencies.append((key, None))
    
    return parsed_dependencies

# Function to update the YAML files with newer versions
def update_yaml_with_new_versions(yaml_file, updated_packages):
    with open(yaml_file, 'r') as f:
        environment = yaml.safe_load(f)
    
    for i, dep in enumerate(environment['dependencies']):
        if isinstance(dep, str):
            package_name = dep.split('=')[0]
            if package_name in updated_packages:
                new_version = updated_packages[package_name]
                environment['dependencies'][i] = f"{package_name}={new_version}"

    with open(yaml_file, 'w') as f:
        yaml.dump(environment, f)

# Main function to check if each dependency is in Bioconda, compare versions, and update YAML if necessary
def main():
    unique_dependencies = collect_unique_dependencies(yaml_dir)
    updated_packages = {}
    
    # Prepare the output file
    for yaml_file in os.listdir(yaml_dir):
        if yaml_file.endswith(".yaml"):
            updated_packages = {}
            for dependency, version in parse_yaml_for_dependencies(os.path.join(yaml_dir, yaml_file)):
                latest_version = check_bioconda_package(dependency, version)
                if latest_version:
                    updated_packages[dependency] = latest_version
            if updated_packages:
                print(f"Updating {yaml_file} with newer versions: {updated_packages}")
                update_yaml_with_new_versions(os.path.join(yaml_dir, yaml_file), updated_packages)
    
    print("Finished updating YAML files with newer Bioconda versions.")

if __name__ == "__main__":
    main()
