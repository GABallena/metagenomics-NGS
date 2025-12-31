# Portfolio-safe copy (shell script)
# Identifiers generalized; secrets removed; paths parameterized.

#!/bin/bash

# Define the absolute path to the Desktop directory
DESKTOP_DIR="/home/<user>/Desktop"

# List of directories to exclude
EXCLUDE_DIRS=("CAMISIM" "mhm2" "upcxx" "contig-extender" "SprayNPray" "PlasSuite")

# Get the current date and time for commit messages
CURRENT_DATETIME=$(date +"%Y-%m-%d %H:%M:%S")

# File size limit in MB (GitHub recommended)
SIZE_LIMIT_MB=100

# Function to display ignored file types from .gitignore
display_ignored_filetypes() {
    local repo_dir=$1
    if [ -f "$repo_dir/.gitignore" ]; then
        echo "Ignored file types in $(basename "$repo_dir"):"
        grep -E '^\*\.[a-zA-Z0-9]+' "$repo_dir/.gitignore" | sort | uniq
    else
        echo "No .gitignore file found in $(basename "$repo_dir")"
    fi
}

# Function to handle large files in the repository
handle_large_files() {
    local repo_dir=$1
    local large_files
    # Use git ls-files to respect .gitignore
    large_files=$(git ls-files | while read -r file; do
        if [ -f "$file" ]; then
            size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file")
            if [ "$size" -gt "$((SIZE_LIMIT_MB * 1024 * 1024))" ]; then
                echo "$file ($(numfmt --to=iec-i --suffix=B $size))"
            fi
        fi
    done)
    
    if [[ -n "$large_files" ]]; then
        echo "Large files detected (>${SIZE_LIMIT_MB}MB):"
        echo "$large_files"
    fi
}

# Function to process each Git repository
process_git_repo() {
    local repo_dir=$1
    echo -e "\n--- Processing Git repository in $(basename "$repo_dir") ---"
    cd "$repo_dir" || exit

    # Display ignored file types
    display_ignored_filetypes "$repo_dir"

    # Check for large files (only tracked files)
    handle_large_files "$repo_dir"

    # Show which files are being ignored
    echo "Files being ignored by .gitignore:"
    git status --ignored --porcelain | grep '^!!' | cut -c4-

    # Show untracked files that aren't ignored
    echo "Untracked files (not ignored):"
    git status --porcelain | grep '^??' | cut -c4-

    # Stage changes (this respects .gitignore by default)
    git add .

    # Double check what's actually staged
    echo "Files staged for commit:"
    git diff --name-only --cached

    # Check if there are changes to commit
    if git diff --staged --quiet; then
        echo "No changes to commit in $(basename "$repo_dir")"
        return
    fi

    # Commit changes with a message
    git commit -m "Automated commit on $CURRENT_DATETIME"

    # Push changes to the remote repository
    git push
}

# Main script logic
for dir in "$DESKTOP_DIR"/*/; do
    dirname=$(basename "$dir")

    # Skip excluded directories
    if [[ " ${EXCLUDE_DIRS[@]} " =~ " $dirname " ]]; then
        echo "Skipping $dirname"
        continue
    fi

    # Check if directory is a git repository
    if [ -d "$dir/.git" ]; then
        process_git_repo "$dir"
    fi
done
