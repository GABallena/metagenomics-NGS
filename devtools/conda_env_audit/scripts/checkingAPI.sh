# Portfolio-safe copy (shell script)
# Identifiers generalized; secrets removed; paths parameterized.

#!/bin/bash

# Fetch WakaTime API key from ~/.wakatime.cfg
echo "Fetching WakaTime API key..."
api_key="${WAKATIME_API_KEY:-$(grep 'api_key' ~/.wakatime.cfg | awk '{print $3}')}"
# Debugging: Echo the fetched API key

# Check if API key is present
if [ -z "$api_key" ]; then
    echo "Error: API key not found!"
    exit 1
else
    echo "API key found!"
fi

# Base64 encode the API key for Basic Authentication
encoded_key=$(echo -n "$api_key" | base64)

# Get the current date in YYYY-MM-DD format
current_date=$(date +%Y-%m-%d)

# Perform API request to WakaTime using Basic Auth, adding the current date parameter
echo "Performing API request to WakaTime..."
wakatime_stats=$(curl -s -H "Authorization: Basic $encoded_key" \
    "https://wakatime.com/api/v1/users/current/heartbeats?date=$current_date")

# Check if the response contains an error
if echo "$wakatime_stats" | grep -q '"error"'; then
    echo "Error in API response: $wakatime_stats"
    exit 1
else
    echo "API response received!"
fi

# Create the output TSV file
output_file="wakatime_project_stats.tsv"

# Check if a TSV file already exists, if not, create it with headers
if [ ! -f "$output_file" ]; then
    echo -e "Project\tLanguage\tLines\tCategory\tCreated At" > $output_file
fi

# Loop through each heartbeat and extract relevant fields
echo "Extracting data for the TSV file..."
for entry in $(echo "$wakatime_stats" | jq -c '.data[]'); do
    project=$(echo "$entry" | jq -r '.project')
    language=$(echo "$entry" | jq -r '.language')
    lines=$(echo "$entry" | jq -r '.lines')
    category=$(echo "$entry" | jq -r '.category')
    created_at=$(echo "$entry" | jq -r '.created_at')

    # Ensure all fields have values
    project=${project:-"N/A"}
    language=${language:-"N/A"}
    lines=${lines:-"N/A"}
    category=${category:-"N/A"}
    created_at=${created_at:-"N/A"}

    # Append to the TSV file
    echo -e "$project\t$language\t$lines\t$category\t$created_at" >> $output_file
done

echo "Data extraction complete. Saved to $output_file."
