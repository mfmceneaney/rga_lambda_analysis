#!/bin/bash

# Default values
DEFAULT_URL="https://cernbox.cern.ch/s/m0SGWybeAxbfseN/download"
DEFAULT_FILENAME="Tree-ep-at-JLAB12-posHelicity.root"

# Print usage function
print_usage() {
    echo "Usage: $0 [URL] [OUTPUT_FILENAME]"
    echo
    echo "Downloads a .root file using wget from a CERNBox download link."
    echo
    echo "Arguments:"
    echo "  URL              Optional. CERNBox download link."
    echo "                   Default: $DEFAULT_URL"
    echo "  OUTPUT_FILENAME  Optional. Output filename for the downloaded .root file."
    echo "                   Default: $DEFAULT_FILENAME"
    echo
    echo "Examples:"
    echo "  $0"
    echo "  $0 \"https://cernbox.cern.ch/s/abc123/download\""
    echo "  $0 \"https://cernbox.cern.ch/s/abc123/download\" \"myfile.root\""
}

# If no arguments provided, print usage
if [ $# -eq 0 ]; then
    print_usage
    exit 1
fi

# Use provided arguments or fallback to defaults
URL="${1:-$DEFAULT_URL}"
FILENAME="${2:-$DEFAULT_FILENAME}"

# Print info
echo "Downloading ROOT file..."
echo "URL: $URL"
echo "Output: $FILENAME"

# Perform download
wget -O "$FILENAME" "$URL"

# Check result
if [[ $? -eq 0 ]]; then
    echo "Download successful!"
else
    echo "Download failed!"
    exit 1
fi
