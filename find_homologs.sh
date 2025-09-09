#!/bin/bash
# Usage: ./find_homologs.sh <query_file> <subject_file> <output_file>

query="$1"
subject="$2"
output="$3"

# Run tblastn: protein query vs nucleotide subject
tblastn -query "$query" -subject "$subject" -outfmt "6 qseqid sseqid pident length" |
awk -v qlen=$(grep -v "^>" "$query" | tr -d '\n' | wc -c) \
    '$3 > 30 && $4/qlen > 0.9 {print}' > "$output"

# Print number of matches
wc -l < "$output"
