#!/bin/bash
# Usage: ./find_homologs.sh <query_file> <subject_file> <output_file>
# Performs tblastn (protein query vs nucleotide subject), keeps hits with
# >30% identity AND alignment length >90% of query length, writes hits to output,
# and prints the number of matches to stdout.

set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <query_file> <subject_file> <output_file>" >&2
  exit 1
fi

query="$1"
subject="$2"
output="$3"

# Run tblastn with fields including qlen so we can filter per-query correctly.
tblastn -query "$query" -subject "$subject" -outfmt "6 qseqid sseqid pident length qlen" \
  | awk '$3 > 30 && ($4 / $5) > 0.9 {print}' > "$output"

# Print number of matches
wc -l < "$output"
