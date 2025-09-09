#!/bin/bash
# Usage: ./extract_spacers.sh <query file> <assembly file> <output file>

query=$1
assembly=$2
outfile=$3

# 1. Find repeats with blastn (query = repeat sequence, subject = genome assembly)
blastn -query "$query" -subject "$assembly" -outfmt "6 sseqid sstart send" |
    sort -k1,1 -k2,2n | \
    awk '{
        start = ($2 < $3 ? $2 : $3);
        end   = ($2 < $3 ? $3 : $2);
        print $1, start, end
    }' OFS="\t" > repeats.bed

# 2. Convert repeat BED to spacer BED (intervals between repeats)
awk 'NR>1 {print $1, prev_end, $2} {prev_end=$3}' OFS="\t" repeats.bed > spacers.bed

# 3. Extract spacers with seqtk
seqtk subseq "$assembly" spacers.bed > "$outfile"

# 4. Report number of spacers
grep -c "^>" "$outfile"

