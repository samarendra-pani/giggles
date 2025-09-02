# this script is for taking gaf alignments and extracting reads from them
# it uses the following files:
# - GAF alignments
# - Read sequences which we get from gaftools find_path using the paths from the GAF alignments

BEGIN {
    FS="\t"
}

# This block now correctly executes ONLY for the first file (GAF alignments)
FNR == NR {
    start[FNR] = $8
    end[FNR] = $9
    name[FNR] = $1
    next # CRITICAL: Skip to the next line to avoid the block below
}

# This block now correctly executes ONLY for the second file (read sequences)
{
    print ">" name[FNR]
    print substr($0, start[FNR] + 1, end[FNR] - start[FNR] + 1)
}