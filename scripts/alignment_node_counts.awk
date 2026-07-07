#!/usr/bin/awk -f

# Set the field separator to < or >
BEGIN {
    FS = "[<>]"
}

{
    nodes_in_alignment = 0
    
    # Count valid nodes in the current alignment (line)
    for(i=1; i<=NF; i++) {
        if ($i != "") {
            nodes_in_alignment++
        }
    }
    
    # Only process lines that actually had nodes
    if (nodes_in_alignment > 0) {
        counts[nodes_in_alignment]++
        
        # Track the maximum number of nodes seen in a single alignment
        if (nodes_in_alignment > max_nodes) {
            max_nodes = nodes_in_alignment
        }
    }
}

END {
    # Print the distribution from 1 to max_nodes
    for (x = 1; x <= max_nodes; x++) {
        print "Number of alignments with " x " nodes: " (counts[x] ? counts[x] : 0)
    }
}