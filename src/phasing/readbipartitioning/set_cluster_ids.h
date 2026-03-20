#ifndef SET_CLUSTER_IDS_H
#define SET_CLUSTER_IDS_H

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <functional>

#include "../../readset.h"

/**
 * A pair hash for PS + HP tags together.
 */
struct pair_hash;

/**
 * From the tags created for all the reads, this function creates a consistent tag
 * across the read set and assigns it as ClusterID.
 */
void set_read_cluster_ids(ReadSet* read_set);

#endif // SET_CLUSTER_IDS_H