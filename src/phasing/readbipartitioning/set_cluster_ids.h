#ifndef SET_CLUSTER_IDS_H
#define SET_CLUSTER_IDS_H

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <boost/functional/hash.hpp>

#include "../../readset.h"

/**
 * From the tags created for all the reads, this function creates a consistent tag
 * across the read set and assigns it as ClusterID.
 */
void set_read_cluster_ids(ReadSet* read_set);

#endif // SET_CLUSTER_IDS_H