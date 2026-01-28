#ifndef SET_CLUSTER_IDS_H
#define SET_CLUSTER_IDS_H

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <utility>
#include <boost/functional/hash.hpp>

#include "../../readset.h"

void set_read_cluster_ids(ReadSet* read_set) {
    assert(read_set != nullptr);
    /**
     * Creating a hash map from (phaseset, haplotag) to cluster ID.
     * Each unique (phaseset, haplotag) pair gets a unique cluster ID which is the ID of the first read in that cluster.
     */
    std::unordered_map<std::pair<u_int32_t, bool>, uint32_t, boost::hash<std::pair<u_int32_t, bool>>> cluster_map;
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (read->hasPhaseSet() && read->hasHaplotag()) {
            std::pair<u_int32_t, bool> key = std::make_pair(read->getPhaseSet(), read->getHaplotag());
            if (cluster_map.find(key) == cluster_map.end()) {
                cluster_map[key] = read->getID();
            }
            read->setClusterID(cluster_map[key]);
            read->setClusterStatus(true);
        } else {
            // Untagged reads get their own read ID as cluster ID.
            read->setClusterID(read->getID());
            read->setClusterStatus(false);
        }
    }

    /**
     * Setting constrained cluster IDs for reads.
     */
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (read->getClusterStatus()) {
            uint32_t cluster_ps = read->getPhaseSet();
            bool cluster_hp = read->getHaplotag();
            // Find other haplotag in the same phaseset
            uint32_t constrained_cluster_id = cluster_map[std::make_pair(cluster_ps, !cluster_hp)];
            read->setConstrainedClusterID(constrained_cluster_id);
        }
        else {
            // untagged reads do not have constrained clusters
        }
    }
}


#endif // SET_CLUSTER_IDS_H