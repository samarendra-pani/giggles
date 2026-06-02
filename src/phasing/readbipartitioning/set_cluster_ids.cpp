#include "set_cluster_ids.h"

struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        // Standard "hash_combine" logic
        return h1 ^ (h2 + 0x9e3779b9 + (h1 << 6) + (h1 >> 2));
    }
};

void set_read_cluster_ids(ReadSet* read_set) {
    assert(read_set != nullptr);
    /**
     * Creating a hash map from (phaseset, haplotag) to cluster ID.
     * Each unique (phaseset, haplotag) pair gets a unique cluster ID which is the ID of the first read in that cluster.
     */
    std::unordered_map<std::pair<uint32_t, bool>, uint32_t, pair_hash> cluster_map;
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (read->hasPhaseSet() && read->hasHaplotag()) {
            std::pair<uint32_t, bool> key = std::make_pair(read->getPhaseSet(), read->getHaplotag());
            auto [it, inserted] = cluster_map.try_emplace(key, read->getID());
            read->setClusterID(it->second);
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
        if (read->isClustered()) {
            uint32_t cluster_ps = read->getPhaseSet();
            bool cluster_hp = read->getHaplotag();
            // Find other haplotag in the same phaseset
            auto key = std::make_pair(cluster_ps, !cluster_hp);
            auto it = cluster_map.find(key);
            if (it != cluster_map.end()) {
                uint32_t constrained_cluster_id = it->second;
                read->setConstrainedClusterID(constrained_cluster_id);
            }
            else {
                // read might be clustered but the constrained cluster might not exist.
            }
        }
        else {
            // unclustered reads do not have constrained clusters
        }
    }
}