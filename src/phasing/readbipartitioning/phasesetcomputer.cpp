#include "phasesetcomputer.h"

void find_phasesets_tag_reads(ComponentFinder<uint32_t>* component_finder, const std::unordered_set<uint32_t>* accessible_positions_set, ReadSet* read_set, const std::unordered_set<uint32_t>& heterozygous_positions) {
    
    bool filter_by_het = !heterozygous_positions.empty();
    
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (!read->isSelected()) {
            continue; // Skip unselected reads.
        }
        std::vector<uint32_t> read_positions;
        read_positions.reserve(read->getVariantCount());
        for (uint32_t j = 0; j < read->getVariantCount(); ++j) {
            uint32_t pos = read->getPosition(j);

            // Check if it is a phased position
            if (accessible_positions_set->find(pos) == accessible_positions_set->end()) { continue; }
            // Check heterozygous constraint
            // If filter_by_het is false, we skip this check
            if (filter_by_het && heterozygous_positions.find(pos) == heterozygous_positions.end()) { continue; }
            read_positions.push_back(pos);
        }

        // Merge components for all pairs of positions in the read
        if (read_positions.size() > 1) {
            uint32_t first = read_positions[0];
            for (size_t m = 1; m < read_positions.size(); ++m) {
                component_finder->merge(first, read_positions[m]);
            }
        }
    }

    std::cerr << "[Core::Phasing] Found " << component_finder->size()-(accessible_positions_set->size()-heterozygous_positions.size()) << " phaseblocks (from " << heterozygous_positions.size() << " heterozygous positions), and " << accessible_positions_set->size()-heterozygous_positions.size() << " positions are homozygous variants (can't be phased)." << std::endl;

    // tagging reads with their phaseset ID (representative position)
    size_t count1 = 0;
    size_t count2 = 0;
    for (uint32_t i = 0; i < read_set->size(); ++i) {
        Read* read = read_set->get(i);
        if (read->isSelected()) {
            // finding first heterozygous position in the read
            uint32_t first_het_pos = 0;
            for (uint32_t j = 0; j < read->getVariantCount(); ++j) {
                uint32_t pos = read->getPosition(j);
                if (heterozygous_positions.find(pos) != heterozygous_positions.end()) {
                    first_het_pos = pos;
                    break;
                }
            }
            if (first_het_pos == 0) {
                // no heterozygous position found in the read
                continue;
            }
            uint32_t ps = component_finder->find(first_het_pos);
            count1++;
            read->setPhaseSet(ps);
            continue;
        }
        else {
            // finding heterozygous positions in the read and checking their components
            bool all_same_component = true;
            bool has_het_position = false;
            uint32_t rep = 0;
            for (uint32_t j = 0; j < read->getVariantCount(); ++j) {
                uint32_t pos = read->getPosition(j);
                if (heterozygous_positions.find(pos) != heterozygous_positions.end())
                {
                    has_het_position = true;
                    uint32_t comp = component_finder->find(pos);
                    if (rep == 0) {
                        rep = comp;
                    }
                    else if (comp != rep) {
                        all_same_component = false;
                        break;
                    }
                }
                else {
                }
            }
            if (all_same_component && has_het_position) {
                count2++;
                read->setPhaseSet(rep);
            }
        }  
    }
    std::cerr << "[Core::Phasing] Out of " << read_set->size() << " total reads, " << count1 << " selected reads and " << count2 << " unselected reads have phasesets. " << read_set->size()-count1-count2 << " don't have a phaseset." << std::endl;
}

void compute_phasesets(const std::vector<uint32_t>* accessible_positions, ReadSet* read_set, ReadSet* superreads) {
    
    std::unordered_set<uint32_t> heterozygous_positions;
    std::unordered_set<uint32_t> accessible_positions_set(accessible_positions->begin(), accessible_positions->end());

    assert (superreads->size() == 2); // two superreads represeting the two haplotypes
    Read* superread0 = superreads->get(0);
    Read* superread1 = superreads->get(1);

    assert (superread0->getVariantCount() == superread1->getVariantCount()); // both superreads should have same number of variants

    for (uint32_t i = 0; i < superread0->getVariantCount(); ++i) {
        assert (superread0->getPosition(i) == superread1->getPosition(i)); // both superreads should have variants at same positions
        // skip positions that are not accessible (that have more than 2 active alleles)
        if (accessible_positions_set.find(superread0->getPosition(i)) == accessible_positions_set.end()) { continue; }
        
        Entry::allele_t allele0 = superread0->getEntry(i)->get_allele_type();
        Entry::allele_t allele1 = superread1->getEntry(i)->get_allele_type();

        if ((allele0 == Entry::ALLELE1 && allele1 == Entry::ALLELE2) ||
            (allele0 == Entry::ALLELE2 && allele1 == Entry::ALLELE1)) {
            heterozygous_positions.insert(superread0->getPosition(i));
        }
    }
    ComponentFinder<uint32_t> component_finder(*accessible_positions);
    find_phasesets_tag_reads(&component_finder, &accessible_positions_set, read_set, heterozygous_positions);
}