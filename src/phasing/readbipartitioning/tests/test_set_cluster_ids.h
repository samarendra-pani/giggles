#ifndef TEST_SET_CLUSTER_IDS_H
#define TEST_SET_CLUSTER_IDS_H

#include "../set_cluster_ids.h"
#include "../haplotagcomputer.h"
#include "../phasesetcomputer.h"
#include "../../../tests_data.h"
#include <cassert>

/** Reads with multiple entries */
void test_multiposition_entries();

/** Reads with only one entry that is BLANK */
void test_singleposition_blankentries();

void test_set_cluster_ids();

#endif // TEST_SET_CLUSTER_IDS_H