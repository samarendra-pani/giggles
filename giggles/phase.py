import logging

from collections import defaultdict

from giggles import __version__
from giggles.core import (
    readselection
)
from giggles.graph import ComponentFinder

from giggles.utils import plural_s, warn_once

logger = logging.getLogger(__name__)

def select_reads(readset, max_coverage, preferred_source_ids):
    logger.info(
        "Reducing coverage to at most %dX by selecting most informative reads ...", max_coverage
    )
    selected_indices = readselection(readset, max_coverage, preferred_source_ids)
    selected_reads = readset.subset(selected_indices)
    logger.info(
        "Selected %d reads covering %d variants",
        len(selected_reads),
        len(selected_reads.get_positions()),
    )

    return selected_reads

def setup_families(samples, max_coverage):
    """
    Return families, family_trios pair.

    families maps a family representative to a list of family members

    family_trios maps a family representative to a list of trios in this family
    """

    # list of all trios across all families
    all_trios = dict()

    # Keep track of connected components (aka families) in the pedigree
    family_finder = ComponentFinder(samples)

    # map family representatives to lists of family members
    families = defaultdict(list)
    for sample in samples:
        families[family_finder.find(sample)].append(sample)

    # map family representatives to lists of trios for this family
    family_trios = defaultdict(list)
    for trio in all_trios:
        family_trios[family_finder.find(trio.child)].append(trio)
    logger.info(
        "Working on %d%s samples from %d famil%s",
        len(samples),
        plural_s(len(samples)),
        len(families),
        "y" if len(families) == 1 else "ies",
    )

    largest_trio_count = max([0] + [len(trio_list) for trio_list in family_trios.values()])
    if max_coverage + 2 * largest_trio_count > 23:
        logger.warning(
            "The maximum coverage is too high! "
            "Giggles may take a long time to finish and require a huge amount of memory."
        )
    return families, family_trios
