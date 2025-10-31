This folder contains scripts from whatshap (version 2.8) for implementing the phasing algorithm internally.

The following changes have been made:
    - Since we do not use pedigree in Giggles, that support has been removed. In some places, the variables have been hardcoded to keep things simple (TODO: remove them)
    - Since phasing will now be done iteratively with positions that are deemed to be clearly bi-allelic by the genotyping step, implementation has been update to reflect this.