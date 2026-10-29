# Rapid Wastewater Metagenomics

Weekly community-level viral surveillance from wastewater sequencing: a genome k-mer-based screen reports first, and a genome mapping-based pass follows.

## Calls

**Genome k-mer-based call**:
A detection of a taxon in a sample by matching the sample's k-mers against reference genomes (sylph). It reports first and is tuned for a very low false-positive rate.
_Avoid_: rapid call, hit, preliminary result, sylph result

**Genome mapping-based call**:
A detection of a taxon in a sample by mapping reads onto reference genomes (EsViritu). It reports after the k-mer-based call.
_Avoid_: sensitive call, confirmed call, validated hit, EsViritu result

**Species-confirmed call**:
A genome k-mer-based call whose species also received a genome mapping-based call in the same sample. It is worked out in the LIMS backend, not the pipeline; a call that is not species-confirmed is kept and labelled, never withdrawn.
_Avoid_: questionable call, false positive

**Inconsistency**:
A sample–taxon pair where the two kinds of call disagree: one detected it and the other did not.
_Avoid_: discrepancy, conflict, mismatch

## Reference and watchlist

**Reference subspecies label**:
VPD's own name for the reference genome a call landed on (for example "SARS-CoV-2_Delta"). It names the nearest reference, not a strain or variant detected in the sample.
_Avoid_: strain call, variant call, lineage detected

**VPD**:
EsViritu's Virus Pathogen Database. Both callers search against it, so calls from either tool resolve to the same taxonomy.
_Avoid_: viral DB, reference DB

**Watchlist**:
A user-supplied set of categories of pathogens of concern. Each category names taxa whose descendants are included, minus any excluded taxa.
_Avoid_: risk-group table, pathogen list, lookup TSV
