# VPD test subset (7 assemblies, 18 records, 76 KB)

Real records from EsViritu's VPD v3.2.4 with the vector-only mask, extracted with
`samtools faidx` from `virus_pathogen_database.fna` plus the matching rows of
`virus_pathogen_database.all_metadata.tsv`. Used by the synthetic test data and the
end-to-end test, so tests need no network and no full VPD.

| Assembly | Species | Why it is here |
| --- | --- | --- |
| `GCF_000859645.1` | Tobamovirus capsici (PMMoV) | normalisation genus, both callers |
| `GCF_001461485.1` | Tobamovirus fructirugosum (ToBRFV) | second Tobamovirus, so the genus sum is tested |
| `acc:ON897592.1` | Betacoronavirus pandemicum (SARS-CoV-2) | ordinary both-callers call |
| `set:43d19931109e46…` | Rotavirus alphagastroenteritidis | 11 segments count as one genome |
| `PSVA/S2-2__Sadwavirus_alphananas…` | Sadwavirus alphananas | Assembly name contains `/` (sanitised genome file) |
| `set:MF678311` | Enterovirus B | k-mer-only case (short reads EsViritu drops) |
| `GCF_001595675.1` | Norovirus norwalkense | mapping-only case (diverged reads, genus-level call) |
