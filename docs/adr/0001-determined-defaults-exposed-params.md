# Determined defaults, exposed as params

Every tool setting and call threshold in the pipeline (Deacon, fastp, sylph, EsViritu and our own call cutoffs) has a **determined default** chosen by the team and recorded in the spec, and **every one is exposed as a param** in the YAML params file (validated by nf-schema). Nothing is hard-coded in a module's `script:` block or in `bin/`. We want one reproducible, agreed configuration that every weekly run uses unchanged. We also want a run to be re-tuned without a code change, as silly-fast-bfx does with its per-reference `build` / `sample_sketch` / `profile` blocks.

## Consequences

- **A default is a decision, not a suggestion.** Changing a default means changing the spec (and a new ADR if the reason is non-obvious). Overriding a default for one run happens in that run's params YAML, and the effective values are written into the run's outputs so results stay traceable.
- **Defaults live in one place.** They live in the params defaults (`nextflow.config` / the schema `default`), not repeated in modules. A module reads its settings from params and passes them through.
- **Current defaults:** the "Parameters and defaults" table in `docs/v1-spec.md` is the one list.
  (Read QC moved from fastp to NVD's BBTools chain on 2026-10-29; sensitive/rapid calls were renamed
  genome mapping-based / genome k-mer-based the same day.)
