# Frozen Marine/S2 workflow (development, not release-certified)

ADDENDUM005 admits the six historical Phase13 screening-map routes with exact
source binding and mandatory mismatch disclosure. Their original responses and
folds are unchanged; ordinary user APIs still require coherent keyed states.
Computational availability is not TGT-009, M2 or release acceptance.

Use a new output directory, outside the input bundle:

```r
library(LaTerra)
run <- lt_run_marine_m2(
  data_root = "/path/to/marine_inputs",
  output_dir = "/path/to/results/marine_run_001",
  cores = 4L
)
```

Choose an explicit worker count. The same historical fold seeds are used with
one or multiple workers. A full run includes 12 grouped-validation jobs; it can
take substantially longer than the M1-only `lt_run_marine()` entry. No GPU,
Python runtime, historical absolute directory or manual intermediate fit is needed.

## Input bundle

Keep the original three M1 inputs, shown by `lt_marine_profile()`, together with
these ten approved M2 inputs. Do not put archived GBI, coefficients, predictions
or screen results into the computation. Exact snapshot checks belong only to
this explicitly frozen historical workflow.

```
marine_inputs/
  17_large_matrix_archives/branch_coordinate_matrices.tar.gz
  00_traits_and_species/trait_table.mammal302.active_TY_NK_final_18pt.tsv
  supplemental_authority/
    marine_2026_fixed_coord_state_17432_v1.tsv.gz
    trait_input.marine_binary_drop_sets.combined.tsv
    trait_input.aquatic_v2_drop_sets.combined.tsv
    marine_NC_submitted_turnover_module_lookup_131x4.tsv
    marine_NC_submitted_legacy_AUC_reference_2x2.tsv
    phase13/
      endpointfix_stage2_deterministic_manifest.tsv
      mammal302.anno.nwk
      mammal302.anno.BL_support.nwk
      mammal.branch.txt
  10_permutation_controls/Fig5B_positive_count_matched_permutation/
    endpointfix_permutation_label_sets_long.tsv
  14_supplementary_tables/TableS5_tsv_exports/
    Supplementary_Table_S5_Figure4C_predictor_annotation.tsv
```

The complete A2U state table supplements provenance; it does not replace the
immutable historical partial classified table or modify raw numerical input.
The two legacy AUC cells are reporting-only references. All actual nested and
Phase11 comparator predictions are computed by this run.

## Read the results

- `M1/`: raw provenance, S2 baseline/GBI, baseline branch states and screens.
- `screens/`: all fourteen explicitly admitted global screens and states.
- `phase13_authority/`: separate complete literal screening, physical-state and
  terminal-response ledgers for six stage2 runs, with source and mismatch audits.
  Historical absolute paths in the manifest are inert provenance, not fallbacks.
- `Fig5B*`: literal historical permutation replay, with its direction qualifier.
- `nested_*`, `phase11_*`, `Fig5A*`: grouped validation, computed auxiliary
  comparisons, complete per-fold objects, and frozen sensitivity tables.
- `full_*`, `display_*`: full-data architecture and admitted annotation joins.
- `turnover_*`: observed sets, exact ordered universes and historical null stream.
- `internal_*`, `focal_projections.tsv`, `nc_*`: descriptive projections and
  source-table assembly. These are not held-out or ancestral-habitat inference.
- `run.rds`, `RUN_STATUS.yaml`, `stages.tsv`, `*_ledger.tsv`, `*_hashes.tsv`,
  `SHA256SUMS`, `QUALIFIERS.txt`: execution status and provenance.

`COMPUTED_M2_NOT_REPLAY_CERTIFIED` means computation completed, not that scientific
review passed. Historical optimizer warnings remain visible and are retained in
fitted objects. Failure preserves partial evidence and reports the failed stage;
the program does not silently skip that stage or pick a replacement model.

## Export the final table views without fitting again

After the scientific run completes, use the same input root and a separate new
report directory. Install the optional R packages `openxlsx`, `xml2` and `zip`
if you want XLSX export; they are not numerical-analysis dependencies.

```r
tables <- lt_export_marine_tables(
  data_root = "/path/to/marine_inputs",
  run_dir = "/path/to/results/marine_run_001",
  output_dir = "/path/to/results/marine_tables_001"
)
tables$files
tables$run$results$target_status
```

This produces Tables S1, S2, S4, S5 and S6 with their eighteen historical
worksheet schemas. It reads fresh output objects and the already admitted
reporting metadata, not archived workbooks. Every dynamic row must already exist
in this run before the final worksheet ordering is applied. It never refits,
regenerates worlds, changes an input, overwrites a scientific run or silently
substitutes a legacy result. A report-only failure can be retried in a **new**
report directory after fixing that reporting issue; no model rerun is needed.

The final S5 E:G view uses the separately admitted six-significant-digit numeric
presentation only after full-precision selection, marker/color and list/order
derivation. These worksheet numbers are not replacement scientific coefficients.
The original generating expression was not recovered; the now-frozen rule and
its exact authority are disclosed in the companions. S4's six historical
discussion-reference rows remain explicitly imported, not recomputed nulls.

`run.rds` combines the source-scoped scientific ledgers with the reporting
version/environment/configuration. `target_status.tsv` and
`target_evidence_ledger.tsv` link all 22 required targets to hashed outputs.
`logical_sheets.rds` is the **final worksheet view**, not a fitting input;
`reporting_provenance.rds`, `REPORTING_QUALIFIERS.txt` and `SHA256SUMS` retain
the complete report audit trail. Keep the scientific directory alongside the
reports: the report does not duplicate the large fitted objects. Actual
dependency roots are recorded, and relocated replay uses its own real paths.

All 22 computation statuses remain `COMPUTED_NOT_REPLAY_VALIDATED` until separate
output-only comparison and independent acceptance. The exporter cannot grant a
TGT019/TGT021, milestone, scientific or release PASS to itself.

## Other data and traits

`lt_run_s2()` accepts a coordinate-valid branch-length `lt_matrix`, a supplied
annotation tree/coordinate ledger, an explicit binary `lt_trait`, terminal groups
and an ordered fold ledger. It runs the explicitly selected historical S2 recipe
from raw input through GBI, ASR, screen, grouped fitting, outputs and provenance.
It does not require Marine hashes or dimensions. `lt_s2_gloocv()` also remains
available when the user already has a GBI matrix and branch states.

The S2 recipe must be chosen explicitly; its historical trimming is not a
recommended generic QC policy. Continuous traits, a general null generator,
automatic recoding and inferred grouping/folds are not supplied by these entries.
This is an initial supported other-data path, not a claim that every future
La Terra workflow or first-priority UX requirement is complete.
