# Marine/S2 historical replay guide

This guide records the inputs and outputs for the historical Marine replay.
For analyses with your own binary trait, see the [binary workflow guide](BINARY.md).

## Baseline and gene screens

The study-specific entry point now computes the frozen S2 arithmetic baseline
and GBI, deterministic branch annotations, and the marine/aquatic baseline
Welch screens from raw inputs. It writes outputs and an audited `lt_run`:

```r
library(LaTerra)
run <- lt_run_marine(data_root = "/path/to/marine_inputs",
                     output_dir = "/path/to/new_marine_run")
run$results$screens
```

Supply these three files beneath `data_root` (exact paths and snapshot hashes
are available through `lt_marine_profile()`):

```text
marine_inputs/
  00_traits_and_species/trait_table.mammal302.active_TY_NK_final_18pt.tsv
  17_large_matrix_archives/branch_coordinate_matrices.tar.gz
  supplemental_authority/marine_2026_fixed_coord_state_17432_v1.tsv.gz
```

The first two are the frozen Dryad inputs. The third is the pre-existing A2U
artifact admitted as TARGET001-AUTHORITY-ADDENDUM001. It supplements
state provenance only; the original raw values remain numerical authority,
and the historical 17,428-gene classified member retains its own partial axis.
No precomputed GBI or screen results are needed or read by this entry point.

ADDENDUM006 restores the historical execution path automatically: use the
supplied crosswalk as an inverse computational-order view, compute the unchanged
S2 recipe, write with the original R writer, reorder only the text tokens back
to the public axes, then read this newly generated GBI. Raw values and public
branch identities do not change. No manual conversion, Python runtime, archived
GBI input or precision option is involved. This historical order/serialization
rule does not apply to the ordinary other-data S2 entry point.

The output parent must exist and the run directory must be new and outside
the input root. No files are overwritten and no historical absolute-path
fallback is used. The default frozen environment is R 4.4.2, ape 5.8-1 and
castor 1.8.4 on aarch64-apple-darwin20. An explicit
`environment = "compatibility"` run is warned and non-certifying.

Major outputs: `raw_matrix.rds`, `upstream_tokens.rds`,
`S2_baseline_GBI.rds` (its `gbi` is the generated/read-back representation),
`GBI.generated_native.tsv`, `GBI.generated_consumed.oldlabels.tsv`,
`GBI_computational_order.tsv`, `GBI_consumption_boundary.rds`,
per-trait branch states/tested/significant/gene-ledger
tables, `run.rds`, `input_hashes.tsv`, `frozen_profile.yaml`, `sessionInfo.txt`,
`RUN_STATUS.yaml` and `SHA256SUMS`. Failures after output reservation retain
their actual stage in `RUN_STATUS.yaml`; a completed run is not automatically
a certified replay. The separate development validator compares results to
the frozen TARGET001 outputs after execution.

This historical study profile includes its frozen inclusive upper-2.5% trim
and strict-prefix FDR rule; neither is a recommended generic default. SHA
checks belong to this explicit frozen replay route, not normal user objects.
M1 alone does not execute the later models or confer full 22-target V1 certification.

## Full computation and reporting

The M2 entry extends the same raw-input workflow through permutation replay,
grouped validation, full-data fits, turnover and descriptive projections. It
requires the thirteen-file admitted bundle, whose layout is documented in
[the installed Marine workflow guide](MARINE_M2.md).

```r
run <- lt_run_marine_m2(
  data_root = "/path/to/marine_inputs",
  output_dir = "/path/to/results/marine_run_001",
  cores = 4L
)
```

`lt_export_marine_tables()` is a separate terminal reporting entry for a
completed run. It does not rerun the expensive scientific stages or import
old result workbooks. The five table files preserve their historical worksheet
schemas; separate companions retain authority, transformations, historical
references and the 22-target evidence ledger. Reporting-only dependencies are
the R packages `openxlsx`, `xml2` and `zip`.

```r
tables <- lt_export_marine_tables(
  data_root = "/path/to/marine_inputs",
  run_dir = "/path/to/results/marine_run_001",
  output_dir = "/path/to/results/marine_tables_001"
)
tables$run$results$target_status
```

Only the final S5 coefficient display uses the explicitly admitted six-
significant-digit numeric view. Fits, signs, markers and ordering retain full
precision. Original-generator-unlocated presentation rules are disclosed;
no old numeric answers are copied. TGT021 table reproduction remains under
independent development review. Neither a completed model run nor assembled
tables certify the development package.
