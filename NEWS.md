# LaTerra 0.0.0.9022 public wording repair (development; 2026-10-08)

- Made two C3 software-evidence messages require independent scientific
  adjudication without an internal reviewer name; scientific_pass stays FALSE.
- Removed the internal reviewer name from one README authority-admission
  sentence. Audit IDs, authority fields, frozen objects and scientific
  computations are unchanged. This is not a new scientific or release gate.

# LaTerra 0.0.0.9021 saved-result integrity repair (development; 2026-10-06)

- Preserve group and scope identity strings on saved-table read/report export,
  including distinct leading-zero keys and literal NA/TRUE keys. Reject
  unexpected internal-branch terminal identities instead of erasing them.
- Mark available plots from failed runs explicitly PARTIAL / FAILED with the
  failed stage, without filling absent dependencies or rerunning analyses.
- Check the loaded namespace's installed runtime before S2 computation/output
  creation. Unsupported load_all sessions receive installation guidance rather
  than failing at provenance after fitting; the runtime hash definition stays
  the installed R/LaTerra.rdb SHA-256.
- Document snapshot indexing versus original file identities. Existing
  scientific operators, writers and saved run files are unchanged.

# LaTerra 0.0.0.9020 OOF plotting repair (development; 2026-10-06)

- Fixed the saved OOF view losing prediction points when more than six
  supplied groups exceeded ggplot2's default shape palette. All groups now
  share a point shape; deterministic horizontal offsets retain grouping,
  and the subtitle reports eligible row and group counts.
- Added actual rendered-point checks for an eight-group public run and
  24, 128 and 512-group plotting fixtures. Scientific output files remain
  unchanged by plotting; group identities remain in editable plot data.

# LaTerra 0.0.0.9019 binary usability increment (development; 2026-10-06)

- Added a generic `lt_import()` entry retaining the complete accepted
  SplitAlignerR exchange import, plus explicit trait, grouping and ordered
  fold TSV readers. No recoding, key inference or fold generation is added.
- Added `lt_run_binary()` as an input-friendly route to the unchanged,
  explicitly selected historical S2 computation and validation recipes.
- Added recipe inspection and a saved-run reader with required-file and hash
  verification, concise result summaries, bounded D1-D4 diagnostics, reports
  and saved-output plots. Reporting does not refit or change the run.
- Added a runnable 24-tip, 32-gene synthetic example exported through the real
  producer. Its algebraic data support a software workflow demonstration,
  not biological, scientific-acceptance or performance claims.
- Positioned La Terra as a matrix-first framework and evolving platform;
  Marine/S2 remains its first exact-replay historical example. This increment
  does not claim release certification or add a future-work placeholder API.

# LaTerra 0.0.0.9018 (development; M1 review candidate)

- Added `lt_run_marine(data_root, output_dir)` for the frozen Marine/S2
  baseline-to-screen vertical slice. It computes from raw inputs without
  reading output oracles and retains outputs, stages, hashes and provenance.
- Admitted the pre-existing A2U full-state artifact only as supplemental
  upstream provenance under TARGET001-AUTHORITY-ADDENDUM001. Historical
  numerical authority and the partial classified table remain unchanged.
- Preserved the exact S2 arithmetic recipe, deterministic historical ASR and
  strict-prefix screen rule. No C3, null generator or inference substitution.
- Preserved the caller's RNG state around deterministic ASR, including an
  initially absent seed, and rejected unexpected consumption of seeded RNG.
- 9017 remains the certified source base. M1 computation is not full V1
  certification; M2 and the 22-target terminal review remain future stages.

# LaTerra 0.0.0.9017

- Removed the legacy finite-logGBI acknowledgement/seal mechanism. Ordinary
  La Terra inference now refuses `lt_logGBI_encoded` regardless of caller
  flags, provenance edits, configuration edits, or direct internal calls.
- Made every ordinary domain, association, calibration, adjustment, and
  internal score route fail closed on `lt_logGBI_encoded` with the dedicated
  `lt_error_legacy_encoded_inference_not_supported` condition.
- Split local structural validation from source-bound validation. The new
  `validate_lt_logGBI_encoded_source()` derives positive values, boundary
  state, axes, source identities, and full-source global-k10 coordinates from
  a separately supplied `lt_log_relative` authority.
- Added `validate_legacy_logGBI_replay()` as a statistic-free comparison to a
  separately frozen historical replay authority. Self-hashes remain accidental
  corruption checks and are not provenance authentication.
- Preserved historical 9016 evidence as structurally readable inert evidence;
  it cannot be reused to create a new association or calibration.

# LaTerra 0.0.0.9016

- Closed the configured-run bypass: every non-empty configuration attached to
  an `lt_run` is now validated by the authoritative configuration validator at
  construction, validation, serialization, and restoration boundaries.
- Replaced internally synthesized legacy encoded acknowledgement with one
  ephemeral authorization gate created only from an explicit public caller
  argument. Association and calibration provenance now records that exact
  authorization source.
- Enforced the mathematical and keyed identity invariants of the historical
  `legacy_global_k10` derived view, including source `m_ref`, positive-log
  minimum, sentinel, positive/boundary state, and unavailable-cell behavior.
- Added exact-key public subsetting for the derived legacy replay view while
  retaining the full source-domain `m_ref` and sentinel provenance.
- Added a fail-closed negative documentation scanner for accidental
  reintroduction of superseded finite-sentinel claims in current materials;
  explicitly versioned historical NEWS sections remain historical evidence.

# LaTerra 0.0.0.9015

- Closed every ordinary association-domain, univariate-association, and
  empirical-calibration route to historical finite-sentinel `lt_rate` objects
  whose representation is `logGBI`. Such objects remain readable and
  explicitly migratable, but ordinary inference now fails with classed
  guidance.
- Historical deserialized `lt_rate_spec(method = "logGBI")`, analysis specs,
  configuration requests, and direct execution are read/migrate-only and fail
  closed outside the explicit migration/replay APIs.
- Made `LaTerra_R_serialization_v3_raw_code_order_v1` executable for new
  boundary semantic and audit identities by explicitly serializing with R
  serialization version 3. Historical package-wide hashes are unchanged.
- Preserved original source-authority `m_ref`, positive-log minimum, and
  global-k10 sentinel metadata in explicit legacy views regardless of later
  analysis-domain filtering.
- Added fail-closed ratio-baseline assertions: positive observed Y with zero
  denominator is contradictory, while zero-denominator ratio cells are
  unavailable with `baseline_zero` and can never become infinity.
- Updated runtime documentation and schemas to distinguish the current
  `lt_log_relative` / `log_relative_boundary` contract from historical
  read/migrate-only finite `logGBI` objects.

# LaTerra 0.0.0.9013

- Added first-class externally supplied branch-state null ensembles. One
  replicate is one coherent branch-state world applied to all genes; La Terra
  does not generate or scientifically adjudicate the null worlds.
- Added streaming empirical calibration with inclusive ties, plus-one p-values,
  explicit non-estimability states, Monte-Carlo precision, and optional
  joint-null maximum-absolute-statistic summaries.
- Added separately invoked none/BH/BY/Holm/Bonferroni adjustment with explicit
  family identities and dependence-boundary documentation.
- Production ADD_LT, GBI, and zero-preserving logGBI values are consumed
  unchanged. The logGBI sentinel and zero mask are never recalculated from
  null traits.
- Serialized object writes now fail closed on existing destinations unless
  replacement is explicitly requested, and use same-directory atomic
  finalization where supported.

# LaTerra 0.0.0.9012

- Changed the V1 production meaning of `method = "logGBI"` to retain exact
  `Y = 0, mu > 0` cells through one global finite zero-state coordinate,
  `min(log(GBI[GBI > 0])) - 10`, plus a first-class `zero_encoded` mask.
- Added `method = "logGBI_strict"` for the former positive-only conditional
  log-relative geometry. Positive-cell numerical values are unchanged.
- Added zero-state provenance and trait-free diagnostic counts. The encoding
  is not a pseudocount and does not modify Y, C3 mu, or GBI.
- Legacy development objects whose representation was named `logGBI` but whose
  provenance identifies the old positive-only contract remain readable and
  are not silently reinterpreted. Recompute them with `logGBI_strict` to obtain
  a first-class current companion object, or with `logGBI` to opt into the new
  default zero-state semantics.

The development version skips 0.0.0.9011 because that number belongs to a
paused, unvalidated T2B draft and is not part of this source line.
