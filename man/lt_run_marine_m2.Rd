\name{lt_run_marine_m2}
\alias{lt_run_marine_m2}
\title{Execute the frozen Marine/S2 M1-M2 workflow}
\usage{
lt_run_marine_m2(data_root, output_dir,
  environment = c("frozen", "compatibility"), cores = 1L)
}
\arguments{
\item{data_root}{Frozen input bundle with M1 and admitted M2 supplements.}
\item{output_dir}{New output directory outside the input bundle. Never overwritten.}
\item{environment}{Reference R/ape/castor/glmnet stack, or explicitly warned non-certifying compatibility mode.}
\item{cores}{Positive integer for ordered forked validation batches. Windows requires 1. Seeds and fold order do not change.}
}
\description{
Computes S2 GBI from raw matrix inputs, deterministic trait annotation, fourteen
global screens, the exact historical 200 Fig5B worlds, nested/Phase11 grouped
validation, eight nested sensitivities, fourteen full-data fits, turnover and
descriptive projections. No archived fitted answers are runtime inputs.
}
\details{
ADDENDUM005 admits a source-bound literal historical screening map for six
Phase13 stage2 runs only. Exact input bytes, parser identities, full screening
vectors, terminal responses and original run ordinals are checked. Separate
screening/response/disagreement ledgers and the mandatory historical qualifier
are retained. Ordinary \code{lt_s2_gloocv}, including its Phase13 recipe, still
rejects terminal/branch disagreement. No public override or replacement hash
profile can activate the exception.

The only imported fitted-summary references are two explicitly admitted legacy
AUC cells, marked historical-reference-only. Phase11 comparisons are recomputed.
The historical Fig5B direction mismatch is preserved and disclosed. Full-data
fits and projections are descriptive, not held-out validation or ancestral
habitat inference. Scientific and numerical acceptance requires independent
TARGET001 replay review. Computation does not certify M2 or a release.

All stages run from the supplied input root. There is no implicit checkpoint
reuse, hidden historical path fallback, skipped required stage, automatically
generated replacement world, or caller-supplied substitute frozen hash profile.
Partial outputs and the failed stage are preserved on error. Completed stage
objects retain feature order, folds, fit diagnostics and preprocessing parameters.
The installed \file{workflow/MARINE_M2.md} describes the input layout and outputs.

This is an explicitly study-specific reproduction profile. Ordinary keyed user
data should use supported structural operators such as \code{lt_s2_gloocv}; no
Marine SHA identity is imposed on those objects.
}
\value{A completed \code{lt_run}. Tables, stage RDS objects and input, configuration,
software, gene, branch, seed, fold and output provenance are also written to disk.}
\seealso{\code{\link{lt_run_marine}}, \code{\link{lt_marine_profile}}, \code{\link{lt_s2_gloocv}}}
