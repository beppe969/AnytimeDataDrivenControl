# Reproduction contract

## Scientific inputs

This package generates the numerical results for the paper
**A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design**,
by **G. Calafiore**. Its scientific scripts are in `code/`, with comparison records
in `reference/`. The runner copies the scripts unchanged into a fresh workspace,
using their output-path convention `Path(__file__).resolve().parents[1]`.
The source-integrity manifest is checked before any numerical command starts.

The direct dependency list is recorded in `requirements.txt`. The additional
lock file records the dependency closure present in the tested environment.
It pins Python-package versions; it does not pin operating-system libraries or
provide cryptographic hashes for downloaded dependency wheels. Each run records
the actual interpreter and installed dependency versions.

## Isolation and failure handling

`full` copies no reference outputs into its workspace. `verify`, `figures` and
selected experiments copy reference data and original numerical audit records;
those copied inputs are enumerated with SHA-256 hashes in `run_manifest.json`.
Old run logs and aggregate verification reports are never seeded into a run.

Existing nonempty workspaces are rejected. The runner also rejects repository
roots, their ancestors and protected source/reference directories. Each child
process receives an isolated Matplotlib cache and Numba cache. Thread counts are
fixed to one. The `MPLBACKEND=Agg` setting permits headless figure generation.

Each script has a dedicated log. The first nonzero process exit fails the run;
subsequent steps are skipped. The manifest records completed steps and marks the
run as failed or interrupted. Python optimization (`-O`/`-OO`) is rejected because
the original verification scripts rely on assertions.

## Reference comparison

`tools/compare_reference.py` checks every retained numerical file in
`reference/data`, `reference/checks` and `reference/tables` against the full run.
It checks file presence before comparing content.

| Type | Rule |
| --- | --- |
| CSV | Same row/column structure; exact integer counts and text fields; numerical tolerance for floating-point columns |
| JSON | Same object structure; exact integers, booleans and text; numerical tolerance for floating-point values |
| NPZ | Same array names and exactly matching array values |
| Source LaTeX tables | Byte equality |
| Main figures | PDF/PNG presence; PNG pixel identity reported separately |

The default floating-point tolerances are `rtol=1e-9` and `atol=1e-11`.
Mechanical files use `rtol=1e-6` and `atol=1e-7` for floating-point solver outputs.
Discrete qualification indices and counts remain exact. JSON keys named
`elapsed_seconds` and `software` are excluded from numerical equivalence;
these values remain available in the underlying run records. Figure pixel
identity is informative and does not govern the numerical pass/fail decision.
The source script `verify_reported_results.py` additionally verifies the printed
paper values and the mechanical audit thresholds.

Cross-platform arithmetic can vary with BLAS/LAPACK implementations. A successful
full run with the matched versions is stronger evidence than a smoke check.
A successful `verify` run checks the reference results without generating the cohorts.
The latest concrete test evidence for this release is in
[validation/VALIDATION.md](../validation/VALIDATION.md).

## Figures and the certificate refinement

`final_figures.py` and `third_review_figures.py` generate the three figures used
in the paper. The numerical helpers also produce diagnostic figures. Output
PDF files can differ in creation metadata even when the images match exactly.
No TeX installation is needed because the plotting code uses
Matplotlib's own math-text renderer.

The persistence-window minimum described in the paper is a further
valid certificate option. The experiments here retain the original certificate.
In Figure 2, the global running minimum belongs to the archive of stored designs;
it is distinct from the minimum restricted to the current design's persistence
interval.

## Updating the repository

Keep `reference/` unchanged when investigating a discrepancy. New outputs belong
under `runs/`. Preserve the file-integrity metadata when refactoring code, and describe
scientific changes explicitly. Update baseline hashes only when adopting an
intentional new baseline; a new hash alone does not validate a scientific change.
To regenerate the release-wide checksum manifest after documentation or metadata
edits, run `python tools/build_release.py --manifest-only`.
