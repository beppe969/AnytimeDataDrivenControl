# Validation of the reproducibility package

**Release:** 1.0.1. **Validation date:** 11 September 2026.

This package generates the numerical results for the paper
**A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design**,
by **G. Calafiore**.

## End-to-end numerical validation

A complete `python reproduce.py full` execution passed in a new workspace.
It began with empty scientific-output directories and seeded **zero reference
outputs**. All 14 scientific-script executions passed. The two additional steps
exported the paper tables and compared regenerated outputs with `reference/`.

| Check | Result |
| --- | --- |
| Scientific input hashes | All 73 match the included manifest |
| Scientific Python files | All 16 unchanged |
| Reference files | All 55 unchanged |
| Numerical records | All 35 match the regenerated results |
| CSV records | All 24 are byte-identical |
| Stored NumPy arrays | All array names and values match exactly |
| Main PNG figures | All three are pixel-identical |
| Main PDF figures | All three render identically at 160 dpi |
| Numerical verification scripts | All five pass |
| Unit and regression tests | All 30 pass |
| Additional command modes | `smoke`, `verify` and `figures` pass |
| Full-run duration on this host | 137.36 seconds |

The run used CPython 3.13.5 and matched all packages in `requirements-lock.txt`.
Numerical thread counts were set to one. The recorded duration excludes
installation and applies to the validation host.

## Workloads and interpretation

The complete run regenerated the 10,000-stream first-target cohort and the
independent 10,000-stream continued-inspection cohort. It also regenerated the
8,000-observation disturbance study and re-solved the nominal tube-MPC checks.
The mechanical study used 800 design plants and 50,000 independent validation
plants. Random seeds and replication counts are unchanged.

The reported experiments use the original anytime certificate. The optional
persistence-window-minimum refinement is left separate. The mechanical example
retains its conditional support-reconstruction interpretation and its count of
34 contraction/LMI violations among 50,000 validation plants. The alternative
performance evaluator retains its exploratory status.

The numerical checks establish reproducibility of these computations. Finite
instance checks retain the qualifications stated in the paper's assumptions.

## Reproduce the checks

From the repository root with the reference Python environment installed:

```bash
python reproduce.py integrity --release
python -m unittest discover -s tests -v
python reproduce.py verify
python reproduce.py full --run-dir runs/paper
```

The new run directory must be empty or absent. Both checksum manifests are
included in the repository root and must be uploaded with the remaining files.
The release builder rejects missing required files.

## Scope of validation

Numerical execution was tested in the recorded Linux environment. A fresh
network installation and hosted GitHub Actions execution were not performed
during this validation. Both workflow definitions were parsed and their action
release pages were checked. Cross-platform execution remains a separate check.

## Evidence

| File or directory | Contents |
| --- | --- |
| `validation_summary.json` | Machine-readable validation summary |
| `full_run_manifest.json` | Environment, input hashes and execution results |
| `source_integrity.json` | Check of the 73 scientific input hashes |
| `reference_comparison.json` | Per-file numerical comparison |
| `pdf_render_comparison.json` | Independent rendered-PDF comparison |
| `checks/` | Numerical verification reports |
| `logs/` | Full-pipeline execution logs |
| `unit_tests.log` | Results of all 30 tests |
| `mode_checks.json` | Additional command-mode results |
| `metadata_checks.json` | Metadata and source consistency checks |
| `paper_results.md` | Recomputed paper tables and linked figures |
| `figures/` | Regenerated main figures |
| `tables/` | CSV and LaTeX table exports |
