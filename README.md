# Anytime Control — Reproducibility Package

This package generates the numerical results for the paper
**A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design**,
by **G. Calafiore**.

It reproduces the numerical studies in Section 6, including Tables 2–4 and
Figures 1–3. The repository contains 16 scientific Python files and the
reference outputs used for numerical comparison. An execution runner and table
exporters provide isolated workspaces with recorded verification results.

**Start here:** [experiment map](docs/EXPERIMENTS.md) ·
[protocols and seeds](docs/PROTOCOLS.md) ·
[validation report](validation/VALIDATION.md) ·
[data dictionary](docs/DATA_DICTIONARY.md).

## Install

The reference environment is **CPython 3.13.5**. No MATLAB, commercial solver or
TeX installation is needed to execute the numerical studies or generate figures.
The mechanical convex program uses SciPy's SLSQP implementation with its original
instance-level checks.

From the repository root, with Python 3.13.5 available:

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-lock.txt
python reproduce.py doctor
```

On Windows, create the environment with `py -3.13 -m venv .venv` and use
`.venv\Scripts\Activate.ps1` in PowerShell. Check the installed Python patch
version with `python --version`.

`requirements.txt` lists the seven direct dependencies.
`requirements-lock.txt` additionally pins the transitive Python packages in the
tested environment. The runner checks these versions and the
interpreter before execution. `--allow-environment-drift` permits a different
environment and records that choice; exact replication should use the reference
versions.

## Reproduce

```bash
# Validate installation with a small QP path and scalar certificate checks.
python reproduce.py smoke

# Recheck the archived paper results using all five numerical verification scripts.
python reproduce.py verify

# Recreate the paper figures and tables from the reference records.
python reproduce.py figures

# Recompute all numerical studies from empty output directories.
python reproduce.py full
```

Each command prints its output directory under `runs/`. Supply a new directory
to obtain a predictable location:

```bash
python reproduce.py full --run-dir runs/paper
```

A full run retains the original seeds and all replication counts. It regenerates
the two 10,000-stream cohorts, the disturbance-set experiments and the mechanical
validation on 50,000 fresh plants. It then checks the reported results and
compares regenerated numerical files against `reference/`. The additional diagnostic
studies documented in the experiment map are also included.

| Command | Uses reference outputs as inputs? | What it establishes |
| --- | --- | --- |
| `smoke` | No | Small execution and implementation checks |
| `verify` | Yes | Consistency of the archived results with the numerical assertions |
| `figures` | Yes | Regeneration of Figures 1–3 and Tables 2–4; also re-solves FH tube-MPC checks |
| `full` | No | End-to-end numerical regeneration followed by reference comparison |
| `experiment NAME` | Yes, for prerequisites | Regeneration of a selected study group |

Existing run directories are never overwritten. Reference outputs remain
unchanged. The smoke test uses a reduced workload and makes no claim to reproduce
the paper's Monte Carlo statistics. The `figures` command also creates several
additional diagnostic plots.

### Selected studies

```bash
python reproduce.py experiment controller-mpc
python reproduce.py experiment stopping-archive
python reproduce.py experiment monitoring
python reproduce.py experiment horizon-comparison
python reproduce.py experiment mechanical
```

The study groups follow the original script boundaries. In particular,
`horizon-comparison` includes the FH/summable cohort comparisons and their replay
checks. See [the experiment map](docs/EXPERIMENTS.md) for exact dependencies.

### Inspect outputs

A full workspace contains:

```text
runs/paper/
├── data/                 # Regenerated stream records and numerical summaries
├── figures/              # PDF and PNG figures
├── tables/               # Source tables plus table2/3/4 in CSV and LaTeX
├── checks/               # Numerical audits and reference_comparison.json
├── logs/                 # Separate stdout/stderr log for every execution step
├── paper_results.md      # Recomputed paper tables and linked figures
└── run_manifest.json     # Environment, input provenance, exit codes and hashes
```

For a second comparison or an integrity check:

```bash
python reproduce.py compare --run-dir runs/paper
python reproduce.py integrity
python reproduce.py integrity --release
python -m unittest discover -s tests -v
```

`SOURCE_MANIFEST.sha256` checks the scientific code and reference inputs.
`SHA256SUMS` checks the complete prepared release, excluding itself. A full
comparison treats integer counts as exact and reports floating-point tolerances
explicitly. Execution times are excluded from numerical comparison. PNG pixel
identity is reported separately from numerical equivalence because rendering can
vary across platforms.

## Scientific interpretation

The experiments use the paper's **original anytime certificate**. The later
persistence-window running-minimum refinement is not applied to these results.
The certified-archive figure tracks the best certificate among stored designs;
that archive is distinct from a running minimum restricted to one persistence
window.

The six-state mechanical study is an **instance-level audit**. Its candidate
bound has a distribution-free interpretation conditional on almost-sure support
reconstruction. The finite numerical checks retain that qualification. Its
50,000-plant validation counts **contraction/LMI violations**. A violation of this
sufficient contraction condition is distinct from a direct instability test.

Performance-selection sensitivity retains both evaluators described in the paper. The alternative
evaluator was introduced after inspecting the aligned result; its comparison
remains exploratory. Exact evaluation risks are used to assess experiments,
while the stopping procedures use their stated certificates and targets.

## Repository layout and file integrity

| Location | Contents |
| --- | --- |
| `code/` | Scientific scripts for the experiments and numerical checks |
| `reference/` | Reference numerical outputs and figure assets |
| `tools/` | Smoke check, table exporter and reference-comparison tools |
| `docs/` | Protocols, experiment map and publication instructions |
| `tests/` | Tests of the reproduction wrapper and its numerical helpers |
| `validation/` | Evidence from testing this prepared repository |
| `.github/workflows/` | Automated checks and a manually dispatched full reproduction |

[SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json) identifies the paper by title
and author and lists the scientific code and reference files with their hashes.
Both checksum manifests are included in the repository root and must be uploaded
with the rest of the files. `certificates_original.py` supports an auxiliary
diagnostic; `certificates.py` implements the certificate used for the main results.

## Citation and release

Please cite **A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design**,
by **G. Calafiore**, when using these results. [CITATION.cff](CITATION.cff) contains
the paper citation and repository address.

The release retains the MIT license text in [LICENSE](LICENSE); the repository
owner should confirm the licensing terms before publication.
[Publication instructions](docs/PUBLISHING.md) explain how to update the existing
repository and check the uploaded files.
