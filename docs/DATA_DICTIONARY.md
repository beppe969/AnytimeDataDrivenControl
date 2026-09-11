# Data dictionary

All original data files live under `reference/`; regenerated versions live under
an individual `runs/.../` workspace. The data are synthetic model observations
or outputs of the stated numerical procedures.

[Machine-readable schemas](data_schema.json) record every CSV column and its
inferred dtype, together with missing-value counts. The schemas preserve the
inherited column names. A column stored as floating point can still represent an
integer-valued sample index or support count.

## CSV inventory

| File relative to `reference/` | Rows | Meaning |
| --- | ---: | --- |
| `data/archive_audit_replicates.csv` | 20,000 | Per-stream archive qualification gaps, for both target values. |
| `data/complexity_sweep.csv` | 13 | Auxiliary comparison of current and original epoch schedules. |
| `data/disturbance_path.csv` | 8,000 | Full hull/box trajectory with exact missing mass and anytime certificates. |
| `data/final_fixed_time_comparison.csv` | 5 | Auxiliary five-complexity fixed-time comparison. |
| `data/fixed_time_comparison.csv` | 12 | Auxiliary fixed-time/schedule comparisons with inherited diagnostic labels. |
| `data/horizon_crossover_sweep.csv` | 724 | Deterministic endpoint grid for Figure 1; no streams of length N are generated. |
| `data/horizon_disturbance_path.csv` | 8,000 | Disturbance trajectory augmented with finite-horizon certificates. |
| `data/horizon_endpoint_comparison.csv` | 14 | Scaled endpoint bounds; the N=100000 rows produce Table 2. |
| `data/horizon_monitoring_800.csv` | 10,000 | First bound-crossing indices in the short first-target cohort. |
| `data/horizon_monitoring_8000.csv` | 10,000 | First bound-crossing indices in the independent long cohort; both Table 4 columns use this file. |
| `data/horizon_stopping_replicates.csv` | 80,000 | FH/summable and replayed CG23/anytime first-target records. |
| `data/lmi_path.csv` | 800 | Scalar-state two-input QP trajectory; this filename refers to the scalar LMI representation. |
| `data/lmi_replicates.csv` | 900 | 150-stream QP/evaluator study, six methods per stream. |
| `data/mechanical_path.csv` | 800 | Six-state mechanical optimization path, including all 18 gain entries. |
| `data/monitoring_example.csv` | 8,000 | Illustrative long-cohort path with current risk and bounds. |
| `data/monitoring_replicates.csv` | 10,000 | Long-cohort replay records for unadjusted CG23 and anytime. |
| `data/recrossing_example.csv` | 800 | First affected short-cohort stream; the plotted example has zero-based index 7. |
| `data/selection_sensitivity.csv` | 150 | Aligned nominal evaluator: first-qualified and retrospectively selected costs. |
| `data/specification_path.csv` | 4,000 | Auxiliary scalar example for post-acquisition specification choice. |
| `data/stopping_audit_replicates.csv` | 80,000 | Fixed-sample, CG23-repeated, GCC and anytime first-target records. |
| `checks/horizon_high_precision_roots.csv` | 60 | Independent mpmath checks of comparator roots. |
| `checks/horizon_root_checks.csv` | 69 | Incomplete-beta roots compared with direct finite summation. |
| `checks/lmi_support_audit.csv` | 1,546 | Scalar QP support-reconstruction checks. |
| `checks/mechanical_audit.csv` | 9 | Mechanical optimizer-update audits of feasibility and reconstruction. |

## Conventions and key fields

| Field or pattern | Meaning |
| --- | --- |
| `n`, `N` | Number of acquired design observations or the prespecified acquisition cap; the file context distinguishes these roles |
| `replicate` | Zero-based Monte Carlo stream index |
| `k`, `k_hull`, `k_box` | Observed boundary/support complexity of the indicated design |
| `eps` | Target violation probability; represented as a fraction |
| `risk`, `risk_hull`, `risk_box` | Exact evaluation risk of the two-dimensional design; represented as a fraction |
| `U`, `U_anytime`, `U_hull`, `U_box` | Original anytime upper risk certificate |
| `U_per_time` | Source diagnostic with summable per-prefix allocation based on the local bound; distinct from the summable CG23 comparator |
| `finite_hull`, `finite_box` | FH certificates with the prescribed cap 8000 |
| `deployed` | Indicator that the procedure returned a design by its simulated horizon |
| `acquired` | Number of observations at stopping; the source convention uses the simulation cap for a censored run |
| `risk_exceeds_eps` | Indicator of a returned design whose exact risk exceeds its target |
| `first` | First anytime-qualified prefix in an archive record |
| `recrossing_intervals` | Number of later contiguous gaps in current-design qualification |
| `unavailable_prefixes` | Total number of prefixes in those qualification gaps |
| `longest_interval` | Longest qualification gap, in numbers of prefixes |
| `first_CG23_failure`, `first_anytime_failure` | First crossing index, or zero if no crossing occurred |
| Comparator columns in `horizon_monitoring_*.csv` | First crossing indices, with zero denoting no crossing; these columns contain indices rather than bound values |
| Comparator columns in `horizon_endpoint_comparison.csv` | Bounds multiplied by the endpoint sample size N |
| `any_union` | Anytime/FH ratio in the horizon-comparison records |
| `L0` through `L17` in `mechanical_path.csv` | Row-major entries of the 3-by-6 feedback gain |
| `objective`, `nominal_cost` | Values of the corresponding optimization/evaluation objective specified in the source protocol |

The manuscript's Table 2 uses **scaled bounds**. The source's auxiliary
`final_fixed_time_comparison.csv` also stores unscaled probabilities in `beta`,
`CG23` and `anytime`, alongside explicitly scaled `n_...` columns. Consult the
specific file schema when combining records.

The two main cohorts use different seeds and sample lengths. To reproduce the
paper, preserve the RNG draw order and the full protocol dimensions. A shorter
smoke path establishes execution only; changing dimensions can alter the sample
sequence produced by a stream-level generator.

## JSON summaries

| File | Contents |
| --- | --- |
| `data/results.json` | Original disturbance/QP studies, including nominal MPC costs and the 150-stream evaluator study |
| `data/revision_results.json` | Auxiliary comparisons, specification choice and evaluator sensitivity |
| `data/final_selection_results.json` | First-target stopping, archive summaries and fast/slow QP cross-checks |
| `data/final_monitoring_results.json` | Independent long-cohort coverage summaries and implementation residuals |
| `data/horizon_baseline_results.json` | Deterministic bounds, FH/summable stopping and continued inspection |
| `data/mechanical_results.json` | Mechanical model/design metadata and independent contraction-violation validation |
| `checks/horizon_mpc_verification.json` | MPC feasibility and cost at FH qualification for hull and box |

The mechanical `certificate_status` field explicitly states
`conditional_on_almost_sure_support_reconstruction`. Its `evaluation.event` is
`contraction_LMI_violation`. The original `CP95` intervals are probabilities,
while the manuscript displays selected intervals as percentages.

## Stored disturbance arrays

`data/disturbance_samples.npz` is readable using `numpy.load(..., allow_pickle=False)`.

| Key | Shape | Meaning |
| --- | --- | --- |
| `W` | 8000 × 2 | Entire sampled disturbance stream |
| `G` | 2 × 2 | Matrix mapping the unit box to the true parallelogram |
| `Box` | 4 × 2 | Box vertices at first anytime qualification |
| `Hull` | 23 × 2 | Hull vertices at first anytime qualification |
| `Hull_per_time` | 22 × 2 | Hull vertices at the auxiliary per-time qualification |

## Generated release tables

The new exporter writes `table2.csv` through `table4.csv` and matching LaTeX
files. Table 3 records the six methods with unrounded means and exact binomial
intervals. Table 4 stores one row per method/inspection horizon. The Markdown
summary applies the manuscript display precision to these recomputed quantities.
