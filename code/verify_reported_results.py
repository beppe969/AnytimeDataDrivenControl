"""Check manuscript-relevant values against the regenerated per-stream records.

Run from any directory. Uses files within this replication package and writes
checks/reported_results_verification.json. Assertions check printed rounding,
counts, and the qualification attached to the mechanical example.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]


def cp(count: int, total: int) -> list[float]:
    return [float(beta.ppf(.025, count, total-count+1)) if count else 0.0,
            float(beta.ppf(.975, count+1, total-count)) if count < total else 1.0]


def main() -> None:
    checks: dict = {}
    endpoint = pd.read_csv(ROOT / 'data/horizon_endpoint_comparison.csv')
    rows = endpoint[endpoint.N == 100000]
    expected = {
        2: '4.74 8.85 23.25 21.38 0.920',
        5: '9.15 14.08 30.16 29.78 0.987',
        23: '31.41 39.81 61.91 68.63 1.109',
        100: '116.99 133.73 168.67 187.03 1.109',
        1000: '1052.29 1109.39 1196.10 1267.44 1.060',
    }
    for r in rows.itertuples():
        formatted = (f'{r.beta:.2f} {r.CG23:.2f} {r.finite_horizon:.2f} '
                     f'{r.anytime:.2f} {r.any_union:.3f}')
        assert formatted == expected[r.k], (r.k, formatted)
    assert len(rows) == len(expected)
    checks['table_2_all_printed_rows_match'] = True

    original = pd.read_csv(ROOT / 'data/stopping_audit_replicates.csv')
    horizon = pd.read_csv(ROOT / 'data/horizon_stopping_replicates.csv')
    specs = [('Fixed N', original, 93.0, 372),
             ('CG23-repeated', original, 165.5, 10),
             ('GCC refined', original, 156.2, 21),
             ('Finite horizon', horizon, 335.0, 0),
             ('Summable CG23', horizon, 401.9, 0),
             ('Anytime', horizon, 344.1, 0)]
    table3 = []
    for name, df, printed_mean, count in specs:
        d = df[(df.method == name) & (df.eps == .05)]
        assert len(d) == 10000 and int(d.deployed.sum()) == 10000
        assert int(d.risk_exceeds_eps.sum()) == count
        assert f'{d.acquired.mean():.1f}' == f'{printed_mean:.1f}'
        table3.append(dict(method=name, mean_n=float(d.acquired.mean()),
                           errors=count, CP95=cp(count, len(d))))
    second = original[(original.method == 'CG23-repeated') & (original.eps == .15)]
    assert int(second.risk_exceeds_eps.sum()) == 7
    checks['table_3_recomputed_from_raw_records'] = table3

    monitoring = pd.read_csv(ROOT / 'data/horizon_monitoring_8000.csv')
    assert len(monitoring) == 10000
    coverage = []
    for n, count in [(800, 488), (8000, 680)]:
        v = monitoring['CG23']
        actual = int(((v > 0) & (v <= n)).sum())
        assert actual == count
        coverage.append(dict(through=n, crossings=actual, CP95=cp(actual, len(v))))
        for name in ['Finite horizon', 'Summable CG23', 'Anytime']:
            v = monitoring[name]
            assert int(((v > 0) & (v <= n)).sum()) == 0
    checks['table_4_recomputed_from_raw_records'] = coverage
    assert coverage[-1]['CP95'][0] > .05

    archive = pd.read_csv(ROOT / 'data/archive_audit_replicates.csv')
    a = archive[archive.eps == .05]
    assert len(a) == 10000
    counts = dict(affected_streams=int((a.recrossing_intervals > 0).sum()),
                  gaps=int(a.recrossing_intervals.sum()),
                  prefixes=int(a.unavailable_prefixes.sum()),
                  longest_gap=int(a.longest_interval.max()))
    assert counts == dict(affected_streams=498, gaps=499, prefixes=23641, longest_gap=90)
    checks['archive'] = counts

    aligned = pd.read_csv(ROOT / 'data/selection_sensitivity.csv')
    ar = 100 * (1-aligned.aligned_selected_cost.mean()/aligned.aligned_first_cost.mean())
    ai = int((aligned.aligned_selected_cost < aligned.aligned_first_cost-1e-10).sum())
    costs = pd.read_csv(ROOT / 'data/lmi_replicates.csv')
    first = costs[costs.method == 'Anytime'].sort_values('replicate').nominal_cost.to_numpy()
    selected = costs[costs.method == 'Posthoc'].sort_values('replicate').nominal_cost.to_numpy()
    br = 100 * (1-selected.mean()/first.mean())
    bi = int((selected < first-1e-10).sum())
    assert len(aligned) == len(first) == len(selected) == 150
    assert f'{ar:.6f}' == '0.000897' and ai == 1
    assert f'{br:.2f}' == '3.38' and bi == 70
    checks['performance_sensitivity'] = dict(aligned_percent=ar, aligned_improved=ai,
        alternative_percent=br, alternative_improved=bi, alternative_was_post_hoc=True)

    results = json.loads((ROOT / 'data/results.json').read_text())
    disturbance = results['disturbance']
    hull = disturbance['selected']['Hull']; box = disturbance['selected']['Box']
    assert hull['n'] == 3267 and hull['k'] == 23 and box['n'] == 1159
    assert disturbance['boundary_changes'] == 116 and disturbance['maximum_boundary'] == 24
    assert f"{hull['risk']:.6f}" == '0.007243' and f"{hull['bound']:.6f}" == '0.019993'
    assert f"{box['risk']:.6f}" == '0.002584'
    assert f"{hull['mpc']['cost']:.5f}" == '2.87185' and f"{box['mpc']['cost']:.5f}" == '2.95335'
    assert hull['mpc']['feasible'] and box['mpc']['feasible']
    fh = json.loads((ROOT / 'checks/horizon_mpc_verification.json').read_text())
    assert {r['scheme']: r['first_qualified'] for r in fh} == {'hull':2868, 'box':1233}
    assert all(r['mpc']['feasible'] and r['mpc']['terminal_residual'] < 1e-10 for r in fh)
    checks['tube_MPC_first_qualification_and_costs_match'] = True

    mechanical = json.loads((ROOT / 'data/mechanical_results.json').read_text())
    assert mechanical['certificate_status'] == 'conditional_on_almost_sure_support_reconstruction'
    assert mechanical['selected_n'] == 322 and mechanical['selected_k'] == 1
    assert mechanical['boundary_updates'] == 9 and mechanical['max_k'] == 2
    assert f"{mechanical['certificate']:.5f}" == '0.04729'
    validation = mechanical['evaluation']
    assert validation['event'] == 'contraction_LMI_violation'
    assert validation['validation_seed'] == 26090762
    assert validation['N'] == 50000 and validation['failures'] == 34
    assert np.max(np.abs(np.array(validation['CP95'])-cp(34, 50000))) < 1e-15
    assert mechanical['worst_reconstruction_error'] < 4.3e-8
    assert mechanical['worst_KKT_residual'] < 1.2e-7
    checks['mechanical_conditional_bound_and_validation_match'] = True

    checks['all_checks_passed'] = True
    target = ROOT / 'checks/reported_results_verification.json'
    target.write_text(json.dumps(checks, indent=2)+'\n')
    print(json.dumps(checks, indent=2))

if __name__ == '__main__':
    main()
