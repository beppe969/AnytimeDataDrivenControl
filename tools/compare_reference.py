#!/usr/bin/env python3
"""Compare regenerated scientific records with immutable reference outputs.

Integer and text fields must match exactly. Floating-point fields use the
explicit tolerances below. Only elapsed_seconds and recorded software versions
are excluded from JSON comparisons; their values remain in the run artifacts.
PDF files are checked for presence; PNG pixel arrays provide figure comparison.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
RTOL, ATOL = 1e-9, 1e-11
IGNORED_KEYS = {"elapsed_seconds", "software"}
PAPER_FIGURES = ("horizon_comparison", "certificate_archive", "disturbance_risk")
DISCRETE_COLUMNS = {"n", "N", "k", "replicate", "k_hull", "k_box", "terminal_k",
                    "deployed", "acquired", "risk_exceeds_eps", "first",
                    "recrossing_intervals", "unavailable_prefixes", "longest_interval",
                    "first_CG23_failure", "first_anytime_failure"}


def tolerances(path: Path) -> tuple[float, float]:
    # Convex-solver diagnostics are more sensitive to BLAS/LAPACK differences.
    return (1e-6, 1e-7) if "mechanical" in path.name else (RTOL, ATOL)


def compare_json(left, right, rtol: float, atol: float, where: str = "$") -> list[str]:
    differences = []
    if isinstance(left, dict) and isinstance(right, dict):
        lk, rk = set(left)-IGNORED_KEYS, set(right)-IGNORED_KEYS
        if lk != rk:
            differences.append(f"{where}: different keys")
        for key in sorted(lk & rk):
            differences.extend(compare_json(left[key], right[key], rtol, atol, f"{where}.{key}"))
    elif isinstance(left, list) and isinstance(right, list):
        if len(left) != len(right):
            differences.append(f"{where}: different list lengths")
        for index, (a,b) in enumerate(zip(left,right)):
            differences.extend(compare_json(a,b,rtol,atol,f"{where}[{index}]"))
    elif isinstance(left, bool) or isinstance(left, int) or left is None or isinstance(left, str):
        if type(left) is not type(right) or left != right:
            differences.append(f"{where}: {left!r} != {right!r}")
    elif isinstance(left, (int,float)) and isinstance(right, (int,float)):
        if not ((math.isnan(left) and math.isnan(right)) or math.isclose(left,right,rel_tol=rtol,abs_tol=atol)):
            differences.append(f"{where}: {left!r} != {right!r}")
    elif left != right:
        differences.append(f"{where}: incompatible values")
    return differences


def compare_csv(a: Path, b: Path) -> dict:
    x, y = pd.read_csv(a), pd.read_csv(b)
    rtol,atol = tolerances(a)
    issues = []
    maxdiff = 0.0
    if list(x.columns) != list(y.columns) or x.shape != y.shape:
        issues.append("Different columns or shape")
    else:
        for name in x.columns:
            u,v = x[name],y[name]
            if pd.api.types.is_numeric_dtype(u) and pd.api.types.is_numeric_dtype(v):
                integer = name in DISCRETE_COLUMNS or (pd.api.types.is_integer_dtype(u) and pd.api.types.is_integer_dtype(v))
                equal = np.array_equal(u.to_numpy(),v.to_numpy(),equal_nan=True) if integer else np.allclose(u,v,rtol=rtol,atol=atol,equal_nan=True)
                if not equal:issues.append(f"Column {name} differs")
                diff = np.abs(u.to_numpy(dtype=float)-v.to_numpy(dtype=float))
                finite = diff[np.isfinite(diff)]
                if finite.size:maxdiff=max(maxdiff,float(finite.max()))
            elif not u.fillna("<NA>").equals(v.fillna("<NA>")):
                issues.append(f"Column {name} differs")
    return dict(passed=not issues, rows=len(x), maximum_absolute_numeric_difference=maxdiff,
                identical_bytes=a.read_bytes()==b.read_bytes(), rtol=rtol,atol=atol,issues=issues)


def compare_workspace(workspace: Path) -> dict:
    records = []
    for folder in ("data", "checks", "tables"):
        for reference in sorted((ROOT / "reference" / folder).iterdir()):
            if not reference.is_file():continue
            relative = f"{folder}/{reference.name}"
            generated = workspace / relative
            item = {"file": relative}
            if not generated.is_file():
                item.update(passed=False,issues=["Missing regenerated file"])
            elif reference.suffix == ".csv":
                item.update(compare_csv(reference, generated))
            elif reference.suffix == ".json":
                rtol,atol = tolerances(reference)
                differences = compare_json(json.loads(reference.read_text()),json.loads(generated.read_text()),rtol,atol)
                item.update(passed=not differences,rtol=rtol,atol=atol,issues=differences[:30])
            elif reference.suffix == ".npz":
                issues=[]
                with np.load(reference,allow_pickle=False) as x, np.load(generated,allow_pickle=False) as y:
                    if set(x.files)!=set(y.files):issues.append("Different array names")
                    for key in sorted(set(x.files)&set(y.files)):
                        if x[key].shape!=y[key].shape or not np.array_equal(x[key],y[key],equal_nan=True):
                            issues.append(f"Array {key} differs")
                    item['arrays']=x.files
                item.update(passed=not issues,issues=issues)
            else:
                same = reference.read_bytes()==generated.read_bytes()
                item.update(passed=same,identical_bytes=same)
            records.append(item)
    figures=[]
    for name in PAPER_FIGURES:
        a=ROOT/'reference/figures'/f'{name}.png';b=workspace/'figures'/f'{name}.png'
        pdf=workspace/'figures'/f'{name}.pdf'
        entry={'figure':name,'pdf_present':pdf.is_file(),'png_present':b.is_file()}
        if b.is_file():
            with Image.open(a) as ia, Image.open(b) as ib:
                x,y=np.asarray(ia.convert('RGBA')),np.asarray(ib.convert('RGBA'))
                same=x.shape==y.shape and np.array_equal(x,y)
                entry.update(pixel_identical=same,reference_shape=list(x.shape),generated_shape=list(y.shape))
                if x.shape==y.shape:entry['maximum_channel_difference']=int(np.abs(x.astype(int)-y.astype(int)).max())
        entry['passed']=entry['pdf_present'] and entry['png_present']
        figures.append(entry)
    return {'all_numerical_records_match':all(x['passed'] for x in records),
            'all_paper_figures_present':all(x['passed'] for x in figures),
            'all_paper_pngs_pixel_identical':all(x.get('pixel_identical',False) for x in figures),
            'figure_policy':'PNG pixel identity is reported separately; font/rendering differences do not fail numerical equivalence. PDF presence is checked; PDF bytes may contain creation metadata.',
            'ignored_json_keys':sorted(IGNORED_KEYS),'records':records,'figures':figures}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace',type=Path,required=True)
    args=parser.parse_args()
    workspace=args.workspace.resolve()
    if not workspace.is_dir():parser.error('Workspace does not exist')
    result=compare_workspace(workspace)
    (workspace/'checks').mkdir(exist_ok=True)
    (workspace/'checks/reference_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return 0 if result['all_numerical_records_match'] and result['all_paper_figures_present'] else 1

if __name__=='__main__':raise SystemExit(main())
