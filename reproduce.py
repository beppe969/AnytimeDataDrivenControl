#!/usr/bin/env python3
"""Run the paper's numerical experiments in isolated, auditable workspaces.

The scientific scripts in code/ are copied unchanged into each workspace. A full
run begins with empty output directories. Only verify/figures and selected
experiments seed explicitly documented inputs from reference/.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                    "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")
VERIFY = ["verify", "verify_revision", "final_checks", "third_review_checks", "verify_reported_results"]
FULL = ["certificates", "verify", "verify_revision", "experiments", "revision_experiments",
        "sdp_scale", "selection_audit", "monitoring_audit", "final_figures", "final_checks",
        "horizon_baselines", "third_review_figures", "third_review_checks", "verify_reported_results"]
EXPERIMENTS = {
    "controller-mpc": ["experiments", "revision_experiments"],
    "stopping-archive": ["selection_audit", "final_figures"],
    "monitoring": ["monitoring_audit", "final_figures"],
    "horizon-comparison": ["horizon_baselines", "third_review_figures", "third_review_checks"],
    "mechanical": ["sdp_scale"],
}


def sha256(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def check_manifest(manifest: Path) -> dict:
    failures = []
    count = 0
    if not manifest.is_file():
        return {"manifest": manifest.name, "checked_files": 0, "passed": False,
                "failures": [{"path": manifest.name, "reason": "missing required manifest; restore it from the repository release"}]}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        expected, relative = line.split("  ", 1)
        path = (ROOT / relative).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError(f"Unsafe path in manifest: {relative}")
        count += 1
        if not path.is_file():
            failures.append({"path": relative, "reason": "missing"})
        elif sha256(path) != expected:
            failures.append({"path": relative, "reason": "SHA-256 mismatch"})
    return {"manifest": manifest.name, "checked_files": count,
            "passed": not failures, "failures": failures}


def check_dependencies() -> dict:
    packages = {}
    issues = []
    for line in (ROOT / "requirements-lock.txt").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        name, expected = line.split("==", 1)
        try:
            observed = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            observed = None
        packages[name] = {"required": expected, "installed": observed}
        if observed != expected:
            issues.append(f"{name}: expected {expected}; found {observed}")
    if sys.version_info[:3] != (3, 13, 5):
        issues.append(f"Python {platform.python_version()}; reference interpreter is 3.13.5")
    return {"python": platform.python_version(), "implementation": platform.python_implementation(),
            "platform": platform.platform(), "machine": platform.machine(),
            "packages": packages, "matches_reference_environment": not issues, "issues": issues}


def prepare_workspace(path: Path, seed_reference: bool) -> dict:
    path = path.resolve()
    # Never permit an invocation to turn the repository or an ancestor into a run.
    if path == ROOT or ROOT.is_relative_to(path):
        raise ValueError("The workspace must differ from the repository and its ancestors.")
    for folder in ("code", "reference", "docs", "tools", "tests", "validation", ".github"):
        if path.is_relative_to(ROOT / folder):
            raise ValueError(f"Choose a workspace outside the protected {folder}/ directory.")
    if path.exists() and any(path.iterdir()):
        raise ValueError(f"Workspace is not empty: {path}. Choose a new --run-dir; existing runs are never overwritten.")
    path.mkdir(parents=True, exist_ok=True)
    shutil.copytree(ROOT / "code", path / "code", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
    for folder in ("data", "figures", "tables", "checks", "logs"):
        (path / folder).mkdir()
    seeded = []
    if seed_reference:
        # No old verification reports or figures are copied into a new run.
        for folder in ("data", "checks"):
            for source in sorted((ROOT / "reference" / folder).iterdir()):
                if source.is_file():
                    target = path / folder / source.name
                    shutil.copyfile(source, target)
                    seeded.append({"path": f"{folder}/{source.name}", "sha256": sha256(source)})
    return {"seeded_reference_inputs": seeded, "starts_from_empty_outputs": not seed_reference}


def child_environment(workspace: Path) -> dict[str, str]:
    env = os.environ.copy()
    env.pop("PYTHONOPTIMIZE", None)  # Source verification scripts use assertions.
    env.pop("PYTHONPATH", None)
    for name in THREAD_VARIABLES:
        env[name] = "1"
    env.update(MPLBACKEND="Agg", MPLCONFIGDIR=str(workspace / ".mplconfig"),
               NUMBA_CACHE_DIR=str(workspace / ".numba_cache"), PYTHONHASHSEED="0",
               PYTHONUTF8="1", PYTHONUNBUFFERED="1")
    return env


def run_step(label: str, argv: list[str], workspace: Path, env: dict[str, str]) -> dict:
    relative_log = f"logs/{label}.log"
    print(f"Running {label}; log: {workspace / relative_log}", flush=True)
    start = time.perf_counter()
    with (workspace / relative_log).open("w", encoding="utf-8") as stream:
        result = subprocess.run(argv, cwd=workspace, env=env, stdout=stream,
                                stderr=subprocess.STDOUT, check=False)
    elapsed = time.perf_counter() - start
    def portable_argument(value: str) -> str:
        if value == sys.executable:
            return "python"
        path = Path(value)
        if path.is_absolute() and path.is_relative_to(workspace):
            return str(path.relative_to(workspace))
        if path.is_absolute() and path.is_relative_to(ROOT):
            return "$REPOSITORY/" + str(path.relative_to(ROOT))
        return value
    item = {"step": label, "command": [portable_argument(x) for x in argv],
            "returncode": result.returncode, "elapsed_seconds": elapsed, "log": relative_log}
    print(f"  {label}: {'PASS' if result.returncode == 0 else 'FAIL'} ({elapsed:.2f} s)", flush=True)
    if result.returncode:
        print("\n".join((workspace / relative_log).read_text(errors="replace").splitlines()[-30:]), file=sys.stderr)
    return item


def run_pipeline(args: argparse.Namespace) -> int:
    if sys.flags.optimize:
        raise ValueError("Use ordinary Python without -O/-OO: numerical checks require assertions.")
    source_audit = check_manifest(ROOT / "SOURCE_MANIFEST.sha256")
    if not source_audit["passed"]:
        raise ValueError("Source/reference files differ from their manifest. Inspect with `python reproduce.py integrity`.")
    environment = check_dependencies()
    if environment["issues"]:
        print("Environment differences:\n  " + "\n  ".join(environment["issues"]), file=sys.stderr)
        if not args.allow_environment_drift:
            raise ValueError("Install requirements-lock.txt with Python 3.13.5, or explicitly use --allow-environment-drift.")
    mode = args.command
    seeded = mode in ("verify", "figures", "experiment")
    tag = args.name if mode == "experiment" else mode
    directory = Path(args.run_dir) if args.run_dir else ROOT / "runs" / (tag + "-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    workspace = directory.expanduser().resolve()
    preparation = prepare_workspace(workspace, seeded)
    env = child_environment(workspace)
    provenance = json.loads((ROOT / "SOURCE_PROVENANCE.json").read_text())
    record = {"schema_version": 2, "mode": mode, "experiment": getattr(args, "name", None),
              "started_utc": datetime.now(timezone.utc).isoformat(), "status": "running",
              "paper": provenance["paper"],
              "package_version": provenance["package_version"],
              "source_manifest_sha256": sha256(ROOT / "SOURCE_MANIFEST.sha256"),
              "orchestration_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in
                  [ROOT / "reproduce.py", *sorted((ROOT / "tools").glob("*.py"))]},
              "source_integrity": source_audit, "environment": environment,
              "execution_environment": {k: env[k] for k in THREAD_VARIABLES},
              "certificate_variant": "original paper certificate; no persistence-window-minimum modification",
              **preparation, "steps": []}
    if mode == "full":
        scripts = FULL
    elif mode == "verify":
        scripts = VERIFY
    elif mode == "figures":
        scripts = ["final_figures", "third_review_figures", "verify_reported_results"]
    elif mode == "smoke":
        scripts = []
    else:
        scripts = EXPERIMENTS[args.name]
    record["planned_scientific_scripts"] = scripts
    manifest_path = workspace / "run_manifest.json"
    write_json(manifest_path, record)
    start = time.perf_counter()
    commands = [(script, [sys.executable, str(workspace / "code" / (script + ".py"))]) for script in scripts]
    if mode == "smoke":
        commands = [("smoke", [sys.executable, str(ROOT / "tools/smoke.py"), "--workspace", str(workspace)])]
    if mode in ("full", "verify", "figures"):
        commands.append(("export_paper_tables", [sys.executable, str(ROOT / "tools/export_paper_tables.py"), "--workspace", str(workspace)]))
    if mode == "full":
        commands.append(("compare_reference", [sys.executable, str(ROOT / "tools/compare_reference.py"), "--workspace", str(workspace)]))
    try:
        for label, command in commands:
            item = run_step(label, command, workspace, env)
            record["steps"].append(item)
            write_json(manifest_path, record)
            if item["returncode"]:
                record["status"] = "failed"
                break
        else:
            record["status"] = "passed"
    except KeyboardInterrupt:
        record["status"] = "interrupted"
        raise
    except Exception as exc:
        record["status"] = "failed"
        record["error"] = str(exc)
        raise
    finally:
        record["elapsed_seconds"] = time.perf_counter() - start
        record["completed_utc"] = datetime.now(timezone.utc).isoformat()
        record["output_inventory"] = [
            {"path": str(p.relative_to(workspace)), "sha256": sha256(p), "bytes": p.stat().st_size}
            for folder in ("data", "figures", "tables", "checks", "logs")
            for p in sorted((workspace / folder).glob("*")) if p.is_file()]
        write_json(manifest_path, record)
    print(f"{record['status'].upper()}: {manifest_path}")
    return 0 if record["status"] == "passed" else 1


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("doctor", help="Report the interpreter and pinned dependency versions; no calculations.")
    integrity = sub.add_parser("integrity", help="Verify SHA-256 hashes of retained source and reference files.")
    integrity.add_argument("--release", action="store_true", help="Also check SHA256SUMS for the packaged release.")
    helptexts = {"verify": "Run all five numerical checks on copies of the archived results.",
                 "figures": "Regenerate the three paper figures and paper tables from archived data.",
                 "full": "Recompute all numerical experiments from empty outputs and compare with reference/.",
                 "smoke": "Test certificate formulas and a small QP path; does not reproduce paper statistics.",
                 "experiment": "Run one study group with archived prerequisite inputs."}
    for name, text in helptexts.items():
        p = sub.add_parser(name, help=text)
        if name == "experiment":
            p.add_argument("name", choices=sorted(EXPERIMENTS))
        p.add_argument("--run-dir", help="New or empty workspace. Default: timestamped directory inside runs/.")
        p.add_argument("--allow-environment-drift", action="store_true", help="Record and permit dependency/interpreter differences.")
    compare = sub.add_parser("compare", help="Compare a completed full workspace against reference/.")
    compare.add_argument("--run-dir", required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "doctor":
            result = check_dependencies()
            print(json.dumps(result, indent=2))
            return 0 if result["matches_reference_environment"] else 1
        if args.command == "integrity":
            results = [check_manifest(ROOT / "SOURCE_MANIFEST.sha256")]
            if args.release:
                results.append(check_manifest(ROOT / "SHA256SUMS"))
            print(json.dumps(results, indent=2))
            return 0 if all(r["passed"] for r in results) else 1
        if args.command == "compare":
            return subprocess.call([sys.executable, str(ROOT / "tools/compare_reference.py"),
                                    "--workspace", str(Path(args.run_dir).resolve())])
        return run_pipeline(args)
    except (ValueError, FileNotFoundError, NotADirectoryError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
