---
name: Reproduction discrepancy
about: Report a failed command or a numerical difference from the reference records
title: "[reproduction] "
labels: ''
assignees: ''
---

## Command and release

Paste the exact command and repository release/commit.

## Expected and observed result

Identify the paper table, figure or output file. Include relevant numerical values.

## Environment

Attach the output of `python reproduce.py doctor` or `run_manifest.json`.

## Logs

Attach the failed step's log and, when available,
`checks/reference_comparison.json`. Remove any unrelated private filesystem paths
before posting. Keep the original `reference/` files unchanged.
