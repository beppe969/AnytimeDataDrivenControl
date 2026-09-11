# Contributing

Preserve the declared scientific baseline when proposing engineering changes.
Run `python -m unittest discover -s tests -v` and `python reproduce.py verify`.
A change to numerical code additionally requires `python reproduce.py full` and
a description of any numerical differences from `reference/`.

Scientific scripts are retained verbatim for the prepared release. Proposed
refactoring should be separated from changes to the certificate or experimental
protocols. Record seeds and confidence allocations for new experiments.

Report reproduction problems with the issue template. Include the exact command
and its run manifest, together with the failed log. The mechanical example's
conditional interpretation and the exploratory evaluator comparison should remain
explicit in derivative analyses.
