# Publishing and updating the GitHub repository

This package generates the numerical results for the paper
**A Universal Iterated-Logarithm Anytime Law for Data-Driven Control Design**,
by **G. Calafiore**.

## Update the existing repository

Copy the contents of this package into the local checkout of
`beppe969/AnytimeDataDrivenControl`, preserving its existing `.git/` directory.
Include the hidden `.github/` directory and both hidden configuration files
`.gitattributes` and `.gitignore`. Replace the package's documentation and
metadata files as well as its execution tools.

Two required files must be present at the repository root:

```text
SOURCE_MANIFEST.sha256
SHA256SUMS
```

`SOURCE_MANIFEST.sha256` is read before a numerical run. `SHA256SUMS` checks the
complete release. Both files must be committed alongside `SOURCE_PROVENANCE.json`
and `reproduce.py`. The release builder fails when a required file is missing.

After copying the files, run from the repository root:

```bash
python reproduce.py integrity --release
python -m unittest discover -s tests -v
python reproduce.py verify
```

To regenerate all experiments from empty output directories:

```bash
python reproduce.py full --run-dir runs/paper
```

Commit and push the updated files using the repository's existing remote.
Generated workspaces under `runs/` and local Python caches are excluded by
`.gitignore`. Verify that both checksum manifests are visible in the root of
the public repository after the update.

## Metadata and licensing

`CITATION.cff` contains the paper title and author, together with the repository
address. Update the bibliographic metadata when a publication identifier becomes
available. The MIT license text is retained in `LICENSE`; the repository owner
should confirm the licensing terms before publication.

## Continuous integration

The ordinary workflow runs the unit tests, numerical verification and figure
generation. The full numerical reproduction workflow is manually dispatched.
Both configurations preserve generated reports as workflow artifacts.
Local validation evidence and its limits are recorded in `validation/`.

## Prepare a versioned code release

After documentation or metadata changes, regenerate the release-wide checksums:

```bash
python tools/build_release.py --manifest-only
python reproduce.py integrity --release
```

Create a distributable code package with:

```bash
python tools/build_release.py
```

The code package and its SHA-256 sidecar are placed under `dist/`. The builder
includes only the declared repository files. New run workspaces and local
configuration remain outside the distributable package.

The release version is `1.0.1`. Keep this value synchronized in `CITATION.cff`
and `tools/build_release.py` when preparing a later release.
