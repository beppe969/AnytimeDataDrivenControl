#!/usr/bin/env python3
"""Build a clean, deterministic release archive from this repository.

Only declared release directories and root files are included. Generated runs,
virtual environments and local metadata never enter the archive.
"""
from __future__ import annotations
import argparse
import hashlib
from pathlib import Path
import stat
import zipfile

ROOT=Path(__file__).resolve().parents[1]
VERSION='1.0.1'
ROOT_FILES={'.gitattributes','.gitignore','.python-version','README.md','LICENSE',
            'CITATION.cff','CHANGELOG.md','CONTRIBUTING.md','Makefile','reproduce.py',
            'requirements.txt','requirements-lock.txt','SOURCE_PROVENANCE.json',
            'SOURCE_MANIFEST.sha256','SHA256SUMS'}
DIRECTORIES={'code','reference','docs','tools','tests','validation','.github'}
EXCLUDED_PARTS={'__pycache__','.pytest_cache','.DS_Store','__MACOSX'}
EXCLUDED_SUFFIXES={'.pyc','.pyo','.nbc','.nbi','.ttf','.otf','.woff','.woff2'}


def release_files():
    required = ROOT_FILES - {'SHA256SUMS'}
    missing = sorted(name for name in required if not (ROOT/name).is_file())
    missing += sorted(folder+'/' for folder in DIRECTORIES if not (ROOT/folder).is_dir())
    if missing:
        raise FileNotFoundError('Required release files or directories are missing: '+', '.join(missing))
    candidates=[ROOT/name for name in ROOT_FILES if (ROOT/name).is_file()]
    for folder in DIRECTORIES:
        candidates.extend(p for p in (ROOT/folder).rglob('*') if p.is_file())
    retained=[]
    for p in candidates:
        relative=p.relative_to(ROOT)
        if p.is_symlink():raise ValueError(f'Symlinks are excluded from releases: {relative}')
        if any(x in EXCLUDED_PARTS or x.startswith('._') for x in relative.parts):continue
        if p.suffix.lower() in EXCLUDED_SUFFIXES:continue
        retained.append(p)
    return sorted(retained,key=lambda p:p.relative_to(ROOT).as_posix())


def digest(path: Path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def write_manifest():
    files=[p for p in release_files() if p.name!='SHA256SUMS']
    (ROOT/'SHA256SUMS').write_text(''.join(f'{digest(p)}  {p.relative_to(ROOT).as_posix()}\n' for p in files),encoding='utf-8')
    return len(files)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-only',action='store_true')
    parser.add_argument('--output',type=Path,help='Archive path; default dist/AnytimeDataDrivenControl-v1.0.1.zip')
    args=parser.parse_args()
    # Validate the scientific manifest before writing any release metadata.
    import sys
    sys.path.insert(0, str(ROOT))
    import reproduce
    audit = reproduce.check_manifest(ROOT/'SOURCE_MANIFEST.sha256')
    if not audit['passed']:
        raise ValueError('Scientific integrity check failed: '+str(audit['failures']))
    count=write_manifest()
    print(f'Wrote SHA256SUMS for {count} release files.')
    if args.manifest_only:return
    output=args.output or ROOT/'dist'/f'AnytimeDataDrivenControl-v{VERSION}.zip'
    output=output.expanduser().resolve();output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError(f'Archive already exists: {output}. Choose another --output.')
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for path in release_files():
            name='AnytimeDataDrivenControl/'+path.relative_to(ROOT).as_posix()
            info=zipfile.ZipInfo(name,date_time=(2026,9,11,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            mode=0o755 if path.suffix=='.py' else 0o644
            info.create_system=3;info.external_attr=(stat.S_IFREG|mode)<<16
            z.writestr(info,path.read_bytes(),compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    checksum=output.with_suffix(output.suffix+'.sha256')
    checksum.write_text(f'{digest(output)}  {output.name}\n')
    print(f'Archive: {output}\nSHA-256: {checksum}')

if __name__=='__main__':main()
