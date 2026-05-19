# Changelog

All notable changes to this project are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] - 2026-05-19

First public release suitable for general OSS distribution.

### Added
- Python-3 FastAPI **web frontend** (`frontend/`) that subprocesses the
  Python-2.7 `chopchop_query.py` CLI; modern UI with per-gene result tabs,
  sortable guide table, TSV download, live run log.
- Interactive **setup wizard** at `bin/zebrachop-setup` that detects
  `bowtie`/`twoBitToFa`/`primer3`, optionally downloads the danRer11 `.2bit`
  reference, and writes `config_local.json`.
- `README.md`, `LICENSE` (MIT, with upstream attribution), `CITATION.cff`,
  `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `CHANGELOG.md`.
- `requirements.txt` (Py2 pins for the engine) and `frontend/requirements.txt` (Py3).
- `.gitignore` for venvs, `config_local.json`, large reference indexes, runtime data.
- GitHub Actions CI: lint + smoke tests on Python 3.9–3.12.
- Dockerfile that bundles Python 2.7 + Python 3 + bowtie + UCSC kentutils so
  users don't have to install the toolchain manually.

### Removed
- `featurization.pyc` (compiled bytecode should never be committed).

[Unreleased]: https://github.com/abachu2005/ZebraCHOP/compare/v1.0.0...HEAD
[1.0.0]: https://github.com/abachu2005/ZebraCHOP/releases/tag/v1.0.0
