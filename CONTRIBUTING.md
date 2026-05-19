# Contributing to ZebraCHOP

Thanks for your interest! Bug reports, feature ideas, and PRs welcome.

## A heads-up about Python versions

- The core CHOPCHOP engine (`chopchop.py`, `chopchop_query.py`, `featurization.py`) runs on **Python 2.7** — it's a fork of upstream CHOPCHOP and depends on a frozen scikit-learn pickle.
- The web frontend (`frontend/`) and the setup wizard (`bin/zebrachop-setup`) are **Python 3.9+** and shell out to the Py2 engine as a subprocess.

This means most user-facing work goes into the Python 3 codebase; touching the Py2 core should be rare and additive.

## Setup

```bash
git clone https://github.com/abachu2005/ZebraCHOP.git
cd ZebraCHOP
python3 bin/zebrachop-setup           # detects bowtie/twoBitToFa/primer3, writes config_local.json
```

For development:

```bash
python3 -m venv frontend/.venv
./frontend/.venv/bin/pip install -e ".[dev]"
./frontend/.venv/bin/pre-commit install
```

## Tests

```bash
./frontend/.venv/bin/pytest -q
```

## Pull requests

1. Fork and branch off `main`.
2. Update `CHANGELOG.md` and tests.
3. `ruff check .` and `pytest` must pass.

## Code of Conduct

Participation governed by [`CODE_OF_CONDUCT.md`](CODE_OF_CONDUCT.md).
