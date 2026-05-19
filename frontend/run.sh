#!/usr/bin/env bash
# Launch the ZebraCHOP web UI. Creates a Python-3 venv at frontend/.venv on
# first run. (The CLI under chopchop.py still runs on Python 2.7.)
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"
PY3=${PY3:-python3}

cd "$ROOT"

if [ ! -d frontend/.venv ]; then
  echo "Creating frontend/.venv …"
  "$PY3" -m venv frontend/.venv
  ./frontend/.venv/bin/pip install --upgrade pip
  ./frontend/.venv/bin/pip install -r frontend/requirements.txt
fi

if [ ! -f config_local.json ]; then
  echo
  echo "  NOTE: config_local.json not found."
  echo "  Run  python3 bin/zebrachop-setup  to configure your bowtie/twoBitToFa paths."
  echo
fi

echo "Open http://127.0.0.1:8000"
exec ./frontend/.venv/bin/python -m uvicorn frontend.main:app --host 127.0.0.1 --port 8000 "$@"
