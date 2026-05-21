"""ZebraCHOP web frontend (Python 3 → wraps Python 2.7 CLI).

Run:

    cd ZebraCHOP
    python3 -m uvicorn frontend.main:app --host 127.0.0.1 --port 8000

The frontend launches ``python2 chopchop_query.py`` as a subprocess. Path to
the Python 2.7 interpreter is read from ``config_local.json["PYTHON2"]`` or
``$PYTHON2``, falling back to ``python2``.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import threading
import uuid
from pathlib import Path

from fastapi import FastAPI, Form, HTTPException
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

ROOT = Path(__file__).resolve().parents[1]
FRONT = Path(__file__).resolve().parent
DATA = FRONT / "data" / "jobs"
DATA.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="ZebraCHOP", version="1.0.0")
if (FRONT / "static").exists():
    app.mount("/static", StaticFiles(directory=str(FRONT / "static")), name="static")

_jobs: dict[str, dict] = {}
_lock = threading.Lock()


def _python2() -> str:
    cfg_file = ROOT / "config_local.json"
    if cfg_file.exists():
        try:
            cfg = json.loads(cfg_file.read_text())
            if cfg.get("PYTHON2"):
                return cfg["PYTHON2"]
        except Exception:
            pass
    return os.environ.get("PYTHON2", "python2")


def _set(job_id: str, **fields):
    with _lock:
        _jobs.setdefault(job_id, {}).update(fields)
        (DATA / job_id / "status.json").write_text(json.dumps(_jobs[job_id], default=str))


def _job_dir(job_id: str) -> Path:
    return DATA / job_id


def _run(job_id: str, genes: list[str], genome: str, mode: int, extra: list[str]):
    job_dir = _job_dir(job_id)
    try:
        _set(job_id, state="running", log="")
        cmd = [_python2(), "chopchop_query.py",
               "--gene_names", ",".join(genes),
               "-o", str(job_dir),
               "--",
               "-G", genome,
               "-T", str(mode)]
        cmd += extra
        log_file = job_dir / "run.log"
        with open(log_file, "w") as lf:
            lf.write("CMD: " + " ".join(cmd) + "\n\n")
            lf.flush()
            r = subprocess.run(cmd, cwd=str(ROOT), stdout=lf, stderr=subprocess.STDOUT)
        if r.returncode == 0:
            _set(job_id, state="done", returncode=0)
        else:
            _set(job_id, state="error", returncode=r.returncode,
                 error=log_file.read_text()[-2000:])
    except FileNotFoundError as e:
        _set(job_id, state="error", error=f"Could not launch python2: {e}. "
             "Set PYTHON2 in config_local.json or your environment.")
    except Exception as e:
        _set(job_id, state="error", error=str(e))


@app.get("/", response_class=HTMLResponse)
def index():
    idx = FRONT / "index.html"
    return HTMLResponse(idx.read_text() if idx.exists() else "<h1>frontend missing</h1>", status_code=200 if idx.exists() else 500)


@app.get("/api/health")
def health():
    return {"ok": True, "version": app.version, "python2": _python2()}


@app.get("/api/config")
def config():
    cfg_file = ROOT / "config_local.json"
    if cfg_file.exists():
        try:
            return JSONResponse(json.loads(cfg_file.read_text()))
        except Exception:
            return JSONResponse({"_error": "config_local.json present but unparseable"}, status_code=500)
    cfg = ROOT / "config.json"
    return JSONResponse(json.loads(cfg.read_text()) if cfg.exists() else {})


@app.post("/api/jobs")
def create_job(
    gene_names: str = Form(...),
    genome: str = Form("danRer11"),
    mode: int = Form(1),
    scoring_method: str = Form("DOENCH_2016"),
    guide_size: int = Form(20),
):
    genes = [g.strip() for g in gene_names.replace("\n", ",").split(",") if g.strip()]
    if not genes:
        raise HTTPException(400, "Provide at least one gene name")
    job_id = uuid.uuid4().hex[:12]
    _job_dir(job_id).mkdir(parents=True, exist_ok=True)
    extra = ["--scoringMethod", scoring_method, "--guideSize", str(guide_size)]
    _set(job_id, id=job_id, state="queued", genes=genes, genome=genome, mode=mode)
    threading.Thread(target=_run, args=(job_id, genes, genome, mode, extra), daemon=True).start()
    return {"job_id": job_id, "genes": genes}


@app.get("/api/jobs/{job_id}")
def status(job_id: str):
    with _lock:
        if job_id in _jobs:
            return _jobs[job_id]
        sf = _job_dir(job_id) / "status.json"
        if sf.exists():
            return json.loads(sf.read_text())
    raise HTTPException(404, "Unknown job")


@app.get("/api/jobs/{job_id}/log")
def job_log(job_id: str):
    lf = _job_dir(job_id) / "run.log"
    if not lf.exists():
        raise HTTPException(404, "no log yet")
    return FileResponse(str(lf), media_type="text/plain")


@app.get("/api/jobs/{job_id}/results")
def job_results(job_id: str):
    """Return parsed TSV results for every gene in the job."""
    job_dir = _job_dir(job_id)
    if not job_dir.exists():
        raise HTTPException(404)
    out = {}
    for tsv in sorted(job_dir.glob("*.tsv")):
        rows = []
        with open(tsv) as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader, None) or []
            for row in reader:
                rows.append(dict(zip(header, row)))
                if len(rows) >= 50:
                    break
        out[tsv.stem] = {"header": header, "rows": rows}
    return out


@app.get("/api/jobs/{job_id}/tsv/{gene}")
def job_tsv(job_id: str, gene: str):
    f = _job_dir(job_id) / f"{gene}.tsv"
    if not f.exists():
        raise HTTPException(404)
    return FileResponse(str(f), media_type="text/tab-separated-values", filename=f.name)


def main() -> None:
    """Console entry point: launch the ZebraCHOP web frontend via uvicorn."""
    import os
    import uvicorn

    host = os.environ.get("ZEBRACHOP_HOST", "127.0.0.1")
    port = int(os.environ.get("ZEBRACHOP_PORT", "8000"))
    uvicorn.run("frontend.main:app", host=host, port=port, reload=False)


if __name__ == "__main__":
    main()
