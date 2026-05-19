"""Frontend smoke tests (no Python 2 engine required)."""
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_frontend_imports():
    import importlib.util
    spec = importlib.util.spec_from_file_location("zc_app", ROOT / "frontend" / "main.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert hasattr(module, "app")
    assert "zebrachop" in module.app.title.lower()


def test_setup_wizard_parses():
    import ast
    src = (ROOT / "bin" / "zebrachop-setup").read_text()
    ast.parse(src)


def test_config_json_template_ok():
    import json
    cfg = json.loads((ROOT / "config.json").read_text())
    assert "PATH" in cfg
    for k in ("BOWTIE", "TWOBITTOFA"):
        assert k in cfg["PATH"]


def test_gene_table_present():
    f = ROOT / "genetable" / "danRer11.gene_table"
    assert f.exists(), "danRer11.gene_table must ship with the repo"
    assert f.stat().st_size > 1_000_000


@pytest.mark.parametrize("path", ["chopchop.py", "chopchop_query.py", "reformat.py", "featurization.py"])
def test_python2_engine_files_present(path):
    assert (ROOT / path).exists()


def test_health_endpoint_responsive():
    import importlib.util
    spec = importlib.util.spec_from_file_location("zc_app", ROOT / "frontend" / "main.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    from fastapi.testclient import TestClient
    client = TestClient(module.app)
    r = client.get("/api/health")
    assert r.status_code == 200
    assert r.json()["ok"] is True
