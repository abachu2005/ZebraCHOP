# ZebraCHOP

> CRISPR guide-RNA design for zebrafish (*Danio rerio*), with Ensembl-aware post-processing and a clean web UI.

ZebraCHOP is a zebrafish-focused fork of [CHOPCHOP](https://chopchop.cbu.uib.no/). It picks high-scoring Cas9/TALEN/Cpf1/Nickase guides for any zebrafish gene, scores them with the published efficiency models (Doench 2016, Xu 2015, …), maps off-targets with Bowtie, and post-processes the output via Ensembl REST to bias selection toward early exons (more likely to produce loss-of-function).

![status: ready](https://img.shields.io/badge/status-ready-brightgreen) ![CLI: Python 2.7](https://img.shields.io/badge/CLI-Python%202.7-orange) ![web: Python 3.9+](https://img.shields.io/badge/web-Python%203.9%2B-blue) ![license: MIT](https://img.shields.io/badge/license-MIT-lightgrey) [![DOI](https://zenodo.org/badge/853115972.svg)](https://zenodo.org/badge/latestdoi/853115972) [![CI](https://github.com/abachu2005/ZebraCHOP/actions/workflows/ci.yml/badge.svg)](https://github.com/abachu2005/ZebraCHOP/actions/workflows/ci.yml)

---

## What it does

```
   Gene name(s)
       │
       ▼
┌────────────────┐   ┌────────────────┐   ┌────────────────┐   ┌────────────────┐
│ Resolve via    │ → │ Extract        │ → │ Score guides   │ → │ Filter via     │
│ danRer11.gene  │   │ sequence with  │   │ + off-targets  │   │ Ensembl exon   │
│ _table         │   │ twoBitToFa     │   │ via Bowtie     │   │ structure      │
└────────────────┘   └────────────────┘   └────────────────┘   └────────────────┘
       │                                                              │
       └──────────────────────────► Per-gene .tsv (ranked) ◄──────────┘
```

Inputs: gene name(s) (or a genePred file). Outputs: one TSV per gene with ranked guides, off-target counts, GC%, self-complementarity, efficiency score, etc., plus an optional Ensembl-aware text summary.

---

## Quick start

```bash
git clone https://github.com/abachu2005/ZebraCHOP.git
cd ZebraCHOP
python3 bin/zebrachop-setup     # interactive: detects bowtie/twoBitToFa/primer3, writes config_local.json
bash frontend/run.sh            # open http://127.0.0.1:8000
```

If you prefer the CLI:

```bash
python2 chopchop_query.py --gene_names rx3,tbx16 -o results/ -- -G danRer11
python  reformat.py -input results/ -output ranked_summary.txt
```

### Why two Python versions?

The core CHOPCHOP engine (`chopchop.py`, `chopchop_query.py`, `featurization.py`) is Python 2.7 — it's a fork of upstream CHOPCHOP and depends on a frozen scikit-learn pickle (`Doench_2016_18.01_model_nopos.pickle`) that doesn't survive a clean Py3 port. The new **web frontend** and **setup wizard** are Python 3 and shell out to the Py2 engine as a subprocess.

The setup wizard records your Python-2 path under `config_local.json["PYTHON2"]` so you only have to set it once.

---

## External dependencies

The setup wizard checks for these and helps you install/configure them:

| Tool | Why |
|---|---|
| `bowtie` | Off-target alignment against the zebrafish genome |
| `twoBitToFa` (UCSC kent utils) | Extract genomic sequence around target |
| `primer3` (optional) | Primer design around the chosen guide |
| `danRer11.2bit` | UCSC two-bit genome file (~770 MB, auto-downloadable) |
| Bowtie index | Built from the zebrafish FASTA via `bowtie-build` |

You'll also want the included `genetable/danRer11.gene_table` (~4.8 MB) which the wizard already points at.

---

## Repository layout

```
.
├── chopchop.py                  # Py2: core guide-design engine
├── chopchop_query.py            # Py2: batch wrapper, one TSV per gene
├── reformat.py                  # Py3 OK: Ensembl REST exon-aware post-processor
├── featurization.py             # Py2: ML feature builder (Doench 2016)
├── config.json                  # default tool paths (empty; copy to config_local.json)
├── frontend/                    # NEW: Python-3 FastAPI web UI (subprocesses chopchop_query.py)
│   ├── main.py
│   ├── index.html
│   ├── requirements.txt
│   └── run.sh
├── bin/zebrachop-setup          # NEW: Python-3 interactive setup wizard
├── genetable/danRer11.gene_table
├── results/                     # example output TSVs (rx3, tbx16, tbxta)
├── ZebraCHOP Command Line Batch Analysis Pipeline.pdf
└── LICENSE
```

## Configuration (`config_local.json`)

The wizard writes this automatically. Manual example:

```json
{
  "PATH": {
    "BOWTIE": "/usr/local/bin/bowtie",
    "TWOBITTOFA": "/usr/local/bin/twoBitToFa",
    "PRIMER3": "/usr/local/bin/primer3_core",
    "TWOBIT_INDEX_DIR": "/data/zebrachop/indexes",
    "BOWTIE_INDEX_DIR":  "/data/zebrachop/indexes",
    "GENE_TABLE_INDEX_DIR": "./genetable"
  },
  "THREADS": 4,
  "PYTHON2": "/usr/local/bin/python2"
}
```

`config_local.json` is git-ignored.

## Web UI

Run `bash frontend/run.sh` and open <http://127.0.0.1:8000>. Paste one or more gene names, pick the scoring model, click **Design guides** — the frontend spawns a Python-2 subprocess per batch, streams a job log, and renders each gene's ranked TSV as a sortable table you can download.

## Citing

If you use ZebraCHOP, please cite the upstream CHOPCHOP paper (Labun *et al.*, NAR 2019). The poison-exon-aware Ensembl post-processing in `reformat.py` is an addition specific to this fork.

## License

MIT — see [`LICENSE`](LICENSE). Portions of the code derive from CHOPCHOP (MIT) and from Microsoft's Azimuth (BSD 3-Clause); see [`NOTICE`](NOTICE) and in-file headers for attributions.
