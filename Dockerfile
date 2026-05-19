# ZebraCHOP container: Python 2.7 engine + Python 3 frontend + bowtie + UCSC tools.
#
# Build:
#     docker build -t zebrachop .
# Run (mount your reference indexes at /refs):
#     docker run -p 8000:8000 -v /path/to/indexes:/refs zebrachop
#
FROM debian:bullseye-slim

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv3/bin:/usr/local/bin:${PATH}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        python2 python2-dev python-setuptools \
        python3 python3-pip python3-venv \
        build-essential curl ca-certificates \
        bowtie \
    && rm -rf /var/lib/apt/lists/*

# UCSC twoBitToFa
RUN curl -fsSL https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -o /usr/local/bin/twoBitToFa \
    && chmod +x /usr/local/bin/twoBitToFa

# get-pip for Python 2.7 (Debian dropped it)
RUN curl -fsSL https://bootstrap.pypa.io/pip/2.7/get-pip.py -o /tmp/get-pip.py \
    && python2 /tmp/get-pip.py \
    && rm /tmp/get-pip.py

WORKDIR /app
COPY requirements.txt ./
RUN python2 -m pip install --no-cache-dir -r requirements.txt

# Python 3 venv for the frontend
RUN python3 -m venv /opt/venv3
COPY frontend/requirements.txt /tmp/req3.txt
RUN /opt/venv3/bin/pip install --no-cache-dir --upgrade pip \
    && /opt/venv3/bin/pip install --no-cache-dir -r /tmp/req3.txt

COPY . .

# Generate a default config_local.json that points at the in-container tools.
RUN echo '{ \
    "PATH": { \
        "BOWTIE": "/usr/bin/bowtie", \
        "TWOBITTOFA": "/usr/local/bin/twoBitToFa", \
        "PRIMER3": "", \
        "TWOBIT_INDEX_DIR": "/refs", \
        "BOWTIE_INDEX_DIR": "/refs", \
        "GENE_TABLE_INDEX_DIR": "/app/genetable" \
    }, \
    "THREADS": 4, \
    "PYTHON2": "/usr/bin/python2" \
}' > /app/config_local.json

VOLUME ["/refs", "/app/frontend/data"]
EXPOSE 8000

CMD ["/opt/venv3/bin/python", "-m", "uvicorn", "frontend.main:app", "--host", "0.0.0.0", "--port", "8000"]
