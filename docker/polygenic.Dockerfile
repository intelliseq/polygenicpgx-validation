# polygenic (polygenic-pgx) image for the validation harness.
# Self-contained: installs the published PyPI package and vendors the PGx models
# (the models are not shipped inside the pip package). Build from the REPO ROOT:
#   docker build -f docker/polygenic.Dockerfile -t pgxbench/polygenic:local .
# Pinned to py3.8 (the package's clean-install target).
FROM python:3.8-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
        tabix bcftools gcc g++ make zlib1g-dev libbz2-dev liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir polygenic-pgx==2.5.24 \
    && python3 -c "import polygenic; from polygenic import pgstk"

# PGx star-allele models (PharmVar-derived + CPIC TPMT + VKORC1 rs9923231).
COPY models/pgx /opt/models/pgx
ENV POLYGENIC_MODELS_DIR=/opt/models/pgx
