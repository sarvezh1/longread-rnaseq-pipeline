FROM python:3.11-slim@sha256:0dd364ba7e10242f07755449e3a3d0e35f9efd987952737b90def6709ab0c5ce

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates curl procps \
 && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir numpy==1.26.4 scipy==1.13.1 pandas==2.2.3 statsmodels==0.14.4 scikit-learn==1.5.2

RUN curl -fsSL https://github.com/comprna/SUPPA/archive/refs/tags/v2.4.tar.gz -o /tmp/suppa.tar.gz \
 && mkdir -p /opt/suppa \
 && tar -xzf /tmp/suppa.tar.gz --strip-components=1 -C /opt/suppa \
 && sed -i 's/from statsmodels.sandbox.stats.multicomp import multipletests/from statsmodels.stats.multitest import multipletests/' /opt/suppa/lib/diff_tools.py \
 && printf '#!/bin/sh\nexec python /opt/suppa/suppa.py "$@"\n' > /usr/local/bin/suppa.py \
 && chmod +x /usr/local/bin/suppa.py \
 && rm /tmp/suppa.tar.gz \
 && suppa.py --version

LABEL org.opencontainers.image.title="SUPPA" \
      org.opencontainers.image.version="2.4" \
      org.opencontainers.image.source="https://github.com/comprna/SUPPA/releases/tag/v2.4"
