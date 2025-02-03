FROM ghcr.io/soedinglab/mmseqs2 AS mmseqs

FROM python:3.11.9-slim-bookworm

WORKDIR /app

# Copy MMseqs2 from official image
COPY --from=mmseqs /usr/local/bin/entrypoint /usr/local/bin/mmseqs

# Install minimal dependencies
RUN apt-get update && apt-get install -y \
    libatomic1 \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements and install Python dependencies
COPY requirements.txt .
RUN pip install -U pip && pip install -r requirements.txt

# Copy service files
COPY tool.py run.sh ./
RUN chmod +x run.sh

# VERSION INFORMATION
ARG VERSION ???
ENV VERSION=$VERSION

# Command to run
ENTRYPOINT ["/app/run.sh"]
