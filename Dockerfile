FROM python:3.11.9-slim-bookworm

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    build-essential \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Install MMseqs2 (using SSE4.1 version instead of AVX2)
RUN wget https://mmseqs.com/latest/mmseqs-linux-sse41.tar.gz \
    && tar xvfz mmseqs-linux-sse41.tar.gz \
    && cp mmseqs/bin/mmseqs /usr/local/bin \
    && rm -rf mmseqs mmseqs-linux-sse41.tar.gz \
    && chmod +x /usr/local/bin/mmseqs

# Set working directory
WORKDIR /app

# Copy requirements first to leverage Docker cache
COPY requirements.txt .

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
COPY . .

# Expose port
EXPOSE 8080

# Run the application
CMD ["uvicorn", "tool:app", "--host", "0.0.0.0", "--port", "8080"]
