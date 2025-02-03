#!/bin/sh
echo "INFO     Unique Gene Finder - $VERSION"
uvicorn tool:app --host 0.0.0.0 --port 8080 --proxy-headers
