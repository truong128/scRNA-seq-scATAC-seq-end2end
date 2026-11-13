#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
DATASET="${1:-geo_pilot}"
if ! command -v nextflow >/dev/null; then echo 'ERROR: Nextflow not on PATH' >&2; exit 1; fi
if ! command -v docker >/dev/null; then echo 'ERROR: Docker not on PATH' >&2; exit 1; fi
if ! docker info >/dev/null 2>&1; then echo 'ERROR: Docker daemon not running' >&2; exit 1; fi
if ! docker image inspect scrna-scatac:end2end >/dev/null 2>&1; then
    docker build -t scrna-scatac:end2end .
fi
mkdir -p results
nextflow run main.nf -profile docker -resume --dataset "$DATASET" "${@:2}"
