#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
MODE="${1:-offline}"
OUTDIR="${2:-$HERE/results}"
RUNTIME="${RUNTIME:-local}"
PYTHON="${PYTHON:-python3}"
PIPE="$HERE/pipeline"

command -v nextflow >/dev/null || { echo "ERROR: nextflow not on PATH" >&2; exit 1; }

cd "$PIPE"
export NXF_DISABLE_CHECK_LATEST=true
export MPLCONFIGDIR="${MPLCONFIGDIR:-$HERE/.cache/mpl}"
export NUMBA_CACHE_DIR="${NUMBA_CACHE_DIR:-$HERE/.cache/numba}"
export PYTHONWARNINGS=ignore
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$OUTDIR"

if [[ "$RUNTIME" == docker ]]; then
  command -v docker >/dev/null && docker info >/dev/null 2>&1 || { echo "ERROR: Docker not running" >&2; exit 1; }
  docker image inspect scrna-scatac:end2end >/dev/null 2>&1 || docker build -t scrna-scatac:end2end .
  PROFILE=docker
else
  VENV="${VENV:-$HERE/pyenv}"
  if [[ ! -x "$VENV/bin/python" ]]; then
    echo ">> Creating venv from pipeline/requirements.txt"
    "$PYTHON" -m venv "$VENV"
    "$VENV/bin/pip" install -q --upgrade pip
    "$VENV/bin/pip" install -q -r requirements.txt
  fi
  export PATH="$VENV/bin:$PATH"
  PROFILE=local
fi

verify_inputs() {
  if command -v shasum >/dev/null; then SUM="shasum -a 256 -c"; else SUM="sha256sum -c"; fi
  ( cd "$1" && cp "$HERE/inputs/SHA256SUMS.txt" /tmp/_sha_$$ && $SUM /tmp/_sha_$$ >/dev/null ) \
    || { echo "ERROR: input checksum mismatch in $1 (see inputs/SHA256SUMS.txt)" >&2; rm -f /tmp/_sha_$$; exit 3; }
  rm -f /tmp/_sha_$$; echo ">> Input checksums OK ($1)"
}

case "$MODE" in
  download)
    echo ">> Mode: download inputs from GEO (GSE274934) + analyse"
    nextflow run main.nf -profile "$PROFILE" -resume --dataset geo_pilot --outdir "$OUTDIR"
    verify_inputs "$OUTDIR/inputs"
    ;;
  offline)
    echo ">> Mode: offline full cohort (RNA+TCR for 9 donors, ATAC where the library exists)"
    "$PYTHON" "$PIPE/bin/verify_local_inputs.py" --root "$HERE" --samples "$PIPE/data/gse274934_cohort.csv"
    nextflow run main.nf -profile "$PROFILE" -resume --dataset existing --run_edger true \
      --samples "$PIPE/data/gse274934_cohort.csv" --local_matrices "$HERE/inputs" --outdir "$OUTDIR"
    ;;
  *) echo "Usage: $0 [download|offline] [OUTDIR]" >&2; exit 2 ;;
esac

echo ">> Done. Open: $OUTDIR/analysis/index.html"
echo ">> Audit:      $OUTDIR/analysis/step_completeness_audit.csv"
