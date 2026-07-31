#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${MPR_CONDA_ENV:-mutation_patrol}"
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ENV_FILE="$PROJECT_DIR/environment.yml"
REQUIREMENTS_FILE="$PROJECT_DIR/backend/requirements.txt"
STATE_DIR="$PROJECT_DIR/app_data"
STATE_FILE="$STATE_DIR/.environment.sha256"

load_conda() {
  if command -v conda >/dev/null 2>&1; then
    # shellcheck disable=SC1091
    source "$(conda info --base)/etc/profile.d/conda.sh"
    return 0
  fi

  local conda_profile
  for conda_profile in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/mambaforge/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh"
  do
    if [[ -f "$conda_profile" ]]; then
      # shellcheck disable=SC1090
      source "$conda_profile"
      return 0
    fi
  done

  return 1
}

environment_exists() {
  conda env list | awk 'NF && $1 !~ /^#/ {print $1}' | grep -Fxq "$ENV_NAME"
}

manifest_hash() {
  sha256sum "$ENV_FILE" "$REQUIREMENTS_FILE" \
    | sha256sum \
    | awk '{print $1}'
}

verify_runtime() {
  python -c "import fastapi, multipart, uvicorn" >/dev/null 2>&1 \
    && command -v npm >/dev/null 2>&1 \
    && command -v samtools >/dev/null 2>&1 \
    && command -v bcftools >/dev/null 2>&1 \
    && command -v minimap2 >/dev/null 2>&1
}

if ! load_conda; then
  echo "[ERROR] Conda was not found inside WSL."
  echo "Install Miniforge or Miniconda, reopen WSL, and start the app again."
  echo "Setup guide: $PROJECT_DIR/WSL_SETUP.md"
  return 1 2>/dev/null || exit 1
fi

if [[ ! -f "$ENV_FILE" ]]; then
  echo "[ERROR] Missing conda environment file: $ENV_FILE"
  return 1 2>/dev/null || exit 1
fi

mkdir -p "$STATE_DIR"
current_hash="$(manifest_hash)"
saved_hash="$(cat "$STATE_FILE" 2>/dev/null || true)"
needs_sync=0

if ! environment_exists; then
  echo "[SETUP] Creating conda environment: $ENV_NAME"
  conda env create -n "$ENV_NAME" -f "$ENV_FILE"
  needs_sync=1
fi

conda activate "$ENV_NAME"

if [[ "$current_hash" != "$saved_hash" ]] \
  || ! verify_runtime \
  || ! python "$PROJECT_DIR/scripts/check_environment.py" >/dev/null 2>&1
then
  needs_sync=1
fi

if [[ "$needs_sync" -eq 1 ]]; then
  echo "[SETUP] Synchronizing conda environment: $ENV_NAME"
  conda env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
  conda activate "$ENV_NAME"
  python -m pip install --disable-pip-version-check -r "$REQUIREMENTS_FILE"
  printf '%s\n' "$current_hash" > "$STATE_FILE"
fi

if ! verify_runtime || ! python "$PROJECT_DIR/scripts/check_environment.py"; then
  echo "[ERROR] Environment verification failed after synchronization."
  echo "Expected: Python/FastAPI, npm, samtools, bcftools, and minimap2."
  return 1 2>/dev/null || exit 1
fi

echo "[OK] Conda environment: $ENV_NAME"
echo "     Python:   $(python --version 2>&1)"
echo "     Node.js:  $(node --version)"
echo "     samtools: $(samtools --version | head -n 1)"
echo "     bcftools: $(bcftools --version | head -n 1)"
echo "     minimap2: $(minimap2 --version 2>&1 | head -n 1)"
