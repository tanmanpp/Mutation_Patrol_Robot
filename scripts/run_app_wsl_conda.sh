#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${MPR_CONDA_ENV:-mutation_patrol}"
PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"

load_conda() {
  if command -v conda >/dev/null 2>&1; then
    return 0
  fi

  for conda_profile in \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/mambaforge/etc/profile.d/conda.sh"
  do
    if [[ -f "$conda_profile" ]]; then
      # shellcheck disable=SC1090
      source "$conda_profile"
      return 0
    fi
  done

  return 1
}

if ! load_conda; then
  echo "Conda was not found in WSL."
  echo "Install Miniconda/Miniforge in WSL first, then run this script again."
  echo "Example: https://docs.conda.io/projects/miniconda/"
  exit 1
fi

if ! conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "Creating conda environment: $ENV_NAME"
  conda create -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda \
    python=3.11 pip nodejs samtools bcftools minimap2
fi

conda activate "$ENV_NAME"

cd "$PROJECT_DIR"

if ! command -v npm >/dev/null 2>&1; then
  echo "npm was not found after activating conda environment: $ENV_NAME"
  echo "Installing nodejs into conda environment: $ENV_NAME"
  conda install -y -n "$ENV_NAME" -c conda-forge nodejs
  conda activate "$ENV_NAME"
fi

if ! command -v npm >/dev/null 2>&1; then
  echo "npm is still unavailable after installing nodejs."
  echo "Please check the conda installation inside WSL."
  exit 1
fi

python -m pip install -r backend/requirements.txt

cd frontend
npm install
npm run build

cd "$PROJECT_DIR"
echo ""
echo "Mutation Patrol Robot is starting."
echo "Open this URL in your Windows browser:"
echo "  http://localhost:8000"
echo ""

uvicorn backend.app:app --host 127.0.0.1 --port 8000
