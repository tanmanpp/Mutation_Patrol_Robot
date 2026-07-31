#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

# shellcheck disable=SC1091
source "$PROJECT_DIR/scripts/ensure_conda_env.sh"

exec uvicorn backend.app:app --host 0.0.0.0 --port 8000 --reload
