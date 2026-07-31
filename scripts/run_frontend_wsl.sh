#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

# shellcheck disable=SC1091
source "$PROJECT_DIR/scripts/ensure_conda_env.sh"

cd frontend
node_platform="$(node -p "process.platform + '-' + process.arch")"
saved_node_platform="$(cat node_modules/.mpr-platform 2>/dev/null || true)"
if [[ ! -d node_modules ]] \
  || [[ package-lock.json -nt node_modules/.package-lock.json ]] \
  || [[ "$node_platform" != "$saved_node_platform" ]]
then
  npm ci
  printf '%s\n' "$node_platform" > node_modules/.mpr-platform
fi

exec npm run dev -- --host 0.0.0.0
