#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

# This script activates the verified environment in the current shell.
# shellcheck disable=SC1091
source "$PROJECT_DIR/scripts/ensure_conda_env.sh"

cd frontend
node_platform="$(node -p "process.platform + '-' + process.arch")"
saved_node_platform="$(cat node_modules/.mpr-platform 2>/dev/null || true)"
if [[ ! -d node_modules ]] \
  || [[ package-lock.json -nt node_modules/.package-lock.json ]] \
  || [[ "$node_platform" != "$saved_node_platform" ]]
then
  echo "[SETUP] Synchronizing frontend packages..."
  npm ci
  printf '%s\n' "$node_platform" > node_modules/.mpr-platform
fi

if [[ ! -f dist/index.html ]] \
  || find src index.html vite.config.ts tsconfig.json -type f -newer dist/index.html -print -quit | grep -q .
then
  echo "[BUILD] Building the web interface..."
  if ! npm run build; then
    echo "[SETUP] Frontend build failed; reinstalling platform dependencies..."
    npm ci
    printf '%s\n' "$node_platform" > node_modules/.mpr-platform
    npm run build
  fi
else
  echo "[OK] Web interface is already up to date."
fi

cd "$PROJECT_DIR"
echo ""
echo "Mutation Patrol Robot is starting."
echo "Open this URL in your Windows browser:"
echo "  http://localhost:8000"
echo ""

APP_URL="${MPR_APP_URL:-http://localhost:8000}"
if command -v cmd.exe >/dev/null 2>&1; then
  (
    for _attempt in $(seq 1 120); do
      if python -c \
        "import urllib.request; urllib.request.urlopen('$APP_URL/api/health', timeout=1).read()" \
        >/dev/null 2>&1
      then
        cmd.exe /c start "" "$APP_URL" >/dev/null 2>&1
        exit 0
      fi
      sleep 1
    done
    echo "[WARN] Server did not become ready within two minutes."
  ) &
fi

exec uvicorn backend.app:app --host 127.0.0.1 --port 8000
