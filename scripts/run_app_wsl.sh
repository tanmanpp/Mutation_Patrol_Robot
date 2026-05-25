#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

if [[ ! -d ".venv" ]]; then
  python3 -m venv .venv
fi

source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r backend/requirements.txt

cd frontend
if [[ ! -d "node_modules" ]]; then
  npm install
else
  npm install
fi
npm run build

cd ..
echo ""
echo "Mutation Patrol Robot is starting."
echo "Open this URL in your Windows browser:"
echo "  http://localhost:8000"
echo ""
uvicorn backend.app:app --host 127.0.0.1 --port 8000
