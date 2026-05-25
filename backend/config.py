# -*- coding: utf-8 -*-

from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
APP_DATA_DIR = PROJECT_ROOT / "app_data"
DATABASES_DIR = APP_DATA_DIR / "databases"
LEGACY_DATABASES_DIR = PROJECT_ROOT / "databases"
UPLOADS_DIR = APP_DATA_DIR / "uploads"
RESULTS_DIR = APP_DATA_DIR / "results"


def ensure_app_dirs():
    for directory in [APP_DATA_DIR, DATABASES_DIR, UPLOADS_DIR, RESULTS_DIR]:
        directory.mkdir(parents=True, exist_ok=True)
