# -*- coding: utf-8 -*-

import json
import subprocess
from pathlib import Path
from datetime import datetime


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def relative_path(path):
    path = Path(path)
    try:
        return path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(path)


def resolve_project_path(path):
    path = Path(path)
    if path.is_absolute():
        return path
    return PROJECT_ROOT / path


def make_json_safe(value):
    if isinstance(value, Path):
        return relative_path(value)
    if isinstance(value, str):
        path = Path(value)
        if path.is_absolute():
            return relative_path(path)
        return value
    if isinstance(value, dict):
        return {key: make_json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [make_json_safe(item) for item in value]
    return value


def setup_logger(log_path: Path):
    # 先用超簡單 logger，之後你要換成 logging 模組也行
    class SimpleLogger:
        def _write(self, level, msg):
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line = f"[{ts}] [{level}] {msg}"
            print(line)
            with open(log_path, "a", encoding="utf-8") as f:
                f.write(line + "\n")

        def info(self, msg): self._write("INFO", msg)
        def warning(self, msg): self._write("WARN", msg)
        def error(self, msg): self._write("ERROR", msg)

    # 清空舊 log（可選）
    log_path.write_text("", encoding="utf-8")
    return SimpleLogger()


def run_cmd(cmd, logger, dry_run=False, cwd=None):
    if isinstance(cmd, (list, tuple)):
        printable = " ".join(map(str, cmd))
    else:
        printable = str(cmd)

    logger.info(f"CMD: {printable}")
    if dry_run:
        return 0

    try:
        subprocess.run(cmd, check=True, cwd=cwd)
        return 0
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}")
        raise


def write_run_metadata(out_json: Path, args):
    payload = {
        "timestamp": datetime.now().isoformat(),
        "args": make_json_safe(vars(args)),
    }
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2, ensure_ascii=False)
