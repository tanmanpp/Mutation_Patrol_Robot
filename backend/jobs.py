# -*- coding: utf-8 -*-

import os
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone


MAX_WORKERS = max(1, int(os.environ.get("MPR_JOB_WORKERS", "1")))
_executor = ThreadPoolExecutor(max_workers=MAX_WORKERS, thread_name_prefix="mpr-job")
_lock = threading.Lock()
_jobs: dict[str, dict] = {}
_futures = {}


def _now():
    return datetime.now(timezone.utc).isoformat()


def _public_job(job: dict):
    return {
        key: value
        for key, value in job.items()
        if key not in {"callable", "args", "kwargs"}
    }


def _set_job(job_id: str, **changes):
    with _lock:
        _jobs[job_id].update(changes)


def _run_job(job_id: str):
    with _lock:
        job = _jobs[job_id]
        function = job.pop("callable")
        args = job.pop("args")
        kwargs = job.pop("kwargs")
        job.update(status="running", started_at=_now(), message="Job is running.")

    try:
        result = function(*args, **kwargs)
    except Exception as exc:
        _set_job(
            job_id,
            status="failed",
            finished_at=_now(),
            message=str(exc) or exc.__class__.__name__,
            error_type=exc.__class__.__name__,
        )
    else:
        _set_job(
            job_id,
            status="completed",
            finished_at=_now(),
            message="Job completed.",
            result=result,
        )


def submit_job(kind: str, label: str, function, *args, **kwargs):
    job_id = uuid.uuid4().hex[:12]
    job = {
        "job_id": job_id,
        "kind": kind,
        "label": label,
        "status": "queued",
        "message": "Job is queued.",
        "created_at": _now(),
        "started_at": None,
        "finished_at": None,
        "result": None,
        "callable": function,
        "args": args,
        "kwargs": kwargs,
    }
    with _lock:
        _jobs[job_id] = job
        _futures[job_id] = _executor.submit(_run_job, job_id)
    return get_job(job_id)


def get_job(job_id: str):
    with _lock:
        job = _jobs.get(job_id)
        return _public_job(dict(job)) if job else None


def list_jobs(limit: int = 50):
    with _lock:
        jobs = list(_jobs.values())[-limit:]
        return [_public_job(dict(job)) for job in reversed(jobs)]


def cancel_job(job_id: str):
    with _lock:
        job = _jobs.get(job_id)
        future = _futures.get(job_id)
        if not job:
            return None
        if job["status"] != "queued" or not future or not future.cancel():
            return _public_job(dict(job))
        job.update(
            status="cancelled",
            finished_at=_now(),
            message="Queued job was cancelled.",
        )
        return _public_job(dict(job))
