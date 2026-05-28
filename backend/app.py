# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Annotated

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .config import UPLOADS_DIR, ensure_app_dirs
from . import services
from .config import PROJECT_ROOT


app = FastAPI(title="Mutation Patrol Robot API", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

def has_upload(upload: UploadFile | None):
    return upload is not None and bool(upload.filename)

FRONTEND_DIST = PROJECT_ROOT / "frontend" / "dist"
if FRONTEND_DIST.exists():
    app.mount("/assets", StaticFiles(directory=FRONTEND_DIST / "assets"), name="assets")


@app.on_event("startup")
def startup():
    ensure_app_dirs()


@app.get("/api/health")
def health():
    return {"ok": True}


@app.get("/api/databases")
def list_databases():
    return {"databases": services.list_databases()}


@app.get("/api/databases/{db_name}/genes")
def database_genes(db_name: str):
    try:
        return {"genes": services.get_database_genes(db_name)}
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.post("/api/databases")
async def create_database(
    name: Annotated[str, Form()],
    ref_fasta: Annotated[UploadFile, File()],
    annotations: Annotated[list[UploadFile], File()],
    genes: Annotated[str | None, Form()] = None,
    force: Annotated[bool, Form()] = True,
):
    try:
        db_name = services.safe_name(name)
        upload_dir = UPLOADS_DIR / "databases" / db_name
        ref_path = await services.save_upload(ref_fasta, upload_dir / ref_fasta.filename)
        annotation_paths = []
        for annotation in annotations:
            if not annotation.filename.lower().endswith(".jsonl"):
                continue
            annotation_paths.append(await services.save_upload(annotation, upload_dir / annotation.filename))
        if not annotation_paths:
            raise ValueError("No .jsonl annotation files were uploaded.")
        gene_list = [item.strip() for item in genes.split(",") if item.strip()] if genes else None
        return services.build_database(db_name, ref_path, annotation_paths, gene_list, force)
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/analyze")
async def analyze_sample(
    run_name: Annotated[str, Form()],
    db_name: Annotated[str, Form()],
    fastq: Annotated[UploadFile | None, File()] = None,
    bam: Annotated[UploadFile | None, File()] = None,
    min_depth: Annotated[int, Form()] = 10,
    threads: Annotated[int, Form()] = 8,
    min_mapq: Annotated[int, Form()] = 20,
    min_alt_count: Annotated[int, Form()] = 1,
    min_alt_freq: Annotated[float, Form()] = 0.05,
    candidate_source: Annotated[str, Form()] = "bam",
    force: Annotated[bool, Form()] = True,
):
    try:
        fastq = fastq if has_upload(fastq) else None
        bam = bam if has_upload(bam) else None
        if (fastq is None) == (bam is None):
            raise ValueError("Upload exactly one of FASTQ or BAM.")
        run_id = services.safe_name(run_name)
        upload_dir = UPLOADS_DIR / "samples" / run_id
        fastq_path = await services.save_upload(fastq, upload_dir / fastq.filename) if fastq else None
        bam_path = await services.save_upload(bam, upload_dir / bam.filename) if bam else None
        return services.analyze_sample(
            run_name=run_id,
            db_name=db_name,
            fastq_path=fastq_path,
            bam_path=bam_path,
            ref_fasta=None,
            min_depth=min_depth,
            threads=threads,
            min_mapq=min_mapq,
            min_alt_count=min_alt_count,
            min_alt_freq=min_alt_freq,
            candidate_source=candidate_source,
            force=force,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.delete("/api/uploads/samples")
def clear_uploaded_samples():
    try:
        return services.clear_uploaded_samples()
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/uploads/samples")
def sample_upload_stats():
    return services.sample_upload_stats()


@app.post("/api/site-query")
async def site_query(
    analysis_run: Annotated[str, Form()],
    db_name: Annotated[str, Form()],
    gene: Annotated[str, Form()],
    query_type: Annotated[str, Form()],
    query_value: Annotated[str, Form()],
    alt: Annotated[str | None, Form()] = None,
    min_alt_freq: Annotated[float, Form()] = 0.05,
    min_mapq: Annotated[int, Form()] = 20,
    force: Annotated[bool, Form()] = True,
):
    try:
        if query_type not in {"genomic_pos", "cds_pos", "aa_pos"}:
            raise ValueError("query_type must be genomic_pos, cds_pos, or aa_pos.")
        return services.query_sites(
            analysis_run=analysis_run,
            db_name=db_name,
            gene=gene,
            query_type=query_type,
            query_value=query_value,
            alt=alt,
            min_alt_freq=min_alt_freq,
            min_mapq=min_mapq,
            force=force,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/results")
def results():
    return {"results": services.list_results()}


@app.get("/api/results/{run_id}")
def result(run_id: str):
    return services.get_result(run_id)


@app.get("/api/results/{run_id}/report")
def result_html_report(run_id: str):
    try:
        report_path = services.generate_html_report(run_id)
        return FileResponse(report_path, media_type="text/html", filename=report_path.name)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/")
def serve_ui_root():
    index_path = FRONTEND_DIST / "index.html"
    if not index_path.exists():
        raise HTTPException(
            status_code=503,
            detail="Frontend is not built yet. Run npm run build in frontend/.",
        )
    return FileResponse(index_path)


@app.get("/{path:path}")
def serve_ui_fallback(path: str):
    if path.startswith("api/"):
        raise HTTPException(status_code=404, detail="API route not found")
    index_path = FRONTEND_DIST / "index.html"
    if not index_path.exists():
        raise HTTPException(
            status_code=503,
            detail="Frontend is not built yet. Run npm run build in frontend/.",
        )
    return FileResponse(index_path)
