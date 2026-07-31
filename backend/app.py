# -*- coding: utf-8 -*-

from pathlib import Path
from typing import Annotated

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

from .config import UPLOADS_DIR, ensure_app_dirs
from . import jobs, services
from .config import PROJECT_ROOT


app = FastAPI(title="Mutation Patrol Robot API", version="0.5.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5173", "http://127.0.0.1:5173"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

def has_upload(upload: UploadFile | None):
    return upload is not None and bool(upload.filename)


def require_upload_extension(upload: UploadFile, extensions: tuple[str, ...], label: str):
    filename = (upload.filename or "").lower()
    if not any(filename.endswith(extension) for extension in extensions):
        expected = ", ".join(extensions)
        raise ValueError(f"{label} must use one of these extensions: {expected}")


def validate_analysis_parameters(
    min_depth: int,
    threads: int,
    min_mapq: int,
    min_baseq: int,
    min_alt_count: int,
    min_alt_freq: float,
    candidate_source: str,
):
    if not 1 <= threads <= 64:
        raise ValueError("threads must be between 1 and 64.")
    if min_depth < 0 or min_alt_count < 0:
        raise ValueError("Depth and alternate allele count must be >= 0.")
    if not 0 <= min_mapq <= 255:
        raise ValueError("min_mapq must be between 0 and 255.")
    if not 0 <= min_baseq <= 93:
        raise ValueError("min_baseq must be between 0 and 93.")
    if not 0 <= min_alt_freq <= 1:
        raise ValueError("min_alt_freq must be between 0 and 1.")
    if candidate_source not in {"bam", "vcf"}:
        raise ValueError("candidate_source must be bam or vcf.")


FRONTEND_DIST = PROJECT_ROOT / "frontend" / "dist"
if FRONTEND_DIST.exists():
    app.mount("/assets", StaticFiles(directory=FRONTEND_DIST / "assets"), name="assets")


@app.on_event("startup")
def startup():
    ensure_app_dirs()


@app.get("/api/health")
def health():
    return {"api_ok": True, **services.runtime_health()}


@app.get("/api/jobs")
def list_jobs():
    return {"jobs": jobs.list_jobs()}


@app.get("/api/jobs/{job_id}")
def get_job(job_id: str):
    job = jobs.get_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    return job


@app.delete("/api/jobs/{job_id}")
def cancel_job(job_id: str):
    job = jobs.cancel_job(job_id)
    if not job:
        raise HTTPException(status_code=404, detail=f"Job not found: {job_id}")
    return job


@app.get("/api/databases")
def list_databases():
    return {"databases": services.list_databases()}


@app.get("/api/databases/{db_name}/genes")
def database_genes(db_name: str):
    try:
        return {"genes": services.get_database_genes(db_name)}
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.post("/api/databases", status_code=202)
async def create_database(
    name: Annotated[str, Form()],
    ref_fasta: Annotated[UploadFile, File()],
    annotations: Annotated[list[UploadFile], File()],
    genes: Annotated[str | None, Form()] = None,
    force: Annotated[bool, Form()] = True,
):
    try:
        require_upload_extension(ref_fasta, (".fa", ".fasta", ".fna"), "Reference FASTA")
        db_name = services.safe_name(name)
        upload_dir = UPLOADS_DIR / "databases" / db_name
        ref_path = await services.save_upload(ref_fasta, upload_dir)
        annotation_paths = []
        for index, annotation in enumerate(annotations):
            if not (annotation.filename or "").lower().endswith(".jsonl"):
                continue
            annotation_paths.append(
                await services.save_upload(
                    annotation,
                    upload_dir,
                    prefix=f"annotation_{index:03d}_",
                )
            )
        if not annotation_paths:
            raise ValueError("No .jsonl annotation files were uploaded.")
        gene_list = [item.strip() for item in genes.split(",") if item.strip()] if genes else None
        return jobs.submit_job(
            "database",
            db_name,
            services.build_database,
            db_name,
            ref_path,
            annotation_paths,
            gene_list,
            force,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/analyze", status_code=202)
async def analyze_sample(
    run_name: Annotated[str, Form()],
    db_name: Annotated[str, Form()],
    fastq: Annotated[UploadFile | None, File()] = None,
    bam: Annotated[UploadFile | None, File()] = None,
    min_depth: Annotated[int, Form()] = 10,
    threads: Annotated[int, Form()] = 8,
    min_mapq: Annotated[int, Form()] = 20,
    min_baseq: Annotated[int, Form()] = 20,
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
        validate_analysis_parameters(
            min_depth,
            threads,
            min_mapq,
            min_baseq,
            min_alt_count,
            min_alt_freq,
            candidate_source,
        )
        if fastq:
            require_upload_extension(fastq, (".fastq", ".fastq.gz", ".fq", ".fq.gz"), "FASTQ")
        if bam:
            require_upload_extension(bam, (".bam",), "BAM")
        run_id = services.safe_name(run_name)
        upload_dir = UPLOADS_DIR / "samples" / run_id
        fastq_path = await services.save_upload(fastq, upload_dir) if fastq else None
        bam_path = await services.save_upload(bam, upload_dir) if bam else None
        return jobs.submit_job(
            "analysis",
            run_id,
            services.analyze_sample,
            run_name=run_id,
            db_name=db_name,
            fastq_path=fastq_path,
            bam_path=bam_path,
            ref_fasta=None,
            min_depth=min_depth,
            threads=threads,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
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


@app.post("/api/site-query", status_code=202)
async def site_query(
    analysis_run: Annotated[str, Form()],
    db_name: Annotated[str, Form()],
    gene: Annotated[str, Form()],
    query_type: Annotated[str, Form()],
    query_value: Annotated[str | None, Form()] = None,
    alt: Annotated[str | None, Form()] = None,
    min_alt_freq: Annotated[float, Form()] = 0.05,
    min_mapq: Annotated[int, Form()] = 20,
    min_depth: Annotated[int, Form()] = 10,
    min_baseq: Annotated[int, Form()] = 20,
    min_allele_count: Annotated[int, Form()] = 2,
    force: Annotated[bool, Form()] = True,
):
    try:
        if query_type not in {"gene_region", "genomic_pos", "cds_pos", "aa_pos"}:
            raise ValueError(
                "query_type must be gene_region, genomic_pos, cds_pos, or aa_pos."
            )
        if not 0 <= min_alt_freq <= 1:
            raise ValueError("min_alt_freq must be between 0 and 1.")
        if not 0 <= min_mapq <= 255:
            raise ValueError("min_mapq must be between 0 and 255.")
        if min_depth < 1:
            raise ValueError("min_depth must be at least 1.")
        if not 0 <= min_baseq <= 93:
            raise ValueError("min_baseq must be between 0 and 93.")
        if min_allele_count < 1:
            raise ValueError("min_allele_count must be at least 1.")
        if query_type != "gene_region":
            services.parse_query_values(query_value or "")
        label = f"{services.safe_name(analysis_run)}:{services.safe_name(gene)}"
        return jobs.submit_job(
            "site_query",
            label,
            services.query_sites,
            analysis_run=analysis_run,
            db_name=db_name,
            gene=gene,
            query_type=query_type,
            query_value=query_value,
            alt=alt,
            min_alt_freq=min_alt_freq,
            min_mapq=min_mapq,
            min_depth=min_depth,
            min_baseq=min_baseq,
            min_allele_count=min_allele_count,
            force=force,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/coverage", status_code=202)
async def coverage(
    analysis_run: Annotated[str, Form()],
    db_name: Annotated[str, Form()],
    minimum_call_depth: Annotated[int, Form()] = 10,
    min_mapq: Annotated[int, Form()] = 20,
    min_baseq: Annotated[int, Form()] = 20,
    force: Annotated[bool, Form()] = True,
):
    try:
        if minimum_call_depth < 1:
            raise ValueError("minimum_call_depth must be at least 1.")
        if not 0 <= min_mapq <= 255:
            raise ValueError("min_mapq must be between 0 and 255.")
        if not 0 <= min_baseq <= 93:
            raise ValueError("min_baseq must be between 0 and 93.")
        label = f"{services.safe_name(analysis_run)}:coverage"
        return jobs.submit_job(
            "coverage",
            label,
            services.generate_coverage,
            analysis_run=analysis_run,
            db_name=db_name,
            minimum_call_depth=minimum_call_depth,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            force=force,
        )
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/results")
def results():
    return {"results": services.list_results()}


@app.get("/api/results/{run_id}")
def result(run_id: str):
    try:
        return services.get_result(run_id)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/results/{run_id}/coverage/{gene}")
def gene_coverage(run_id: str, gene: str):
    try:
        return services.get_gene_coverage(run_id, gene)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/results/{run_id}/site-query/{query_id}/complete-table")
def complete_gene_table(run_id: str, query_id: str):
    try:
        table_path = services.get_complete_gene_table_path(run_id, query_id)
        filename = (
            f"{services.safe_name(run_id)}_"
            f"{services.safe_name(query_id)}_complete_gene_table.csv"
        )
        return FileResponse(
            table_path,
            media_type="text/csv",
            filename=filename,
        )
    except FileNotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


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
