"""
Celery tasks module — imported by both the FastAPI app and the Celery worker,
ensuring the task is registered in both processes.
"""

from celery_app import celery
from core.rmats_parser import parse_rmats_files


@celery.task(name="process_rmats_data")
def process_rmats_task(
    path: str,
    fdr: float,
    dpsi: float,
    include_jc: bool = True,
    include_jcec: bool = True,
):
    results = parse_rmats_files(path, fdr, dpsi, include_jc, include_jcec)
    return {"status": "success", "data": results}
