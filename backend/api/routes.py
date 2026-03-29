from fastapi import APIRouter
from pydantic import BaseModel
from celery_app import celery
from tasks import process_rmats_task

router = APIRouter()


class AnalysisQuery(BaseModel):
    path: str
    fdr: float = 0.05
    dpsi: float = 0.1
    include_jc: bool = True
    include_jcec: bool = True


@router.post("/analyze")
def start_analysis(query: AnalysisQuery):
    task = process_rmats_task.apply_async(
        args=[query.path, query.fdr, query.dpsi, query.include_jc, query.include_jcec]
    )
    return {"task_id": str(task.id), "status": "processing"}


@router.get("/task/{task_id}")
def get_status(task_id: str):
    task_result = celery.AsyncResult(task_id)
    if task_result.state == "FAILURE":
        return {
            "task_id": task_id,
            "status": "FAILED",
            "error": str(task_result.result),
        }
    if task_result.ready():
        result = task_result.result
        return {"task_id": task_id, "status": "COMPLETED", "result": result}
    return {"task_id": task_id, "status": task_result.state}
