from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from celery_app import celery
import os

app = FastAPI(title="SpliceSight API")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

from api.routes import router
app.include_router(router, prefix="/api")

@app.get("/")
def read_root():
    return {"message": "SpliceSight Backend API is running."}
