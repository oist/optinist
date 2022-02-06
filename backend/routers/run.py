import os
import uuid
from glob import glob

from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import List

from workflow.params import get_typecheck_params
from workflow.set_workflow import set_workflow
from workflow.results import get_results
from cui_api.snakemake import run_snakemake
from .const import BASE_DIR


router = APIRouter()


class RunItem(BaseModel):
    nodeList: list = []
    edgeList: list = []
    snakemakeParam: dict = {}
    nwbParam: dict = {}


@router.post("/run")
async def params(runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())

    set_workflow(unique_id, runItem)
    
    snakemake_params = get_typecheck_params(runItem.snakemakeParam, "snakemake")
    background_tasks.add_task(run_snakemake, snakemake_params)

    print("run snakemake")

    return unique_id


@router.post("/run/result/{uid}")
async def params(uid: str):
    results = get_results(uid)
    return results
