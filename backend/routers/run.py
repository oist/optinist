import os
import uuid
from glob import glob

from fastapi import APIRouter, BackgroundTasks
from pydantic import BaseModel
from typing import List

from workflow.params import get_typecheck_params
from workflow.set_workflow import set_workflow
from workflow.results import get_results

from cui_api.run import run_snakemake
from cui_api.const import BASE_DIR


router = APIRouter()

class ForceRun(BaseModel):
    nodeId: str
    name: str


class RunItem(BaseModel):
    name: str = None
    nodeList: list = []
    edgeList: list = []
    snakemakeParam: dict = {}
    nwbParam: dict = {}
    forceRunList: List[ForceRun]


def run_workflow(unique_id, background_tasks, runItem):
    set_workflow(unique_id, runItem)

    snakemake_params = get_typecheck_params(runItem.snakemakeParam, "snakemake")
    snakemake_params["forcerun"] = get_forcerun_list(runItem.forceRunList)
    background_tasks.add_task(run_snakemake, snakemake_params)


@router.post("/run")
async def params(runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())
    run_workflow(unique_id, background_tasks, runItem)
    print("run snakemake")

    return unique_id


@router.post("/run/{uid}")
async def params(uid: str, runItem: RunItem, background_tasks: BackgroundTasks):
    run_workflow(uid, background_tasks, runItem)
    print("run snakemake")

    return uid


class NodeItem(BaseModel):
    pendingNodeIdList: list = []


@router.post("/run/result/{uid}")
async def params(uid: str, nodeList: NodeItem):
    results = get_results(uid, nodeList.pendingNodeIdList)
    return results
