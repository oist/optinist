from fastapi import APIRouter, BackgroundTasks

import uuid
from pydantic import BaseModel
from typing import List
from dataclasses import asdict

from optinist.workflow.set_file import ForceRun
from optinist.workflow.workflow import run_workflow
from optinist.workflow.results import get_results


router = APIRouter()


class RunItem(BaseModel):
    name: str = None
    nodeList: list = []
    edgeList: list = []
    snakemakeParam: dict = {}
    nwbParam: dict = {}
    forceRunList: List[ForceRun]


class NodeItem(BaseModel):
    pendingNodeIdList: list = []


@router.post("/run")
async def run(runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())
    run_workflow(unique_id, background_tasks, runItem)
    print("run snakemake")
    return unique_id


@router.post("/run/{uid}")
async def run_id(uid: str, runItem: RunItem, background_tasks: BackgroundTasks):
    run_workflow(uid, background_tasks, runItem)
    print("run snakemake")
    print("forcerun list: ", runItem.forceRunList)
    return uid


@router.post("/run/result/{uid}")
async def run_result(uid: str, nodeList: NodeItem):
    return get_results(uid, nodeList.pendingNodeIdList)
