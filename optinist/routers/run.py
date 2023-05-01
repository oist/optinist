from typing import Dict
from fastapi import APIRouter, BackgroundTasks
import uuid

from optinist.api.workflow.workflow import NodeItem, RunItem, Message
from optinist.api.workflow.workflow_runner import WorkflowRunner
from optinist.api.workflow.workflow_result import WorkflowResult

router = APIRouter()


@router.post("/run", response_model=str, tags=['run'])
async def run(runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())[:8]
    WorkflowRunner(unique_id, runItem).run_workflow(background_tasks)
    print("run snakemake")
    return unique_id


@router.post("/run/{uid}", response_model=str, tags=['run'])
async def run_id(uid: str, runItem: RunItem, background_tasks: BackgroundTasks):
    WorkflowRunner(uid, runItem).run_workflow(background_tasks)
    print("run snakemake")
    print("forcerun list: ", runItem.forceRunList)
    return uid


@router.post("/run/result/{uid}", response_model=Dict[str, Message], tags=['run'])
async def run_result(uid: str, nodeDict: NodeItem):
    return WorkflowResult(uid).get(nodeDict.pendingNodeIdList)
