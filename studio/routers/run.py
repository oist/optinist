from fastapi import APIRouter, BackgroundTasks
import uuid

from studio.api.workflow.workflow import NodeItem, RunItem
from studio.api.workflow.workflow_runner import WorkflowRunner
from studio.api.workflow.workflow_result import WorkflowResult

router = APIRouter()


@router.post("/run")
async def run(runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())[:8]
    WorkflowRunner(unique_id, runItem).run_workflow(background_tasks)
    print("run snakemake")
    return unique_id


@router.post("/run/{uid}")
async def run_id(uid: str, runItem: RunItem, background_tasks: BackgroundTasks):
    WorkflowRunner(uid, runItem).run_workflow(background_tasks)
    print("run snakemake")
    print("forcerun list: ", runItem.forceRunList)
    return uid


@router.post("/run/result/{uid}")
async def run_result(uid: str, nodeDict: NodeItem):
    return WorkflowResult(uid).get(nodeDict.pendingNodeIdList)
