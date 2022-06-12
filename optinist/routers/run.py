from fastapi import APIRouter, BackgroundTasks
import uuid

from optinist.api.workflow.workflow import NodeItem, RunItem
from optinist.api.workflow.workflow_runner import WorkflowRunner
from optinist.api.workflow.workflow_result import WorkflowResult

router = APIRouter()


@router.post("/run")
async def run(runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())[:8]
    WorkflowRunner.run_workflow(unique_id, background_tasks, runItem)
    print("run snakemake")
    return unique_id


@router.post("/run/{uid}")
async def run_id(uid: str, runItem: RunItem, background_tasks: BackgroundTasks):
    WorkflowRunner.run_workflow(uid, background_tasks, runItem)
    print("run snakemake")
    print("forcerun list: ", runItem.forceRunList)
    return uid


@router.post("/run/result/{uid}")
async def run_result(uid: str, nodeDict: NodeItem):
    return WorkflowResult(uid).get(nodeDict.pendingNodeIdList)
