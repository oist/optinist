import uuid
from typing import Dict

from fastapi import APIRouter, BackgroundTasks

from studio.app.common.core.workflow.workflow import Message, NodeItem, RunItem
from studio.app.common.core.workflow.workflow_result import WorkflowResult
from studio.app.common.core.workflow.workflow_runner import WorkflowRunner

router = APIRouter()


@router.post("/run/{workspace_id}", response_model=str, tags=["run"])
async def run(workspace_id: str, runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())[:8]
    WorkflowRunner(workspace_id, unique_id, runItem).run_workflow(background_tasks)
    print("run snakemake")
    return unique_id


@router.post("/run/{workspace_id}/{uid}", response_model=str, tags=["run"])
async def run_id(
    workspace_id: str, uid: str, runItem: RunItem, background_tasks: BackgroundTasks
):
    WorkflowRunner(workspace_id, uid, runItem).run_workflow(background_tasks)
    print("run snakemake")
    print("forcerun list: ", runItem.forceRunList)
    return uid


@router.post(
    "/run/result/{workspace_id}/{uid}", response_model=Dict[str, Message], tags=["run"]
)
async def run_result(workspace_id: str, uid: str, nodeDict: NodeItem):
    return WorkflowResult(workspace_id, uid).get(nodeDict.pendingNodeIdList)
