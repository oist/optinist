import uuid
from typing import Dict

from fastapi import APIRouter, BackgroundTasks, Depends

from studio.app.common.core.workflow.workflow import Message, NodeItem, RunItem
from studio.app.common.core.workflow.workflow_result import WorkflowResult
from studio.app.common.core.workflow.workflow_runner import WorkflowRunner
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
    is_workspace_owner,
)

router = APIRouter(prefix="/run", tags=["run"])


@router.post(
    "/{workspace_id}",
    response_model=str,
    dependencies=[Depends(is_workspace_owner)],
)
async def run(workspace_id: str, runItem: RunItem, background_tasks: BackgroundTasks):
    unique_id = str(uuid.uuid4())[:8]
    WorkflowRunner(workspace_id, unique_id, runItem).run_workflow(background_tasks)
    print("run snakemake")
    return unique_id


@router.post(
    "/{workspace_id}/{uid}",
    response_model=str,
    dependencies=[Depends(is_workspace_owner)],
)
async def run_id(
    workspace_id: str, uid: str, runItem: RunItem, background_tasks: BackgroundTasks
):
    WorkflowRunner(workspace_id, uid, runItem).run_workflow(background_tasks)
    print("run snakemake")
    print("forcerun list: ", runItem.forceRunList)
    return uid


@router.post(
    "/result/{workspace_id}/{uid}",
    response_model=Dict[str, Message],
    dependencies=[Depends(is_workspace_available)],
)
async def run_result(workspace_id: str, uid: str, nodeDict: NodeItem):
    return WorkflowResult(workspace_id, uid).get(nodeDict.pendingNodeIdList)


@router.post(
    "/cancel/{workspace_id}/{uid}",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def run_cancel(workspace_id: str, uid: str):
    return WorkflowResult(workspace_id, uid).cancel()
