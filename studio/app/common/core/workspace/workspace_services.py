import os
import shutil

from fastapi import HTTPException, status
from sqlmodel import Session

from studio.app.common.core.experiment.experiment_services import ExperimentService
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.models.workspace import Workspace
from studio.app.dir_path import DIRPATH

logger = AppLogger.get_logger()


class WorkspaceService:
    @classmethod
    def delete_workspace_contents(
        cls,
        db: Session,
        ws: Workspace,
    ):
        workspace_id = str(ws.id)
        logger.info(f"Deleting workspace data for workspace '{workspace_id}'")

        workspace_dir = join_filepath([DIRPATH.OUTPUT_DIR, workspace_id])
        deleted_statuses = []

        # Delete experiment folders under workspace
        if os.path.exists(workspace_dir):
            for experiment_id in os.listdir(workspace_dir):
                # Skip hidden files and directories
                if experiment_id.startswith("."):
                    continue

                deleted_status = ExperimentService.delete_experiment(
                    db, workspace_id, experiment_id, auto_commit=False
                )

                deleted_statuses.append(deleted_status)
        else:
            logger.warning(f"Workspace directory '{workspace_dir}' does not exist")

        if all(deleted_statuses):
            # Delete the workspace directory itself
            cls.delete_workspace_files(workspace_id=workspace_id)

            # Delete input directory
            cls.delete_workspace_files(workspace_id=workspace_id, is_input_dir=True)

            # Soft delete the workspace
            ws.deleted = True
        else:
            # Throw Exception if data was not deleted
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Failed to delete workspace '{workspace_id}'",
            )

    @classmethod
    def delete_workspace_files(cls, workspace_id: str, is_input_dir: bool = False):
        if is_input_dir:
            directory = join_filepath([DIRPATH.INPUT_DIR, workspace_id])
        else:
            directory = join_filepath([DIRPATH.OUTPUT_DIR, workspace_id])
        try:
            if os.path.exists(directory):
                shutil.rmtree(directory)
                logger.info(f"Deleted directory: {directory}")
            else:
                logger.warning(f"'{directory}' already deleted or never existed")
        except Exception as e:
            logger.error(
                f"Failed to delete directory '{directory}': {e}",
                exc_info=True,
            )

    @classmethod
    def process_workspace_deletion(cls, db: Session, workspace_id: str, user_id: str):
        try:
            # Search for workspace
            ws: Workspace = (
                db.query(Workspace)
                .filter(
                    Workspace.id == workspace_id,
                    Workspace.user_id == user_id,
                    Workspace.deleted.is_(False),
                )
                .first()
            )

            if not ws:
                raise HTTPException(status_code=404, detail="Workspace not found")

            # Delete workspace storage files
            cls.delete_workspace_contents(db, ws)

            # Commit all DB changes before doing anything irreversible
            db.commit()
        except Exception as e:
            db.rollback()
            logger.error(
                "Error deleting or updating workspace %s: %s",
                workspace_id,
                e,
                exc_info=True,
            )
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail=f"Failed to delete or update workspace {workspace_id}.",
            )
