import os
from pathlib import Path

from sqlalchemy.exc import NoResultFound
from sqlmodel import Session, delete, update

from studio.app.common.core.experiment.experiment_reader import ExptConfigReader
from studio.app.common.core.experiment.experiment_writer import ExptConfigWriter
from studio.app.common.core.logger import AppLogger
from studio.app.common.core.mode import MODE
from studio.app.common.core.utils.file_reader import get_folder_size
from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.common.db.database import session_scope
from studio.app.common.models.experiment import ExperimentRecord
from studio.app.common.models.workspace import Workspace
from studio.app.dir_path import DIRPATH

logger = AppLogger.get_logger()


class WorkspaceDataCapacityService:
    @classmethod
    def is_available(cls) -> bool:
        # The workspace data capcaticy feature is available in multiuser mode
        available = MODE.IS_MULTIUSER
        return available

    @classmethod
    def update_experiment_data_usage(cls, workspace_id: str, unique_id: str):
        workflow_dir = join_filepath([DIRPATH.OUTPUT_DIR, workspace_id, unique_id])
        if not os.path.exists(workflow_dir):
            logger.error(f"'{workflow_dir}' does not exist")
            return

        data_usage = get_folder_size(workflow_dir)

        cls._update_exp_data_usage_yaml(workspace_id, unique_id, data_usage)

        if cls.is_available():
            cls._update_exp_data_usage_db(workspace_id, unique_id, data_usage)

    @classmethod
    def _update_exp_data_usage_yaml(cls, workspace_id: str, unique_id: str, data_usage):
        # Read experiment config
        config = ExptConfigReader.read(workspace_id, unique_id)
        if not config:
            logger.error(f"[{workspace_id}/{unique_id}] does not exist")
            return

        # Make overwrite params
        update_params = {"data_usage": data_usage}

        # Overwrite experiment config
        ExptConfigWriter(workspace_id, unique_id).overwrite(update_params)

    @classmethod
    def _update_exp_data_usage_db(
        cls, workspace_id: str, unique_id: str, data_usage: int
    ):
        with session_scope() as db:
            try:
                exp = (
                    db.query(ExperimentRecord)
                    .filter(
                        ExperimentRecord.workspace_id == workspace_id,
                        ExperimentRecord.uid == unique_id,
                    )
                    .one()
                )
                exp.data_usage = data_usage
            except NoResultFound:
                exp = ExperimentRecord(
                    workspace_id=workspace_id,
                    uid=unique_id,
                    data_usage=data_usage,
                )
                db.add(exp)

    @classmethod
    def update_workspace_data_usage(
        cls, db: Session, workspace_id: str, auto_commit: bool = True
    ):
        workspace_dir = join_filepath([DIRPATH.INPUT_DIR, workspace_id])
        if not os.path.exists(workspace_dir):
            logger.error(f"'{workspace_dir}' does not exist")
            return

        input_data_usage = get_folder_size(workspace_dir)
        db.execute(
            update(Workspace)
            .where(Workspace.id == workspace_id)
            .values(input_data_usage=input_data_usage)
        )

        if auto_commit:
            db.commit()

    @classmethod
    def sync_workspace_data_capacity(
        cls, db: Session, workspace_id: str, delete_existing: bool = False
    ):
        """
        Sync workspace data usage and recalculate data capacity
        This is a convenience method that combines update_workspace_data_usage
        and recalculate_workspace_data_capacity

        Args:
            db: Database session
            workspace_id: The workspace ID to sync
            delete_existing: If True, delete all existing records
              before syncing (default: False)
        """
        cls.update_workspace_data_usage(db, workspace_id)
        cls.recalculate_workspace_data_capacity(
            db, workspace_id, delete_existing=delete_existing
        )

    @classmethod
    def recalculate_workspace_data_capacity(
        cls, db: Session, workspace_id: str, delete_existing: bool = False
    ):
        folder = join_filepath([DIRPATH.OUTPUT_DIR, workspace_id])
        if not os.path.exists(folder):
            logger.error(f"'{folder}' does not exist")
            return
        exp_records = []

        for exp_folder in Path(folder).iterdir():
            try:
                unique_id = exp_folder.name
                data_usage = get_folder_size(exp_folder.as_posix())

                cls._update_exp_data_usage_yaml(workspace_id, unique_id, data_usage)

                exp_records.append(
                    ExperimentRecord(
                        workspace_id=workspace_id,
                        uid=unique_id,
                        data_usage=data_usage,
                    )
                )
            except Exception as e:
                logger.error(
                    f"Failed to update Record information [{exp_folder}] [{e}]"
                )

        if cls.is_available():
            if delete_existing:
                # Delete all existing records and bulk insert new ones
                db.execute(
                    delete(ExperimentRecord).where(
                        ExperimentRecord.workspace_id == workspace_id
                    )
                )
                db.bulk_save_objects(exp_records)
                logger.info(
                    f"Deleted and recreated {len(exp_records)} experiment records "
                    f"for workspace [{workspace_id}]"
                )
            else:
                # Upsert: Update existing records or insert new ones
                for exp_record in exp_records:
                    try:
                        existing_record = (
                            db.query(ExperimentRecord)
                            .filter(
                                ExperimentRecord.workspace_id
                                == exp_record.workspace_id,
                                ExperimentRecord.uid == exp_record.uid,
                            )
                            .one()
                        )
                        # Update only data_usage field
                        existing_record.data_usage = exp_record.data_usage
                    except NoResultFound:
                        # Insert new record if it doesn't exist
                        db.add(exp_record)

                logger.info(
                    f"Updated/inserted {len(exp_records)} experiment records "
                    f"for workspace [{workspace_id}]"
                )

        logger.info(
            "Workspace capacity recalculation succeeded. "
            f"[workspace_id: {workspace_id}]"
        )
