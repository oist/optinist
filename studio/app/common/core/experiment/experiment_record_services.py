from sqlalchemy.exc import NoResultFound
from sqlmodel import Session, delete

from studio.app.common.core.logger import AppLogger
from studio.app.common.core.mode import MODE
from studio.app.common.db.database import session_scope
from studio.app.common.models.experiment import ExperimentRecord

logger = AppLogger.get_logger()


class ExperimentRecordService:
    @classmethod
    def is_available(cls) -> bool:
        # ExperimentRecordService is available in multiuser mode
        available = MODE.IS_MULTIUSER
        return available

    @classmethod
    def regist_record_on_workflow_completed(cls, workspace_id: str, unique_id: str):
        """
        Processing upon workflow completion
        """

        # Update ExperimentRecord to database
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

                # Currently nothing
                # *Add data to be registered in the future as needed
                pass

            except NoResultFound:
                exp = ExperimentRecord(
                    workspace_id=workspace_id,
                    uid=unique_id,
                    # *Add data to be registered in the future as needed
                )
                db.add(exp)

    @classmethod
    def delete_record(
        cls, db: Session, workspace_id: str, unique_id: str, auto_commit: bool = False
    ):
        db.execute(
            delete(ExperimentRecord).where(
                ExperimentRecord.workspace_id == workspace_id,
                ExperimentRecord.uid == unique_id,
            )
        )

        if auto_commit:
            db.commit()

    @classmethod
    def copy_record(
        cls,
        db: Session,
        workspace_id: str,
        unique_id: str,
        new_unique_id: str,
        _new_name: str,
        auto_commit: bool = False,
    ):
        try:
            exp = (
                db.query(ExperimentRecord)
                .filter(
                    ExperimentRecord.workspace_id == workspace_id,
                    ExperimentRecord.uid == unique_id,
                )
                .one()
            )
            new_exp = ExperimentRecord(
                workspace_id=workspace_id,
                uid=new_unique_id,
                data_usage=exp.data_usage,
            )
            db.add(new_exp)

            if auto_commit:
                db.commit()

        except NoResultFound:
            # If it fails roll back the transaction
            logger.error(
                f"Experiment {unique_id} not found in workspace {workspace_id}"
            )
