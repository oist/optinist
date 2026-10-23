import argparse

from sqlmodel import select

from studio.app.common.core.logger import AppLogger
from studio.app.common.core.workspace.workspace_data_capacity_services import (
    WorkspaceDataCapacityService,
)
from studio.app.common.db.database import session_scope
from studio.app.common.models.workspace import Workspace

logger = AppLogger.get_logger()


def main(args):
    if WorkspaceDataCapacityService.is_available():
        if args.delete_existing:
            logger.info(
                "Running with --delete-existing flag"
                " - will delete and recreate all records"
            )
        else:
            logger.info(
                "Running without --delete-existing flag - will update existing records"
            )

        with session_scope() as db:
            workspace_list = db.execute(
                select(Workspace.id).filter(Workspace.deleted.is_(False))
            ).scalars()
            for workspace_id in workspace_list:
                logger.info(
                    f"Syncing workspace data capacity for workspace: [{workspace_id}]"
                )
                WorkspaceDataCapacityService.sync_workspace_data_capacity(
                    db, str(workspace_id), delete_existing=args.delete_existing
                )
    else:
        # Single-user mode always uses delete_existing=True (default)
        WorkspaceDataCapacityService.recalculate_workspace_data_capacity(
            db=None, workspace_id="1", delete_existing=args.delete_existing
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sync workspace data capacity for all workspaces"
    )
    parser.add_argument(
        "--delete-existing",
        action="store_true",
        help="Delete all existing records before syncing. "
        "Without this flag, existing records will be updated (upsert).",
    )

    main(parser.parse_args())
