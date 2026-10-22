import asyncio
import json
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from contextlib import suppress
from pathlib import PurePath
from typing import Callable, Dict, List, Optional
from urllib.parse import urlparse

import requests
import tifffile
from fastapi import APIRouter, BackgroundTasks, Depends, File, HTTPException, UploadFile
from filelock import FileLock, Timeout
from requests.models import Response
from sqlmodel import Session
from tqdm import tqdm

from studio.app.common.core.logger import AppLogger
from studio.app.common.core.utils.file_reader import JsonReader
from studio.app.common.core.utils.filelock_handler import FileLockUtils
from studio.app.common.core.utils.filepath_creater import (
    create_directory,
    join_filepath,
)
from studio.app.common.core.workspace.workspace_data_capacity_services import (
    WorkspaceDataCapacityService,
)
from studio.app.common.core.workspace.workspace_dependencies import (
    is_workspace_available,
    is_workspace_owner,
)
from studio.app.common.db.database import get_db
from studio.app.common.schemas.files import (
    DownloadFileRequest,
    DownloadStatus,
    FilePath,
    TreeNode,
)
from studio.app.const import ACCEPT_FILE_EXT, FILETYPE
from studio.app.dir_path import DIRPATH

router = APIRouter(prefix="/files", tags=["files"])

logger = AppLogger.get_logger()

IMAGE_SHAPE_FILE = ".image_shape.json"
CACHE_LOCK_TIMEOUT = 300

# Dedicated pool so slow NAS walks queue behind each other, not other executor work.
_tree_executor = ThreadPoolExecutor(max_workers=4, thread_name_prefix="filetree")


class DirTreeGetter:
    @classmethod
    def get_tree(
        cls,
        workspace_id,
        file_types: List[str],
        dirname: str = None,
        _shape_cache: dict = None,
    ) -> List[TreeNode]:
        nodes: List[TreeNode] = []

        if dirname is None:
            absolute_dirpath = join_filepath([DIRPATH.INPUT_DIR, workspace_id])
        else:
            absolute_dirpath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, dirname])

        if not os.path.exists(absolute_dirpath):
            return nodes

        is_image = file_types == ACCEPT_FILE_EXT.TIFF_EXT.value
        is_root = _shape_cache is None
        if is_root:
            _shape_cache = {
                "known": get_image_shape_dict(workspace_id) if is_image else {},
                "new": {},
            }

        try:
            sorted_listdir = sorted(
                os.listdir(absolute_dirpath),
                key=lambda x: (
                    not os.path.isdir(join_filepath([absolute_dirpath, x])),
                    x,
                ),
            )

            for node_name in sorted_listdir:
                if dirname is None:
                    relative_path = node_name
                else:
                    relative_path = join_filepath([dirname, node_name])

                search_dirpath = join_filepath([absolute_dirpath, node_name])

                if os.path.isfile(search_dirpath) and node_name.endswith(
                    tuple(file_types)
                ):
                    shape = (
                        cls._cached_shape(search_dirpath, relative_path, _shape_cache)
                        if is_image
                        else None
                    )
                    nodes.append(
                        TreeNode(
                            path=relative_path,
                            name=node_name,
                            isdir=False,
                            nodes=[],
                            shape=shape,
                        )
                    )
                elif os.path.isdir(search_dirpath) and cls.has_accepted_file(
                    search_dirpath, file_types
                ):
                    nodes.append(
                        TreeNode(
                            path=node_name,
                            name=node_name,
                            isdir=True,
                            nodes=cls.get_tree(
                                workspace_id, file_types, relative_path, _shape_cache
                            ),
                        )
                    )
        finally:
            if is_root and _shape_cache["new"]:
                _merge_image_shape_dict(
                    workspace_id, _shape_cache["known"], _shape_cache["new"]
                )

        return nodes

    @classmethod
    def _cached_shape(cls, filepath: str, relative_path: str, cache: dict):
        try:
            st = os.stat(filepath)
        except OSError:
            return []
        entry = cache["known"].get(relative_path)
        shape = _entry_shape(entry)
        # Legacy entries carry no mtime/size and are trusted until rewritten.
        if (
            shape
            and entry.get("mtime", st.st_mtime) == st.st_mtime
            and entry.get("size", st.st_size) == st.st_size
        ):
            return shape
        entry = _shape_entry(filepath, st)
        # A failed read is retried next walk rather than cached.
        if entry["shape"]:
            cache["new"][relative_path] = entry
        return entry["shape"]

    @classmethod
    def has_accepted_file(cls, path: str, file_types: List[str]) -> bool:
        # Mirrors glob("**/*ext"): hidden names skipped, directory names match too.
        exts = tuple(file_types)
        seen = {os.path.realpath(path)}
        for root, dirs, files in os.walk(path, followlinks=True):
            dirs[:] = [d for d in dirs if not d.startswith(".")]
            if any(
                name.endswith(exts) for name in files + dirs if not name.startswith(".")
            ):
                return True
            kept = []
            for d in dirs:
                real = os.path.realpath(os.path.join(root, d))
                if real not in seen:
                    seen.add(real)
                    kept.append(d)
            dirs[:] = kept
        return False


def _image_shape_path(workspace_id: str) -> str:
    return join_filepath([DIRPATH.INPUT_DIR, workspace_id, IMAGE_SHAPE_FILE])


def get_image_shape_dict(workspace_id: str) -> dict:
    return _read_json_dict(_image_shape_path(workspace_id))


def _read_json_dict(filepath: str) -> dict:
    """Missing, unparsable or non-object JSON reads as an empty cache."""
    try:
        data = JsonReader.read(filepath)
    except (FileNotFoundError, ValueError):
        return {}
    return data if isinstance(data, dict) else {}


def _entry_shape(entry) -> Optional[list]:
    shape = entry.get("shape") if isinstance(entry, dict) else None
    return shape if isinstance(shape, list) else None


def _merge_image_shape_dict(workspace_id: str, known: dict, new: dict) -> None:
    # An entry rewritten while the walk ran is newer than the walk's read of it.
    def merge(existing: dict) -> None:
        existing.update(
            {k: v for k, v in new.items() if existing.get(k) == known.get(k)}
        )

    try:
        _locked_json_edit(_image_shape_path(workspace_id), merge)
    except OSError as e:
        logger.warning(f"image shape cache not written ({workspace_id}): {e}")


def read_image_shape(filepath: str) -> list:
    """Shape of the first series from the TIFF header, without decoding pixels."""
    try:
        with tifffile.TiffFile(filepath) as tif:
            return list(tif.series[0].shape)
    except Exception:
        return []


def _shape_entry(filepath: str, st: os.stat_result) -> dict:
    return {
        "shape": read_image_shape(filepath),
        "mtime": st.st_mtime,
        "size": st.st_size,
    }


def update_image_shape(workspace_id: str, relative_file_path: str):
    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, relative_file_path])

    try:
        entry = _shape_entry(filepath, os.stat(filepath))
    except OSError:
        return []

    _locked_json_edit(
        _image_shape_path(workspace_id),
        lambda data: data.update({relative_file_path: entry}),
    )

    return entry["shape"]


def _locked_json_edit(filepath: str, edit: Callable[[dict], None]) -> None:
    """Read-modify-write a JSON cache under its file lock."""
    lock = FileLock(FileLockUtils.get_lockfile_path(filepath))
    try:
        lock.acquire(timeout=CACHE_LOCK_TIMEOUT)
    except Timeout:
        logger.warning(
            f"cache lock timeout after {CACHE_LOCK_TIMEOUT}s for {filepath}; "
            "proceeding without exclusive lock"
        )
    try:
        data = _read_json_dict(filepath)
        edit(data)
        _write_json_atomic(filepath, data)
    finally:
        lock.release()


def _write_json_atomic(filepath: str, data: dict) -> None:
    fd, tmp_path = tempfile.mkstemp(dir=os.path.dirname(filepath), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as f:
            json.dump(data, f, indent=2)
        os.replace(tmp_path, filepath)
    except BaseException:
        with suppress(OSError):
            os.unlink(tmp_path)
        raise


@router.get(
    "/{workspace_id}",
    response_model=List[TreeNode],
    dependencies=[Depends(is_workspace_available)],
)
async def get_files(workspace_id: str, file_type: str = None):
    file_types = {
        FILETYPE.IMAGE: ACCEPT_FILE_EXT.TIFF_EXT.value,
        FILETYPE.CSV: ACCEPT_FILE_EXT.CSV_EXT.value,
        FILETYPE.HDF5: ACCEPT_FILE_EXT.HDF5_EXT.value,
        FILETYPE.MICROSCOPE: ACCEPT_FILE_EXT.MICROSCOPE_EXT.value,
        FILETYPE.MATLAB: ACCEPT_FILE_EXT.MATLAB_EXT.value,
    }.get(file_type)
    if not file_types:
        return []
    return await asyncio.get_running_loop().run_in_executor(
        _tree_executor, DirTreeGetter.get_tree, workspace_id, file_types
    )


@router.post(
    "/{workspace_id}/shape/{filepath}",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def set_shape(workspace_id: str, filepath: str):
    try:
        await asyncio.get_running_loop().run_in_executor(
            None, update_image_shape, workspace_id, filepath
        )
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    return True


@router.post(
    "/{workspace_id}/upload/{filename}",
    response_model=FilePath,
    dependencies=[Depends(is_workspace_owner)],
)
async def create_file(
    workspace_id: str,
    filename: str,
    background_tasks: BackgroundTasks,
    file: UploadFile = File(...),
    db: Session = Depends(get_db),
):
    create_directory(join_filepath([DIRPATH.INPUT_DIR, workspace_id]))

    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, filename])

    with open(filepath, "wb") as f:
        shutil.copyfileobj(file.file, f)

    await asyncio.get_running_loop().run_in_executor(
        None, update_image_shape, workspace_id, filename
    )

    if WorkspaceDataCapacityService.is_available():
        background_tasks.add_task(
            WorkspaceDataCapacityService.update_workspace_data_usage, db, workspace_id
        )

    return {"file_path": filename}


DOWNLOAD_STATUS: Dict[str, DownloadStatus] = {}


@router.delete(
    "/{workspace_id}/delete/{filename:path}",
    response_model=bool,
    dependencies=[Depends(is_workspace_owner)],
)
async def delete_file(
    workspace_id: str,
    filename: str,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db),
):
    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, filename])
    if not os.path.exists(filepath):
        raise HTTPException(status_code=404, detail="File not found.")
    try:
        os.remove(filepath)

        if WorkspaceDataCapacityService.is_available():
            background_tasks.add_task(
                WorkspaceDataCapacityService.update_workspace_data_usage,
                db,
                workspace_id,
            )

        return True
    except Exception as e:
        logger.error(e, exc_info=True)
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/{workspace_id}/download/status",
    response_model=DownloadStatus,
    dependencies=[Depends(is_workspace_available)],
)
async def get_download_status(workspace_id: str, file_name: str):
    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, file_name])
    try:
        return DOWNLOAD_STATUS[filepath]
    except:  # noqa
        raise HTTPException(status_code=404)


@router.post(
    "/{workspace_id}/download",
    dependencies=[Depends(is_workspace_owner)],
)
async def download_file(
    workspace_id: str,
    file: DownloadFileRequest,
    background_tasks: BackgroundTasks,
    db: Session = Depends(get_db),
):
    path = PurePath(urlparse(file.url).path)
    if path.suffix not in ACCEPT_FILE_EXT.ALL_EXT.value:
        raise HTTPException(status_code=400, detail="Invalid url")

    create_directory(join_filepath([DIRPATH.INPUT_DIR, workspace_id]))

    try:
        res = requests.get(file.url, stream=True)
        res.raise_for_status()
    except Exception as e:
        raise HTTPException(status_code=422, detail=str(e))
    background_tasks.add_task(download, db, res, path.name, workspace_id)
    background_tasks.add_task(update_image_shape, workspace_id, path.name)
    return {"file_name": path.name}


def download(
    db: Session, res: Response, file_name: str, workspace_id: str, chunk_size=1024
):
    total = int(res.headers.get("content-length", 0))
    filepath = join_filepath([DIRPATH.INPUT_DIR, workspace_id, file_name])
    current = 0

    try:
        with open(filepath, "wb") as file, tqdm(
            desc=filepath,
            total=total,
            unit="iB",
            unit_scale=True,
            unit_divisor=1024,
        ) as bar:
            for data in res.iter_content(chunk_size=chunk_size):
                size = file.write(data)
                current += size
                DOWNLOAD_STATUS[filepath] = DownloadStatus(total=total, current=current)
                bar.update(size)
    except Exception as e:
        DOWNLOAD_STATUS[filepath] = DownloadStatus(error=str(e))

    if WorkspaceDataCapacityService.is_available():
        WorkspaceDataCapacityService.update_workspace_data_usage(db, workspace_id)
