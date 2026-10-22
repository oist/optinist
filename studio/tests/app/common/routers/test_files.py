import asyncio
import dataclasses
import io
import json
import os
import threading
import time
from glob import glob

import numpy as np
import pytest
import tifffile
from fastapi import BackgroundTasks, UploadFile
from filelock import FileLock

from studio.app.common.core.utils.filelock_handler import FileLockUtils
from studio.app.common.routers import files as files_module
from studio.app.common.routers.files import (
    IMAGE_SHAPE_FILE,
    DirTreeGetter,
    get_files,
    get_image_shape_dict,
    read_image_shape,
    update_image_shape,
)
from studio.app.common.schemas.files import TreeNode
from studio.app.const import ACCEPT_FILE_EXT
from studio.app.dir_path import DIRPATH

workspace_id = "1"
TIFF = ACCEPT_FILE_EXT.TIFF_EXT.value
FIXTURE_INPUT_DIR = DIRPATH.INPUT_DIR


def test_create_files(client):
    response = client.get(f"/files/{workspace_id}?file_type=image")
    data = response.json()

    assert response.status_code == 200
    assert isinstance(data, list)
    assert len(data) > 0


def test_DirTreeGetter_tif():
    output = DirTreeGetter.get_tree(
        workspace_id, [".tif", ".tiff", ".TIF", ".TIFF"], "files"
    )
    assert len(output) == 4
    assert isinstance(output[0], TreeNode)


@pytest.fixture
def ws(tmp_path, monkeypatch):
    monkeypatch.setattr(DIRPATH, "INPUT_DIR", str(tmp_path))
    root = tmp_path / "ws"
    root.mkdir()
    return root


def _tif(path, shape, **kwargs):
    path.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        path, np.zeros(shape, dtype=np.uint8), photometric="minisblack", **kwargs
    )
    return path


def _as_dicts(nodes):
    return [dataclasses.asdict(n) for n in nodes]


def _cache(root):
    return json.loads((root / IMAGE_SHAPE_FILE).read_text())


def _glob_has_match(path, exts):
    """Reference implementation: one recursive glob per extension."""
    return any(glob(os.path.join(path, "**", f"*{e}"), recursive=True) for e in exts)


@pytest.mark.parametrize(
    "shape,kwargs",
    [
        ((4, 5), {}),
        ((3, 4, 5), {}),
        ((3, 4, 5), {"bigtiff": True}),
        ((2, 3, 4, 5), {"imagej": True}),
        ((2, 3, 4, 5), {"ome": True}),
    ],
)
def test_read_image_shape_matches_imread(tmp_path, shape, kwargs):
    p = _tif(tmp_path / "a.tif", shape, **kwargs)
    assert read_image_shape(str(p)) == list(tifffile.imread(str(p)).shape)


def test_read_image_shape_fixture_files():
    sample = os.path.join(FIXTURE_INPUT_DIR, workspace_id, "test.tif")
    empty = os.path.join(FIXTURE_INPUT_DIR, workspace_id, "files", "test1.tif")
    assert read_image_shape(sample) == list(tifffile.imread(sample).shape)
    assert read_image_shape(empty) == []
    assert read_image_shape(os.path.join(FIXTURE_INPUT_DIR, "missing.tif")) == []


def test_has_accepted_file_agrees_with_glob(tmp_path):
    root = tmp_path / "tree"
    (root / "empty" / "deeper").mkdir(parents=True)
    _tif(root / "deep" / "a" / "b" / "c" / "d" / "x.TIF", (2, 2))
    _tif(root / "mixed" / "y.tiff", (2, 2))
    (root / "none").mkdir()
    (root / "none" / "z.csv").write_text("a,b\n")
    _tif(root / "hidden" / ".h.tif", (2, 2))
    _tif(root / "hiddendir" / ".d" / "q.tif", (2, 2))
    (root / "dirnamed" / "foo.tif").mkdir(parents=True)
    outside = tmp_path / "outside"
    _tif(outside / "o.tif", (2, 2))
    (root / "linked").mkdir()
    os.symlink(outside, root / "linked" / "nas")
    (root / "broken").mkdir()
    os.symlink(tmp_path / "nowhere", root / "broken" / "dangling.tif")
    (root / "loop").mkdir()
    os.symlink(root / "loop", root / "loop" / "self")
    (root / "micro").mkdir()
    (root / "micro" / "a.thor.zip").write_bytes(b"")

    for sub in sorted(os.listdir(root)):
        path = str(root / sub)
        for exts in (TIFF, ACCEPT_FILE_EXT.MICROSCOPE_EXT.value):
            expected = _glob_has_match(path, exts)
            assert DirTreeGetter.has_accepted_file(path, exts) is expected, (sub, exts)


def test_has_accepted_file_ignores_glob_metacharacters(tmp_path):
    _tif(tmp_path / "br[ack]et" / "a.tif", (2, 2))
    assert DirTreeGetter.has_accepted_file(str(tmp_path / "br[ack]et"), TIFF)


def test_has_accepted_file_visits_symlink_cycle_once(tmp_path, monkeypatch):
    (tmp_path / "loop" / "inner").mkdir(parents=True)
    os.symlink(tmp_path / "loop", tmp_path / "loop" / "inner" / "back")
    scanned = []
    real_scandir = os.scandir
    monkeypatch.setattr(
        os, "scandir", lambda p: (scanned.append(p), real_scandir(p))[1]
    )

    assert DirTreeGetter.has_accepted_file(str(tmp_path / "loop"), TIFF) is False
    assert len(scanned) <= 2, scanned


def test_has_accepted_file_stops_at_first_match(tmp_path, monkeypatch):
    for i in range(20):
        (tmp_path / f"d{i}" / "deep" / "deeper").mkdir(parents=True)
    _tif(tmp_path / "top.tif", (2, 2))
    scanned = []
    real_scandir = os.scandir
    monkeypatch.setattr(
        os, "scandir", lambda p: (scanned.append(p), real_scandir(p))[1]
    )

    assert DirTreeGetter.has_accepted_file(str(tmp_path), TIFF) is True
    assert len(scanned) == 1, scanned


def test_get_tree_response_shape(ws):
    _tif(ws / "b.tif", (2, 3))
    _tif(ws / "sub" / "a.tif", (4, 5, 6))
    (ws / "empty").mkdir()
    (ws / "other").mkdir()
    (ws / "other" / "c.csv").write_text("a\n")

    assert _as_dicts(DirTreeGetter.get_tree("ws", TIFF)) == [
        {
            "path": "sub",
            "name": "sub",
            "isdir": True,
            "shape": None,
            "nodes": [
                {
                    "path": "sub/a.tif",
                    "name": "a.tif",
                    "isdir": False,
                    "nodes": [],
                    "shape": [4, 5, 6],
                }
            ],
        },
        {
            "path": "b.tif",
            "name": "b.tif",
            "isdir": False,
            "nodes": [],
            "shape": [2, 3],
        },
    ]
    assert _as_dicts(DirTreeGetter.get_tree("ws", ACCEPT_FILE_EXT.CSV_EXT.value)) == [
        {
            "path": "other",
            "name": "other",
            "isdir": True,
            "shape": None,
            "nodes": [
                {
                    "path": "other/c.csv",
                    "name": "c.csv",
                    "isdir": False,
                    "nodes": [],
                    "shape": None,
                }
            ],
        }
    ]
    assert not (ws / "other" / IMAGE_SHAPE_FILE).exists()


def test_get_tree_follows_symlinked_directory(ws, tmp_path):
    outside = tmp_path / "nas"
    _tif(outside / "o.tif", (2, 3))
    os.symlink(outside, ws / "link")

    tree = _as_dicts(DirTreeGetter.get_tree("ws", TIFF))
    assert tree == [
        {
            "path": "link",
            "name": "link",
            "isdir": True,
            "shape": None,
            "nodes": [
                {
                    "path": "link/o.tif",
                    "name": "o.tif",
                    "isdir": False,
                    "nodes": [],
                    "shape": [2, 3],
                }
            ],
        }
    ]


def test_get_tree_never_decodes_pixels(ws, monkeypatch):
    _tif(ws / "a.tif", (2, 3))

    def boom(*args, **kwargs):
        raise AssertionError("tifffile.imread called during tree walk")

    monkeypatch.setattr(tifffile, "imread", boom)
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]


def test_get_tree_writes_cache_once_per_request(ws, monkeypatch):
    for i in range(5):
        _tif(ws / f"d{i}" / f"f{i}.tif", (2, i + 1))
    writes = []
    real = files_module._write_json_atomic
    monkeypatch.setattr(
        files_module,
        "_write_json_atomic",
        lambda p, d: (writes.append(p), real(p, d)),
    )

    DirTreeGetter.get_tree("ws", TIFF)
    assert len(writes) == 1
    cache = _cache(ws)
    assert sorted(cache) == [f"d{i}/f{i}.tif" for i in range(5)]
    assert all({"shape", "mtime", "size"} <= set(v) for v in cache.values())

    DirTreeGetter.get_tree("ws", TIFF)
    assert len(writes) == 1


def test_get_tree_merges_into_entries_written_during_walk(ws, monkeypatch):
    _tif(ws / "a.tif", (2, 3))
    _tif(ws / "late.tif", (4, 5))
    (ws / IMAGE_SHAPE_FILE).write_text(
        json.dumps({"late.tif": {"shape": [4, 5]}, "gone.tif": {"shape": [1]}})
    )
    real = files_module.read_image_shape

    def read_and_interleave(path):
        if path.endswith("a.tif"):
            update_image_shape("ws", "late.tif")
        return real(path)

    monkeypatch.setattr(files_module, "read_image_shape", read_and_interleave)
    DirTreeGetter.get_tree("ws", TIFF)
    cache = _cache(ws)
    assert cache["a.tif"]["shape"] == [2, 3]
    assert cache["gone.tif"] == {"shape": [1]}
    assert cache["late.tif"]["shape"] == [4, 5] and "mtime" in cache["late.tif"]


def test_get_tree_does_not_overwrite_entry_replaced_during_walk(ws, monkeypatch):
    p = _tif(ws / "a.tif", (2, 3))
    real = files_module.read_image_shape

    def read_then_replace(path):
        shape = real(path)
        monkeypatch.setattr(files_module, "read_image_shape", real)
        _tif(p, (7, 8, 9))
        update_image_shape("ws", "a.tif")
        return shape

    monkeypatch.setattr(files_module, "read_image_shape", read_then_replace)
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]
    assert _cache(ws)["a.tif"]["shape"] == [7, 8, 9]


def test_get_tree_survives_cache_write_failure(ws, monkeypatch):
    _tif(ws / "a.tif", (2, 3))

    def fail(*args, **kwargs):
        raise PermissionError("read-only")

    monkeypatch.setattr(files_module, "_write_json_atomic", fail)
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]
    assert not (ws / IMAGE_SHAPE_FILE).exists()
    assert list(ws.glob("*.tmp")) == []


def test_get_tree_persists_shapes_read_before_a_failure(ws, monkeypatch):
    _tif(ws / "aa" / "a.tif", (2, 3))
    _tif(ws / "zz" / "b.tif", (4, 5))
    real = os.listdir
    monkeypatch.setattr(
        os,
        "listdir",
        lambda p: (_ for _ in ()).throw(OSError("stale handle"))
        if p.endswith("zz")
        else real(p),
    )

    with pytest.raises(OSError):
        DirTreeGetter.get_tree("ws", TIFF)
    assert _cache(ws)["aa/a.tif"]["shape"] == [2, 3]


def test_get_tree_detects_same_size_replacement(ws):
    p = _tif(ws / "a.tif", (2, 6))
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 6]
    size = p.stat().st_size

    _tif(p, (3, 4))
    assert p.stat().st_size == size
    os.utime(p, (p.stat().st_atime, p.stat().st_mtime + 10))
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [3, 4]


def test_get_tree_detects_same_mtime_replacement(ws):
    p = _tif(ws / "a.tif", (2, 3))
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]
    st = p.stat()

    _tif(p, (4, 5))
    os.utime(p, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [4, 5]


def test_concurrent_cache_updates_are_serialized(ws, monkeypatch):
    real = files_module._write_json_atomic
    monkeypatch.setattr(
        files_module,
        "_write_json_atomic",
        lambda p, d: (time.sleep(0.2), real(p, d))[1],
    )
    path = files_module._image_shape_path("ws")
    threads = [
        threading.Thread(
            target=files_module._locked_json_edit,
            args=(path, lambda d, k=k, i=i: d.update({k: {"shape": [i]}})),
        )
        for i, k in enumerate(("a.tif", "b.tif", "c.tif"))
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert sorted(_cache(ws)) == ["a.tif", "b.tif", "c.tif"]
    assert os.path.exists(FileLockUtils.get_lockfile_path(path))


def _image_shape_lock():
    path = files_module._image_shape_path("ws")
    return FileLock(FileLockUtils.get_lockfile_path(path))


def test_cache_write_waits_for_the_file_lock(ws):
    _tif(ws / "a.tif", (2, 3))
    writer = threading.Thread(target=update_image_shape, args=("ws", "a.tif"))

    with _image_shape_lock():
        writer.start()
        time.sleep(0.3)
        assert writer.is_alive()
        assert not (ws / IMAGE_SHAPE_FILE).exists()

    writer.join(timeout=5)
    assert not writer.is_alive()
    assert _cache(ws)["a.tif"]["shape"] == [2, 3]


def test_cache_write_gives_up_on_a_wedged_lock(ws, monkeypatch):
    monkeypatch.setattr(files_module, "CACHE_LOCK_TIMEOUT", 0.1)
    warnings = []
    monkeypatch.setattr(files_module.logger, "warning", warnings.append)
    _tif(ws / "a.tif", (2, 3))

    with _image_shape_lock():
        writer = threading.Thread(target=update_image_shape, args=("ws", "a.tif"))
        writer.start()
        writer.join(timeout=5)
        assert not writer.is_alive()

    assert _cache(ws)["a.tif"]["shape"] == [2, 3]
    assert any("cache lock timeout" in w for w in warnings)


def test_write_json_atomic_cleans_up_on_failure(tmp_path, monkeypatch):
    target = tmp_path / "cache.json"
    target.write_text("{}")
    monkeypatch.setattr(
        os, "replace", lambda *a: (_ for _ in ()).throw(OSError("no space"))
    )

    with pytest.raises(OSError):
        files_module._write_json_atomic(str(target), {"a": 1})
    assert target.read_text() == "{}"
    assert list(tmp_path.glob("*.tmp")) == []


def test_get_tree_trusts_legacy_cache_entries(ws, monkeypatch):
    _tif(ws / "a.tif", (2, 3))
    (ws / IMAGE_SHAPE_FILE).write_text(
        json.dumps({"a.tif": {"shape": [9, 9]}, "gone.tif": {"shape": [1]}})
    )
    monkeypatch.setattr(
        files_module, "_merge_image_shape_dict", lambda *a: pytest.fail("rewrote")
    )

    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [9, 9]


def test_get_tree_rereads_legacy_empty_shape(ws):
    _tif(ws / "a.tif", (2, 3))
    (ws / IMAGE_SHAPE_FILE).write_text(json.dumps({"a.tif": {"shape": []}}))

    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]
    assert _cache(ws)["a.tif"]["shape"] == [2, 3]


def test_get_tree_does_not_cache_an_unreadable_tiff(ws, monkeypatch):
    (ws / "bad.tif").write_bytes(b"not a tiff")

    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == []
    assert not (ws / IMAGE_SHAPE_FILE).exists()

    monkeypatch.setattr(
        files_module, "_merge_image_shape_dict", lambda *a: pytest.fail("rewrote")
    )
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == []


@pytest.mark.parametrize(
    "payload",
    ["{not json", "[]", '{"a.tif": "shape"}', '{"a.tif": {"shape": "x"}}'],
)
def test_get_tree_repairs_malformed_cache(ws, payload):
    _tif(ws / "a.tif", (2, 3))
    (ws / IMAGE_SHAPE_FILE).write_text(payload)

    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]
    assert _cache(ws)["a.tif"]["shape"] == [2, 3]


def test_update_image_shape_repairs_non_dict_cache(ws):
    _tif(ws / "a.tif", (2, 3))
    (ws / IMAGE_SHAPE_FILE).write_text("[]")

    assert update_image_shape("ws", "a.tif") == [2, 3]
    assert _cache(ws)["a.tif"]["shape"] == [2, 3]


def test_update_image_shape_entry_is_accepted_by_tree(ws, monkeypatch):
    _tif(ws / "a.tif", (2, 3))
    assert update_image_shape("ws", "a.tif") == [2, 3]
    assert {"shape", "mtime", "size"} <= set(get_image_shape_dict("ws")["a.tif"])

    monkeypatch.setattr(
        files_module, "_merge_image_shape_dict", lambda *a: pytest.fail("rewrote")
    )
    assert DirTreeGetter.get_tree("ws", TIFF)[0].shape == [2, 3]


def test_update_image_shape_missing_file_writes_nothing(ws):
    assert update_image_shape("ws", "nope.tif") == []
    assert not (ws / IMAGE_SHAPE_FILE).exists()


def test_get_files_does_not_block_event_loop(monkeypatch):
    def slow(*args, **kwargs):
        time.sleep(0.5)
        return []

    monkeypatch.setattr(DirTreeGetter, "get_tree", slow)

    async def run():
        task = asyncio.create_task(get_files("ws", "image"))
        t0 = time.monotonic()
        await asyncio.sleep(0.05)
        assert time.monotonic() - t0 < 0.3
        assert await task == []

    asyncio.run(run())


@pytest.mark.parametrize("handler", ["set_shape", "create_file"])
def test_shape_update_does_not_block_event_loop(ws, monkeypatch, handler):
    def slow(*args):
        time.sleep(0.5)
        return []

    monkeypatch.setattr(files_module, "update_image_shape", slow)
    monkeypatch.setattr(
        files_module.WorkspaceDataCapacityService, "is_available", lambda: False
    )
    call = {
        "set_shape": lambda: files_module.set_shape("ws", "a.tif"),
        "create_file": lambda: files_module.create_file(
            "ws",
            "a.tif",
            BackgroundTasks(),
            UploadFile(io.BytesIO(b"x"), filename="a.tif"),
            None,
        ),
    }[handler]

    async def run():
        task = asyncio.create_task(call())
        t0 = time.monotonic()
        await asyncio.sleep(0.05)
        assert time.monotonic() - t0 < 0.3
        await task

    asyncio.run(run())


def test_get_files_walks_on_the_dedicated_pool(monkeypatch):
    names = []
    monkeypatch.setattr(
        DirTreeGetter,
        "get_tree",
        lambda *a: names.append(threading.current_thread().name) or [],
    )
    assert asyncio.run(get_files("ws", "image")) == []
    assert names[0].startswith("filetree"), names


def test_get_files_unknown_type(monkeypatch):
    monkeypatch.setattr(
        DirTreeGetter, "get_tree", lambda *a, **k: pytest.fail("walked")
    )
    assert asyncio.run(get_files("ws", "bogus")) == []
    assert asyncio.run(get_files("ws", None)) == []
