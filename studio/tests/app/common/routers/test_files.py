from studio.app.common.routers.files import DirTreeGetter
from studio.app.common.schemas.files import TreeNode

workspace_id = "1"


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
