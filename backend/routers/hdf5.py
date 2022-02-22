import h5py
from fastapi import APIRouter

router = APIRouter()


@router.get("/hdf5/{file_path:path}")
async def get_files(file_path: str):

    ds_list = []
    def get_ds_dictionaries(name, node):
        fullname = node.name
        if isinstance(node, h5py.Dataset):
            path_list = name.split('/')
            if len(node.shape) != 0:
                get_dir_tree(ds_list, path_list, node, "")

    def get_dir_tree(nodes, path_list, node, path):
        name = path_list[0]
        if path == "":
            path = name
        else:
            path += "/" + name
        is_exists = False
        # 既にkeyがある
        for i, value in enumerate(nodes):
            if "name" in value.keys() and value["name"] == name:
                is_exists = True
                if len(path_list) > 1:
                    get_dir_tree(nodes[i]["nodes"], path_list[1:], node, path)

        if not is_exists:
            if len(path_list) > 1:
                nodes.append({
                    "isDir": True,
                    "name": name,
                    "path": path,
                    "nodes": [],
                })
                get_dir_tree(nodes[-1]["nodes"], path_list[1:], node, path)
            else:
                nodes.append({
                    "isDir": False,
                    "name": name,
                    "shape": node.shape,
                    "path": path,
                    "nbytes": f"{int(node.nbytes / (1000**2))} M",
                })

    with h5py.File(file_path, "r") as f:
        f.visititems(get_ds_dictionaries)

    return ds_list
