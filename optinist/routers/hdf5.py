import h5py
from fastapi import APIRouter

from optinist.routers.model import HDF5Node

router = APIRouter()


def get_hdf5_file(filepath):
    hdf5_list = []
    def get_ds_dictionaries(path, node):
        if isinstance(node, h5py.Dataset):
            if len(node.shape) != 0:
                recursive_dir_tree(hdf5_list, path.split('/'), node, "")

    def recursive_dir_tree(nodes, path_list, node, parent_path):
        name = path_list[0]
        path = name if parent_path == "" else f"{parent_path}/{name}"

        is_exists = False
        # 既にkeyがある
        for i, value in enumerate(nodes):
            if value.name == name:
                is_exists = True
                if len(path_list) > 1:
                    recursive_dir_tree(nodes[i].nodes, path_list[1:], node, path)

        if not is_exists:
            if len(path_list) > 1:
                nodes.append(
                    HDF5Node(
                        isDir=True,
                        name=name,
                        path=path,
                        nodes=[],
                    )
                )
                recursive_dir_tree(nodes[-1].nodes, path_list[1:], node, path)
            else:
                nodes.append(
                    HDF5Node(
                        isDir=False,
                        name=name,
                        path=path,
                        shape=node.shape,
                        nbytes=f"{int(node.nbytes / (1000**2))} M",
                    )
                )

    with h5py.File(filepath, "r") as f:
        f.visititems(get_ds_dictionaries)

    return hdf5_list


@router.get("/hdf5/{file_path:path}")
async def get_files(file_path: str):
    print(file_path)
    return get_hdf5_file(file_path)
