import h5py
from fastapi import APIRouter

router = APIRouter()


@router.get("/hdf5")
async def get_files():
    file_path = '/tmp/optinist/20220129-2131/20220129-2131.hdf5'

    # print_list = {}
    ds_list = []
    def get_ds_dictionaries(name, node):
        fullname = node.name
        if isinstance(node, h5py.Dataset):
            path_list = name.split('/')
            get_dir_tree(ds_list, path_list, node)

    def get_dir_tree(nodes, path_list, node):
        name = path_list[0]
        is_exists = False
        # 既にkeyがある
        for i, value in enumerate(nodes):
            if "name" in value.keys() and value["name"] == name:
                is_exists = True
                if len(path_list) > 1:
                    get_dir_tree(nodes[i]["nodes"], path_list[1:], node)

        if not is_exists:
            if len(path_list) > 1:
                nodes.append({
                    "isDir": True,
                    "name": name,
                    "nodes": [],
                })
                get_dir_tree(nodes[-1]["nodes"], path_list[1:], node)
            else:
                nodes.append({
                    "isDir": False,
                    "name": name,
                    "shape": node.shape,
                })

    with h5py.File(file_path, "r") as f:
        f.visititems(get_ds_dictionaries)

    return ds_list
