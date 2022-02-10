import h5py
from fastapi import APIRouter

router = APIRouter()


@router.get("/hdf5")
async def get_files():
    file_path = '/tmp/optinist/20220129-2131/20220129-2131.hdf5'

    # print_list = {}
    ds_dict = {}
    def get_ds_dictionaries(name, node):
        fullname = node.name
        if isinstance(node, h5py.Dataset):
            path_list = name.split('/')
            dict2leaf(ds_dict, path_list)

    def dict2leaf(ds_dict, path_list):
        name = path_list[0]
        if len(path_list) > 1:
            if name not in ds_dict.keys():
                ds_dict[name] = {
                    "isDir": True,
                }
            dict2leaf(ds_dict[name], path_list[1:])
        else:
            ds_dict[name] = {
                "isDir": False
            }

    with h5py.File(file_path, "r") as f:
        f.visititems(get_ds_dictionaries)

    return ds_dict
