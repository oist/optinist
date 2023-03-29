import os

from fastapi import HTTPException, status
from snakemake import snakemake

from optinist.api.config.config_writer import ConfigWriter
from optinist.api.dir_path import DIRPATH
from optinist.api.edit_ROI.wrappers import edit_roi_wrapper_dict


class EditROI:
    def __init__(self, action: str, filepath: str, params: dict) -> None:
        self.action = action
        self.filepath = filepath
        self.params = params
        self.algo = self.get_algo()
        self.node_dirpath = os.path.dirname(self.filepath)


    def get_algo(self):
        filepath = self.filepath
        
        if (
            not os.path.exists(filepath)
            and os.path.commonpath([DIRPATH.OPTINIST_DIR, filepath])
            != DIRPATH.OPTINIST_DIR
        ):
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND)
        if os.path.basename(filepath) != 'cell_roi.json':
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)

        algo_list = edit_roi_wrapper_dict.keys()

        algo = next((algo for algo in algo_list if algo in filepath), None)
        if not algo:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST)
        return algo


    def excute(self):
        self.set_smk_config()

        snakemake(
            DIRPATH.SNAKEMAKE_FILEPATH,
            use_conda=True,
            cores=2,
            workdir=f"{os.path.dirname(DIRPATH.ROOT_DIR)}",
        )
        
    def set_smk_config(self):
        config = {
            'action': self.action,
            'node_dirpath': self.node_dirpath,
            'algo': self.algo,
            'params': self.params
        }
        ConfigWriter.write(
            dirname=DIRPATH.ROOT_DIR,
            filename=DIRPATH.SNAKEMAKE_CONFIG_YML,
            config=config
        )
