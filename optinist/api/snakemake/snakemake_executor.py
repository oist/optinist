import os
from collections import deque
from typing import Dict

from snakemake import snakemake

from optinist.api.dir_path import DIRPATH
from optinist.api.logger import Logger
from optinist.api.snakemake.smk import SmkParam
from optinist.api.utils.filepath_creater import get_pickle_file, join_filepath
from optinist.api.workflow.workflow import Edge, Node


def snakemake_execute(unique_id: str, params: SmkParam):
    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=params.forceall,
        cores=params.cores,
        use_conda=params.use_conda,
        workdir=f"{os.path.dirname(DIRPATH.ROOT_DIR)}",
        log_handler=[Logger(unique_id).smk_logger],
    )


def delete_dependencies(
    unique_id: str,
    smk_params: SmkParam,
    nodeDict: Dict[str, Node],
    edgeDict: Dict[str, Edge],
):
    queue = deque()

    for param in smk_params.forcerun:
        queue.append(param.nodeId)

    while True:
        # terminate condition
        if len(queue) == 0:
            break

        # delete pickle
        node_id = queue.pop()
        algo_name = nodeDict[node_id].data.label

        pickle_filepath = join_filepath(
            [
                DIRPATH.OUTPUT_DIR,
                get_pickle_file(
                    unique_id=unique_id,
                    node_id=node_id,
                    algo_name=algo_name,
                ),
            ]
        )
        # print(pickle_filepath)

        if os.path.exists(pickle_filepath):
            os.remove(pickle_filepath)

        # 全てのedgeを見て、node_idがsourceならtargetをqueueに追加する
        for edge in edgeDict.values():
            if node_id == edge.source:
                queue.append(edge.target)
