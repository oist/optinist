import os
from collections import deque
from typing import Dict

from snakemake import snakemake

from studio.app.common.core.logger import Logger
from studio.app.common.core.snakemake.smk import SmkParam
from studio.app.common.core.utils.filepath_creater import get_pickle_file, join_filepath
from studio.app.common.core.workflow.workflow import Edge, Node
from studio.app.dir_path import DIRPATH


def snakemake_execute(workspace_id: str, unique_id: str, params: SmkParam):
    logger = Logger(workspace_id, unique_id)
    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=params.forceall,
        cores=params.cores,
        use_conda=params.use_conda,
        workdir=f"{os.path.dirname(DIRPATH.STUDIO_DIR)}",
        configfiles=[
            join_filepath(
                [
                    DIRPATH.OUTPUT_DIR,
                    workspace_id,
                    unique_id,
                    DIRPATH.SNAKEMAKE_CONFIG_YML,
                ]
            )
        ],
        log_handler=[logger.smk_logger],
    )
    logger.clean_up()


def delete_dependencies(
    workspace_id: str,
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
                    workspace_id=workspace_id,
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
