import os
from collections import deque

from snakemake import snakemake
from snakemake.exceptions import print_exception
from snakemake.logging import logger
from snakemake.workflow import Workflow
from snakemake.dag import DAG
from snakemake.persistence import Persistence

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam
from optinist.api.utils.filepath_creater import join_filepath


class SmkExecutor:
    def __init__(self, snakefile, forceall=False, cores=2, use_conda=False):
        self.snakefile = os.path.abspath(snakefile)
        self.logger = logger
        self.logger.setup_logfile()
        self.forceall = forceall
        self.cores = cores
        self.use_conda = use_conda

    def init_workflow(self):
        self.workflow = Workflow(
            snakefile=self.snakefile,
            cores=self.cores,
            use_conda=self.use_conda,
            overwrite_configfiles=[],
        )

        self.workflow.include(
            self.snakefile,
            overwrite_default_target=True,
            print_compilation=False,
        )
        self.workflow.check()

    def init_dag(self):
        self.init_workflow()
        targetrules = map(
            self.workflow._rules.__getitem__,
            filter(self.workflow.is_rule, ["all"])
        )

        self.dag = DAG(
            self.workflow,
            self.workflow.rules,
            targetfiles=[],
            targetrules=targetrules,
            forceall=self.forceall,
        )

        self.workflow.persistence = Persistence(
            nolock=True,
            dag=self.dag,
            conda_prefix=self.workflow.conda_prefix,
            singularity_prefix=self.workflow.singularity_prefix,
            shadow_prefix=self.workflow.shadow_prefix,
        )
        self.dag.init()

    def init_graph(self):
        self.init_dag()
        graph = {}
        for job in self.dag.jobs:
            graph[job.rule] = [
                dep.rule for dep in self.dag.dependencies[job]
            ]

        # node ids
        ids = {node: i for i, node in enumerate(graph)}

        file_graph = {}
        for job in self.dag.jobs:
            file_graph[ids[job.rule]] = [
                "".join(files) for files in self.dag.dependencies[job].values()
            ]

        # calculate edges
        edge_list = [
            [ids[dep], ids[node]]
            for node, deps in graph.items()
            for dep in deps
        ]

        edge_dict = create_edge_dict(edge_list)
        return edge_dict, file_graph

    def execute(self, forcerun):
        self.init_graph()
        try:
            success = self.workflow.execute(
                forceall=self.forceall,
                forcerun=forcerun,
                nolock=True,
                unlock=False,
                updated_files=[],
                
            )
        except BrokenPipeError:
            success = False
        except (Exception, BaseException) as ex:
            if "workflow" in locals():
                print_exception(ex, self.workflow.linemaps)
            else:
                print_exception(ex, dict())
            success = False
            assert False, "Snakemake execute error"

        if "workflow" in locals() and self.workflow.persistence:
            self.workflow.persistence.unlock()

        self.logger.cleanup()

        return success


def snakemake_execute(params: SmkParam):
    # run snakemake
    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=params.forceall,
        cores=params.cores,
        use_conda=params.use_conda,
        workdir=f"{os.path.dirname(DIRPATH.ROOT_DIR)}",
    )


def delete_dependencies(params):
    """
        [[1, 0], [2, 1], [3, 2]]
        1  →  0
        ['/Users/shogoakiyama/Desktop/optinist/optinist/test_data/snakemake/1/suite2p_file_convert.pkl']
        →  ['/Users/shogoakiyama/Desktop/optinist/optinist/test_data/snakemake/2/suite2p_roi.pkl']
        2  →  1
        ['/Users/shogoakiyama/Desktop/optinist/optinist/test_data/snakemake/0/data_endoscope.pkl']
        →  ['/Users/shogoakiyama/Desktop/optinist/optinist/test_data/snakemake/1/suite2p_file_convert.pkl']
        3  →  2
        []  →  ['/Users/shogoakiyama/Desktop/optinist/optinist/test_data/snakemake/0/data_endoscope.pkl']
    """

    del_filepath_list = [
        join_filepath([DIRPATH.OUTPUT_DIR, x])
        for x in params.forcerun
    ]

    smk_executor = SmkExecutor(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=params.forceall,
        cores=params.cores,
    )
    edge_dict, file_graph = smk_executor.init_graph()

    queue = deque()
    for key, value_list in file_graph.items():
        for value in value_list:
            if value in del_filepath_list:
                queue.append(key)

    while True:
        # terminate condition
        if len(queue) == 0:
            break
        
        # delete item path
        del_key = queue.pop()
        del_filepath = file_graph[del_key]

        for filepath in del_filepath:
            print(filepath)
            if os.path.exists(filepath):
                os.remove(filepath)

        # push
        if del_key in edge_dict:
            for push_key in edge_dict[del_key]:
                queue.append(push_key)


def create_edge_dict(edges):
    edge_dict = {}
    for e in edges:
        if e[0] not in edge_dict:
            edge_dict[e[0]] = []
        edge_dict[e[0]].append(e[1])
    return edge_dict
