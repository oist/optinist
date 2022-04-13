import os

from snakemake.exceptions import print_exception
from snakemake.logging import logger
from snakemake.workflow import Workflow
from snakemake.dag import DAG
from snakemake.persistence import Persistence

from optinist.api.dir_path import DIRPATH
from optinist.api.snakemake.smk import SmkParam


class SmkExecutor:
    def __init__(self, snakefile, forceall=False):
        self.snakefile = os.path.abspath(snakefile)
        self.logger = logger
        self.logger.setup_logfile()
        self.forceall = forceall

    def init_workflow(self, cores):
        self.workflow = Workflow(
            snakefile=self.snakefile,
            cores=cores,
        )

        self.workflow.include(
            self.snakefile,
            overwrite_default_target=True,
            print_compilation=False,
        )
        self.workflow.check()

    def init_dag(self):
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
                files for files in self.dag.dependencies[job].values()
            ]

        # calculate edges
        edges = [
            [ids[dep], ids[node]]
            for node, deps in graph.items()
            for dep in deps
        ]
        print(edges)
        for x in edges:
            print(x[0], " → ", x[1])
            print(file_graph[x[0]], " → ", file_graph[x[1]])

    def execute(self, forcerun):
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

        if "workflow" in locals() and self.workflow.persistence:
            self.workflow.persistence.unlock()

        self.logger.cleanup()

        return success


def get_dependencies_graph(params: SmkParam):
    smk_executor = SmkExecutor(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=params.forceall,
    )
    smk_executor.init_workflow(cores=params.cores)
    smk_executor.init_dag()
    smk_executor.init_graph()
    return smk_executor


def snakemake_execute(params: SmkParam):
    smk_executor = get_dependencies_graph(params)
    success = smk_executor.execute(params.forcerun)
    print("success: ", success)


if __name__ == '__main__':
    get_dependencies_graph()
