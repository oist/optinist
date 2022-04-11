import os
from snakemake.exceptions import print_exception
from snakemake.logging import setup_logger, logger
from snakemake.common import Mode
from snakemake.workflow import Workflow
from snakemake.dag import DAG
from snakemake.persistence import Persistence

from optinist.api.dir_path import DIRPATH


def snakemake(
    snakefile,
    forcetargets=False,
    forceall=False,
    cores=1,
    lock=True,
    targets=None,
    dryrun=False,
    printreason=False,
    printshellcmds=False,
    debug_dag=False,
    nocolor=False,
    quiet=False,
    print_compilation=False,
    keep_target_files=False,
    log_handler=[],
    keep_logger=False,
    verbose=False,
    mode=Mode.default,
    show_failed_logs=False,
):

    if not keep_logger:
        setup_logger(
            handler=log_handler,
            quiet=quiet,
            printreason=printreason,
            printshellcmds=printshellcmds,
            debug_dag=debug_dag,
            nocolor=nocolor,
            stdout=False,
            debug=verbose,
            mode=mode,
            show_failed_logs=show_failed_logs,
        )

    snakefile = os.path.abspath(snakefile)

    logger.setup_logfile()

    try:
        workflow = Workflow(
            snakefile=snakefile,
            cores=cores,
        )
        success = True

        workflow.include(
            snakefile,
            overwrite_default_target=True,
            print_compilation=print_compilation,
        )
        workflow.check()

        targetrules = map(
            workflow._rules.__getitem__,
            filter(workflow.is_rule, ["all"])
        )
        import pdb; pdb.set_trace()

        dag = DAG(
            workflow,
            workflow.rules,
            targetfiles=[],
            targetrules=targetrules,#set(workflow.rules(["all"])),
            forceall=True,
        )

        workflow.persistence = Persistence(
            nolock=True,
            dag=dag,
            conda_prefix=workflow.conda_prefix,
            singularity_prefix=workflow.singularity_prefix,
            shadow_prefix=workflow.shadow_prefix,
        )
        dag.init()

        graph = {}
        for i, job in enumerate(dag.jobs):
            graph[job.rule] = [
                dep.rule for dep in dag.dependencies[job]
            ]

        # node ids
        ids = {node: i for i, node in enumerate(graph)}

        file_graph = {}
        for i, job in enumerate(dag.jobs):
            file_graph[ids[job.rule]] = [files for files in dag.dependencies[job].values()]

        # calculate edges
        edges = [
            [ids[dep], ids[node]]
            for node, deps in graph.items()
            for dep in deps
        ]
        print(edges)

        print([[file_graph[x[0]], file_graph[x[1]]] for x in edges])

        success = workflow.execute(
            targets=targets,
            dryrun=True,
            forcetargets=forcetargets,
            forceall=forceall,
            keep_target_files=keep_target_files,
            updated_files=[],
            unlock=not lock,
            # printrulegraph=True,
        )
    except BrokenPipeError:
        success = False
    except (Exception, BaseException) as ex:
        if "workflow" in locals():
            print_exception(ex, workflow.linemaps)
        else:
            print_exception(ex, dict())
        success = False

    if "workflow" in locals() and workflow.persistence:
        workflow.persistence.unlock()
    if not keep_logger:
        logger.cleanup()
    return success


if __name__ == '__main__':
    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=True,
        # unlock=True,
    )
