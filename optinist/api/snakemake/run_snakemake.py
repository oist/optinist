from snakemake import snakemake

from optinist.api.dir_path import DIRPATH


def run_snakemake(params):
    # run snakemake
    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=params.forceall,
        cores=params.cores,
        use_conda=params.use_conda,
    )


# import os
# from snakemake.exceptions import print_exception
# from snakemake.logging import setup_logger, logger
# from snakemake.common import Mode
# from snakemake.workflow import Workflow
# from snakemake.dag import DAG
# from snakemake.persistence import Persistence

# from optinist.api.dir_path import DIRPATH


# def snakemake(
#     snakefile,
#     forcetargets=False,
#     forceall=False,
#     lock=False,
#     targets=None,
#     dryrun=False,
#     print_compilation=False,
#     keep_logger=False,
# ):
#     logger.setup_logfile()
#     workflow = Workflow(
#         snakefile=snakefile,
#     )
#     success = True

#     workflow.include(
#         snakefile,
#         overwrite_default_target=True,
#         print_compilation=print_compilation,
#     )
#     workflow.check()

#     success = workflow.execute(
#         targets=targets,
#         dryrun=dryrun,
#         forcetargets=forcetargets,
#         forceall=forceall,
#         nolock=not lock,
#         updated_files=[],
#     )

#     if "workflow" in locals() and workflow.persistence:
#         workflow.persistence.unlock()
#     if not keep_logger:
#         logger.cleanup()


# if __name__ == '__main__':
#     snakemake(
#         DIRPATH.SNAKEMAKE_FILEPATH,
#         forceall=True,
#     )
