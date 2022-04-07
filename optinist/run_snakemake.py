from snakemake import snakemake
from snakemake.workflow import Workflow
from snakemake.rules import Rule
from snakemake.dag import DAG
from snakemake.persistence import Persistence

from cui_api.const import OPTINIST_DIR
from cui_api.utils import join_file_path

if __name__ == '__main__':
    snakefile = join_file_path([OPTINIST_DIR, 'Snakefile'])
    workflow = Workflow(
        snakefile=snakefile,
        
    )
    workflow.include(
        snakefile,
        overwrite_default_target=True,
    )
    workflow.check()
    # success = workflow.execute(
    #     forceall=True,
    #     # unlock=True,
    # )
    workflow.check_localrules()
    dag = DAG(
        workflow=workflow,
        rules=workflow.rules,
        forceall=True,
    )
    
    workflow.persistence = Persistence(
        dag=dag,
    )

    dag.init()
    dag.update_checkpoint_dependencies()
    dag.check_dynamic()
    workflow.persistence.lock()
    # # snakemake(
    # #     join_file_path([OPTINIST_DIR, 'Snakefile']),
    # # )
    # import pdb; pdb.set_trace()
