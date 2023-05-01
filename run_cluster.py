import os
import argparse
import shutil
from snakemake import snakemake

from optinist.api.dir_path import DIRPATH
from optinist.api.utils.filepath_creater import join_filepath

def main(args):
    # copy config file
    if args.config is not None:
        shutil.copyfile(
            args.config,
            join_filepath([DIRPATH.ROOT_DIR, DIRPATH.SNAKEMAKE_CONFIG_YML])
        )

    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=args.forceall,
        cores=args.cores,
        use_conda=args.use_conda,
        workdir=f"{os.path.dirname(DIRPATH.ROOT_DIR)}",
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="optinist")
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument("--forceall", action="store_true")
    parser.add_argument("--use_conda", action="store_true")
    parser.add_argument("--config", type=str, default=None)

    main(parser.parse_args())
