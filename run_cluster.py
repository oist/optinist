import argparse
import os
import shutil

from snakemake import snakemake

from studio.app.common.core.utils.filepath_creater import join_filepath
from studio.app.dir_path import DIRPATH


def main(args):
    # copy config file
    if args.config is not None:
        shutil.copyfile(
            args.config,
            join_filepath([DIRPATH.STUDIO_DIR, DIRPATH.SNAKEMAKE_CONFIG_YML]),
        )

    snakemake(
        DIRPATH.SNAKEMAKE_FILEPATH,
        forceall=args.forceall,
        cores=args.cores,
        use_conda=args.use_conda,
        workdir=f"{os.path.dirname(DIRPATH.STUDIO_DIR)}",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="optinist")
    parser.add_argument("--cores", type=int, default=2)
    parser.add_argument("--forceall", action="store_true")
    parser.add_argument("--use_conda", action="store_true")
    parser.add_argument("--config", type=str, default=None)

    main(parser.parse_args())
