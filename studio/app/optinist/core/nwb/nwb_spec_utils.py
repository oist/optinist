import os
import time
from typing import List, Union

from filelock import FileLock
from pynwb.spec import NWBGroupSpec, NWBNamespaceBuilder

from studio.app.common.core.utils.filelock_handler import FileLockUtils
from studio.app.dir_path import DIRPATH
from studio.app.version import Version

NWB_SPEC_FILE_EXPORT_DIR = os.path.join(os.path.dirname(__file__), "specs")
root_dir = DIRPATH.ROOT_DIR


def get_namespace_file_path(ns_name: str) -> str:
    ns_path = os.path.join(NWB_SPEC_FILE_EXPORT_DIR, f"{ns_name}.namespace.yaml")
    return ns_path


def export_spec_files(
    ns_name: str,
    group_spec: Union[NWBGroupSpec, List[NWBGroupSpec]],
    file_cache_interval: int = 60,
):
    """
    Generation of NWB namespace files
      (*.namespace.yaml, *.extensions.yaml)
    """

    ns_path = get_namespace_file_path(ns_name)

    ext_source = (
        f"{ns_name}.extensions.yaml"  # This path must be specified in basename.
    )

    # Note:
    # Considering calls from multi-process exclusive processing
    # is performed (using FileLock)
    lock_path = FileLockUtils.get_lockfile_path(ns_path)
    with FileLock(lock_path, timeout=10):
        flle_update_elapsed_time = (
            (time.time() - os.path.getmtime(ns_path)) if os.path.exists(ns_path) else 0
        )

        if (not os.path.exists(ns_path)) or (
            flle_update_elapsed_time > file_cache_interval
        ):
            current_version = Version.APP_VERSION
            ns_builder = NWBNamespaceBuilder(
                f"{ns_name} extensions", ns_name, version=current_version
            )

            if isinstance(group_spec, list):
                group_specs = group_spec
            else:
                group_specs = [group_spec]

            for spec in group_specs:
                ns_builder.add_spec(ext_source, spec)

            previous_dir = os.getcwd()

            try:
                # ATTENTION:
                #  There seems to be a restriction that the extensions file
                #  is automatically created in the current directory.
                #  Therefore, before generating the extensions file,
                # move the current directory to the output path.
                # *The extensions file is generated in NWBNamespaceBuilder.export.
                os.chdir(NWB_SPEC_FILE_EXPORT_DIR)

                # Export nwb spec files
                ns_builder.export(ns_path)
            finally:
                # Return current directory
                os.chdir(previous_dir)
