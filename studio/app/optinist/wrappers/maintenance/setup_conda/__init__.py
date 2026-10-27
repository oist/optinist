from studio.app.optinist.wrappers.maintenance.setup_conda.setup_conda_caiman import (
    setup_conda_caiman,
)
from studio.app.optinist.wrappers.maintenance.setup_conda.setup_conda_custom import (
    setup_conda_custom,
)
from studio.app.optinist.wrappers.maintenance.setup_conda.setup_conda_lccd import (
    setup_conda_lccd,
)
from studio.app.optinist.wrappers.maintenance.setup_conda.setup_conda_microscope import (  # noqa: E501
    setup_conda_microscope,
)
from studio.app.optinist.wrappers.maintenance.setup_conda.setup_conda_optinist import (
    setup_conda_optinist,
)
from studio.app.optinist.wrappers.maintenance.setup_conda.setup_conda_suite2p import (
    setup_conda_suite2p,
)

setup_conda_wrapper_dict = {
    "setup_conda_caiman": {
        "function": setup_conda_caiman,
        "conda_name": "caiman",
    },
    "setup_conda_suite2p": {
        "function": setup_conda_suite2p,
        "conda_name": "suite2p",
    },
    "setup_conda_lccd": {
        "function": setup_conda_lccd,
        "conda_name": "lccd",
    },
    "setup_conda_microscope": {
        "function": setup_conda_microscope,
        "conda_name": "microscope",
    },
    "setup_conda_optinist": {
        "function": setup_conda_optinist,
        "conda_name": "optinist",
    },
    "setup_conda_custom": {
        "function": setup_conda_custom,
        "conda_name": "custom",
    },
}
