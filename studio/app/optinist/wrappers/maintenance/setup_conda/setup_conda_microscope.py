from studio.app.common.core.logger import AppLogger
from studio.app.common.dataclass import ImageData

logger = AppLogger.get_logger()


def setup_conda_microscope(
    image: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict():
    """
    Note:
        This wrapper function is a mock function prepared to instruct snakemake
        to create a conda env, and there is no specific processing within this function.
    """

    logger.info("setup conda env done.")

    info = {"image": image}  # set dummy data

    return info
