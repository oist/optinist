from studio.app.common.core.experiment.experiment import ExptOutputPathIds
from studio.app.common.core.logger import AppLogger
from studio.app.common.dataclass import ImageData
from studio.app.optinist.dataclass import Suite2pData

logger = AppLogger.get_logger()


def suite2p_registration(
    ops: Suite2pData, output_dir: str, params: dict = None, **kwargs
) -> dict(ops=Suite2pData, mc_images=ImageData):
    from suite2p import default_ops, io, registration

    function_id = ExptOutputPathIds(output_dir).function_id
    logger.info("start suite2p registration: %s", function_id)

    ops = ops.data
    refImg = ops["meanImg"]

    # REGISTRATION
    if len(refImg.shape) == 3:
        refImg = refImg[0]

    ops = {**default_ops(), **ops, **params}

    # register binary
    ops = registration.register_binary(ops, refImg=refImg)

    # compute metrics for registration
    if ops.get("do_regmetrics", True) and ops["nframes"] >= 1500:
        ops = registration.get_pc_metrics(ops)

    mv = io.BinaryFile(
        Lx=ops["Lx"], Ly=ops["Ly"], read_filename=ops["reg_file"]
    ).data.copy()

    info = {
        "refImg": ImageData(ops["refImg"], output_dir=output_dir, file_name="refImg"),
        "meanImgE": ImageData(
            ops["meanImgE"], output_dir=output_dir, file_name="meanImgE"
        ),
        "mc_images": ImageData(mv, output_dir=output_dir, file_name="mc_images"),
        "ops": Suite2pData(ops, file_name="ops"),
    }

    return info
