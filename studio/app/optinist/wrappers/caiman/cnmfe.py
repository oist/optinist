from studio.app.common.dataclass import ImageData
from studio.app.optinist.dataclass import FluoData, IscellData
from studio.app.optinist.wrappers.caiman.cnmf import caiman_cnmf


def caiman_cnmfe(
    images: ImageData, output_dir: str, params: dict = None, **kwargs
) -> dict(fluorescence=FluoData, iscell=IscellData):
    cnmfe_fixed_params = {
        "center_psf": True,
        "method_init": "corr_pnr",  # use this for 1 photon
        "only_init": True,  # set it to True to run CNMF-E
        "normalize_init": False,
    }
    params["cnmfe_fixed_params"] = cnmfe_fixed_params

    return caiman_cnmf(images=images, output_dir=output_dir, params=params, **kwargs)
