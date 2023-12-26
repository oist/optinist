import sys

import isx
from MicroscopeDataReaderBase import MicroscopeDataReaderBase, OMEDataModel


class IsxdReader(MicroscopeDataReaderBase):
    """Inscopix isxd data reader"""

    SDK_MODULE_NAME = "isx"

    @staticmethod
    def get_library_path() -> str:
        """Returns the path of the library (dll) file"""
        # Note: In inscopix SDK, library (ddl) files are not used directly.
        return None  # do nothing.

    @staticmethod
    def is_available() -> bool:
        """Determine if library is available"""
        return __class__.SDK_MODULE_NAME in sys.modules

    def _init_library(self):
        # Note: in inscopix sdk, there is no library (ddl)
        #         initialization process. (using pip module)
        pass  # do nothing.

    def _load_data_file(self, data_file_path: str) -> object:
        return isx.Movie.read(data_file_path)

    def _build_original_metadata(self, handle: object, data_name: str) -> dict:
        movie: isx.Movie = handle
        spacing: isx.Spacing = movie.spacing
        timing: isx.Timing = movie.timing

        original_metadata = {
            "data_name": data_name,
            "spacing": {
                "width": spacing.num_pixels[0],
                "height": spacing.num_pixels[1],
            },
            "timing": {
                "start": timing.start.to_datetime().strftime("%Y-%m-%d %H:%M:%S"),
                "period_msec": timing.period.to_msecs(),
                "num_samples": timing.num_samples,
                "dropped": timing.dropped,
                "cropped": timing.cropped,
                "blank": timing.blank,
            },
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        @link OME/NativeND2Reader
        """

        spacing = original_metadata["spacing"]
        timing = original_metadata["timing"]

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=spacing["width"],
            size_y=spacing["height"],
            size_t=timing["num_samples"],
            size_c=0,
            fps=(1000 / timing["period_msec"]),
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # Note: Not currently supported
        return None

    def _release_resources(self, handle: object) -> None:
        # Note: in inscopix sdk, there is no library (ddl) release process.
        pass  # do nothing.

    def get_images_stack(self) -> list:
        movie: isx.Movie = self.data_handle
        image_frames = []

        for i in range(movie.timing.num_samples):
            image_frames.append(movie.get_frame_data(i))

        return image_frames
