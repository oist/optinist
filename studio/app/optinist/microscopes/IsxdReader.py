import sys

try:
    import isx
except ModuleNotFoundError:
    pass

from studio.app.optinist.microscopes.MicroscopeDataReaderBase import (
    MicroscopeDataReaderBase,
    OMEDataModel,
)


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

    def _load_file(self, data_file_path: str) -> object:
        handle = isx.Movie.read(data_file_path)
        return (handle,)

    def _build_original_metadata(self, data_name: str) -> dict:
        movie: isx.Movie = None
        (movie,) = self.resource_handles

        spacing: isx.Spacing = movie.spacing
        timing: isx.Timing = movie.timing

        # Get acquisition start time
        # Note: Individual support for `isx v1.0.3`
        if not hasattr(timing, "start"):
            from datetime import datetime, timezone

            # Get acquisition start time from metadata
            timing_info_start = (
                movie.footer.get("timingInfo", {})
                .get("start", {})
                .get("secsSinceEpoch", {})
            )
            timing_info_start_num = int(timing_info_start.get("num"))
            timing_info_start_den = int(timing_info_start.get("den"))
            timing_start_unixtime = int(timing_info_start_num / timing_info_start_den)

            timing_start_datetime = datetime.fromtimestamp(
                timing_start_unixtime, timezone.utc
            )
            start_datatime = timing_start_datetime.strftime("%Y-%m-%d %H:%M:%S")

            del (
                timing_info_start,
                timing_info_start_num,
                timing_info_start_den,
                timing_start_unixtime,
                timing_start_datetime,
            )
        else:
            start_datatime = timing.start.to_datetime().strftime("%Y-%m-%d %H:%M:%S")

        original_metadata = {
            "data_name": data_name,
            "movie": {
                "data_type": movie.data_type.__name__,
            },
            "spacing": {
                "width": spacing.num_pixels[0],
                "height": spacing.num_pixels[1],
            },
            "timing": {
                "start": start_datatime,
                "period_msec": round(timing.period.secs_float * 1000),
                "num_samples": timing.num_samples,
                "dropped": (
                    timing.dropped if hasattr(timing, "dropped") else None
                ),  # Note: No attr at `isx v1.0.3`
                "cropped": (
                    timing.cropped if hasattr(timing, "cropped") else None
                ),  # Note: No attr at `isx v1.0.3`
                "blank": (
                    timing.blank if hasattr(timing, "blank") else None
                ),  # Note: No attr at `isx v1.0.3`
            },
        }

        return original_metadata

    def _build_ome_metadata(self, original_metadata: dict) -> OMEDataModel:
        """
        Note: Inscopix format is not supported in OME Bio-Formats
        """

        movie = original_metadata["movie"]
        spacing = original_metadata["spacing"]
        timing = original_metadata["timing"]

        fps = round(1000 / timing["period_msec"], 2)

        omeData = OMEDataModel(
            image_name=original_metadata["data_name"],
            size_x=spacing["width"],
            size_y=spacing["height"],
            size_t=timing["num_samples"],
            size_z=0,  # Note: currently unsettled
            size_c=0,  # Note: currently unsettled
            depth=OMEDataModel.get_depth_from_pixel_type(movie["data_type"]),
            significant_bits=0,  # Note: currently unsettled
            acquisition_date=timing["start"],
            objective_model=None,  # Note: currently unsettled
            fps=fps,
        )

        return omeData

    def _build_lab_specific_metadata(self, original_metadata: dict) -> dict:
        # Note: Not currently supported
        return None

    def _release_resources(self) -> None:
        # Note: in inscopix sdk, there is no library (ddl) release process.
        pass  # do nothing.

    def _get_image_stacks(self) -> list:
        movie: isx.Movie = None
        (movie,) = self.resource_handles

        # get frame data's np.ndarray.dtype
        # Note: Individual support for `isx v1.0.3`
        pixel_np_dtype = self.ome_metadata.pixel_np_dtype

        image_frames = [
            movie.get_frame_data(i).astype(pixel_np_dtype)
            for i in range(movie.timing.num_samples)
        ]

        return image_frames
