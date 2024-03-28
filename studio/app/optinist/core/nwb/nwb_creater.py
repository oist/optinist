import os
import shutil
from datetime import datetime

from dateutil.tz import tzlocal
from pynwb import NWBHDF5IO, NWBFile
from pynwb.ophys import (
    CorrectedImageStack,
    Fluorescence,
    ImageSegmentation,
    ImageSeries,
    MotionCorrection,
    OpticalChannel,
    RoiResponseSeries,
    TimeSeries,
    TwoPhotonSeries,
)

from studio.app.optinist.core.nwb.nwb import NWBDATASET
from studio.app.optinist.core.nwb.optinist_data import PostProcess


class NWBCreater:
    @classmethod
    def acquisition(cls, config):
        nwbfile = NWBFile(
            session_description=config["session_description"],
            identifier=config["identifier"],
            experiment_description=config["experiment_description"],
            session_start_time=datetime.now(tzlocal()),
        )

        # 顕微鏡情報を登録
        device = nwbfile.create_device(
            name=config["device"]["name"],
            description=config["device"]["description"],
            manufacturer=config["device"]["manufacturer"],
        )

        # 光チャネルを登録
        optical_channel = OpticalChannel(
            name=config["optical_channel"]["name"],
            description=config["optical_channel"]["description"],
            emission_lambda=float(config["optical_channel"]["emission_lambda"]),
        )

        # imaging planeを追加
        imaging_plane = nwbfile.create_imaging_plane(
            name=config["imaging_plane"]["name"],
            description=config["imaging_plane"]["description"],
            optical_channel=optical_channel,  # 光チャネル
            device=device,  # 電極デバイス
            imaging_rate=float(config["imaging_plane"]["imaging_rate"]),  # 画像の比率Hz
            excitation_lambda=float(
                config["imaging_plane"]["excitation_lambda"]
            ),  # 励起（れいき）波長
            indicator=config["imaging_plane"]["indicator"],  # カルシウムインジケーター
            location=config["imaging_plane"]["location"],
        )

        # using internal data. this data will be stored inside the NWB file
        if (
            NWBDATASET.IMAGE_SERIES in config
            and "external_file" in config[NWBDATASET.IMAGE_SERIES]
        ):
            external_file = config[NWBDATASET.IMAGE_SERIES]["external_file"]
            image_path = external_file.path

            starting_frames = (
                config[NWBDATASET.IMAGE_SERIES]["starting_frame"]
                if "starting_frame" in config[NWBDATASET.IMAGE_SERIES]
                else None
            )
            save_raw_image_to_nwb = config[NWBDATASET.IMAGE_SERIES][
                "save_raw_image_to_nwb"
            ]

            if isinstance(image_path, list) and len(image_path) > 1:
                if starting_frames == 0 or starting_frames == [0]:
                    starting_frames = [0 for _ in range(len(image_path))]
                elif isinstance(starting_frames, str):
                    starting_frames = starting_frames.split(",")
                    starting_frames = list(map(int, starting_frames))
            elif isinstance(image_path, str):
                image_path = [image_path]

            image_series = TwoPhotonSeries(
                name="TwoPhotonSeries",
                starting_frame=starting_frames,
                external_file=image_path if not save_raw_image_to_nwb else None,
                imaging_plane=imaging_plane,
                starting_time=float(config[NWBDATASET.IMAGE_SERIES]["starting_time"]),
                rate=1.0,
                unit="normalized amplitude",
                data=external_file.data if save_raw_image_to_nwb else None,
            )
            nwbfile.add_acquisition(image_series)

        nwbfile.create_processing_module(
            name="ophys", description="optical physiology processed data"
        )

        nwbfile.create_processing_module(name="optinist", description="description")
        cls.ophys(nwbfile)

        return nwbfile

    @classmethod
    def ophys(cls, nwbfile):
        img_seg = ImageSegmentation()
        nwbfile.processing["ophys"].add(img_seg)

        if "TwoPhotonSeries" in nwbfile.acquisition:
            reference_images = nwbfile.acquisition["TwoPhotonSeries"]

            img_seg.create_plane_segmentation(
                name="PlaneSegmentation",
                description="output",
                imaging_plane=nwbfile.imaging_planes["ImagingPlane"],
                reference_images=reference_images,
            )

        return nwbfile

    @classmethod
    def add_plane_segmentation(cls, nwbfile, function_id):
        image_seg = nwbfile.processing["ophys"].data_interfaces["ImageSegmentation"]
        if "TwoPhotonSeries" in nwbfile.acquisition:
            reference_images = nwbfile.acquisition["TwoPhotonSeries"]

            try:
                image_seg.plane_segmentations.pop(function_id)
                nwbfile.processing["ophys"].data_interfaces.pop(function_id)
            except KeyError:
                pass

            image_seg.create_plane_segmentation(
                name=function_id,
                description="output",
                imaging_plane=nwbfile.imaging_planes["ImagingPlane"],
                reference_images=reference_images,
            )

        return nwbfile

    @classmethod
    def motion_correction(cls, nwbfile, function_id, mc_data, xy_trans_data):
        # image_data = mc_data.data
        image_path = mc_data.path
        corrected = ImageSeries(
            name="corrected",  # this must be named "corrected"
            # data=image_data,
            external_file=image_path,
            unit="na",
            format="external",
            starting_time=0.0,
            rate=1.0,
        )

        xy_translation = TimeSeries(
            name="xy_translation",
            data=xy_trans_data,
            unit="pixels",
            starting_time=0.0,
            rate=1.0,
        )

        corrected_image_stack = CorrectedImageStack(
            corrected=corrected,
            original=nwbfile.acquisition["TwoPhotonSeries"],
            xy_translation=xy_translation,
        )

        motion_correction = MotionCorrection(
            corrected_image_stacks=corrected_image_stack
        )

        try:
            nwbfile.processing.pop(function_id)
        except KeyError:
            pass

        function_process = nwbfile.create_processing_module(
            name=function_id, description="processed by " + function_id
        )
        function_process.add(motion_correction)

        return nwbfile

    @classmethod
    def roi(cls, nwbfile, function_id, roi_list):
        nwbfile = cls.add_plane_segmentation(nwbfile, function_id)
        image_seg = nwbfile.processing["ophys"].data_interfaces["ImageSegmentation"]
        plane_seg = image_seg.plane_segmentations[function_id]

        if roi_list:
            for col in roi_list[0]:
                if col != "pixel_mask" and col not in plane_seg.colnames:
                    plane_seg.add_column(col, f"{col} list")

        for col in roi_list:
            plane_seg.add_roi(**col)

        return nwbfile

    @classmethod
    def column(cls, nwbfile, function_id, name, description, data):
        image_seg = nwbfile.processing["ophys"].data_interfaces["ImageSegmentation"]
        plane_seg = image_seg.plane_segmentations[function_id]
        plane_seg.add_column(name, description, data)

        return nwbfile

    @classmethod
    def fluorescence(cls, nwbfile, function_id, roi_list):
        image_seg = nwbfile.processing["ophys"].data_interfaces["ImageSegmentation"]
        plane_seg = image_seg.plane_segmentations[function_id]
        fluo = Fluorescence(name=function_id)
        for key in roi_list.keys():
            roi = roi_list[key]
            region_roi = plane_seg.create_roi_table_region(
                roi["table_name"], region=roi["region"]
            )

            roi_resp_series = RoiResponseSeries(
                name=roi["name"],
                data=roi["data"],
                rois=region_roi,
                unit=roi["unit"],
                timestamps=roi.get("timestamps"),
                rate=float(roi.get("rate", 0.0)),
            )
            fluo.add_roi_response_series(roi_resp_series)

        nwbfile.processing["ophys"].add(fluo)

        return nwbfile

    @classmethod
    def timeseries(cls, nwbfile, key, value):
        timeseries_data = TimeSeries(
            name=key,
            data=value.data,
            unit="second",
            starting_time=0.0,
            rate=1.0,
        )

        nwbfile.processing["ophys"].add(timeseries_data)

        return nwbfile

    @classmethod
    def behavior(cls, nwbfile, key, value):
        timeseries_data = TimeSeries(
            name=key,
            data=value.data,
            unit="second",
            starting_time=0.0,
            rate=1.0,
        )

        nwbfile.processing["optinist"].add(timeseries_data)

        return nwbfile

    @classmethod
    def postprocess(cls, nwbfile, function_id, data):
        for key, value in data.items():
            process_name = f"{function_id}_{key}"
            postprocess = PostProcess(name=process_name, data=value)

            try:
                nwbfile.processing["optinist"].add_container(postprocess)
            except ValueError:
                nwbfile.processing["optinist"].data_interfaces.pop(process_name)
                nwbfile.processing["optinist"].add_container(postprocess)

        return nwbfile

    @classmethod
    def reaqcuisition(cls, nwbfile):
        new_nwbfile = NWBFile(
            session_description=nwbfile.session_description,
            identifier=nwbfile.identifier,
            experiment_description=nwbfile.experiment_description,
            session_start_time=nwbfile.session_start_time,
        )

        devices = []
        for key in nwbfile.devices.keys():
            device = new_nwbfile.create_device(
                name=key,
                description=nwbfile.devices[key].description,
                manufacturer=nwbfile.devices[key].manufacturer,
            )
            devices.append(device)

        old_optical_channel = nwbfile.imaging_planes["ImagingPlane"].optical_channel[0]
        optical_channels = []
        for old_optical_channel in nwbfile.imaging_planes[
            "ImagingPlane"
        ].optical_channel:
            optical_channel = OpticalChannel(
                name=old_optical_channel.name,
                description=old_optical_channel.description,
                emission_lambda=old_optical_channel.emission_lambda,
            )
            optical_channels.append(optical_channel)

        imaging_planes = []
        for key in nwbfile.imaging_planes.keys():
            imaging_plane = new_nwbfile.create_imaging_plane(
                name=nwbfile.imaging_planes[key].name,
                description=nwbfile.imaging_planes[key].description,
                optical_channel=optical_channels,  # 光チャネル
                device=devices[0],  # 電極デバイス
                imaging_rate=float(
                    nwbfile.imaging_planes["ImagingPlane"].imaging_rate
                ),  # 画像の比率Hz
                excitation_lambda=float(
                    nwbfile.imaging_planes["ImagingPlane"].excitation_lambda
                ),  # 励起（れいき）波長
                indicator=nwbfile.imaging_planes[
                    "ImagingPlane"
                ].indicator,  # カルシウムインディケーター
                location=nwbfile.imaging_planes["ImagingPlane"].location,
            )
            imaging_planes.append(imaging_plane)

        image_series = TwoPhotonSeries(
            name="TwoPhotonSeries",
            starting_frame=nwbfile.acquisition["TwoPhotonSeries"].starting_frame,
            external_file=nwbfile.acquisition["TwoPhotonSeries"].external_file,
            imaging_plane=imaging_planes[0],
            starting_time=nwbfile.acquisition["TwoPhotonSeries"].starting_time,
            rate=nwbfile.acquisition["TwoPhotonSeries"].rate,
            unit=nwbfile.acquisition["TwoPhotonSeries"].unit,
        )

        new_nwbfile.add_acquisition(image_series)

        new_nwbfile.create_processing_module(
            name="ophys", description="optical physiology processed data"
        )
        new_nwbfile.create_processing_module(name="optinist", description="description")
        cls.ophys(new_nwbfile)
        return new_nwbfile


def set_nwbconfig(nwbfile, config):
    if NWBDATASET.POSTPROCESS in config:
        for function_key in config[NWBDATASET.POSTPROCESS]:
            NWBCreater.postprocess(
                nwbfile, function_key, config[NWBDATASET.POSTPROCESS][function_key]
            )

    if NWBDATASET.TIMESERIES in config:
        for key, value in config[NWBDATASET.TIMESERIES].items():
            NWBCreater.timeseries(nwbfile, key, value)

    if NWBDATASET.BEHAVIOR in config:
        for key, value in config[NWBDATASET.BEHAVIOR].items():
            NWBCreater.behavior(nwbfile, key, value)

    if NWBDATASET.MOTION_CORRECTION in config:
        for function_key in config[NWBDATASET.MOTION_CORRECTION]:
            nwbfile = NWBCreater.motion_correction(
                nwbfile,
                function_key,
                **config[NWBDATASET.MOTION_CORRECTION][function_key],
            )

    if NWBDATASET.ROI in config:
        for function_key in config[NWBDATASET.ROI]:
            nwbfile = NWBCreater.roi(
                nwbfile, function_key, config[NWBDATASET.ROI][function_key]
            )

    if NWBDATASET.COLUMN in config:
        for function_key in config[NWBDATASET.COLUMN]:
            nwbfile = NWBCreater.column(
                nwbfile, function_key, **config[NWBDATASET.COLUMN][function_key]
            )

    if NWBDATASET.FLUORESCENCE in config:
        for function_key in config[NWBDATASET.FLUORESCENCE]:
            nwbfile = NWBCreater.fluorescence(
                nwbfile,
                function_key,
                config[NWBDATASET.FLUORESCENCE][function_key],
            )

    return nwbfile


def save_nwb(save_path, input_config, config):
    nwbfile = NWBCreater.acquisition(input_config)

    nwbfile = set_nwbconfig(nwbfile, config)

    with NWBHDF5IO(save_path, "w") as f:
        f.write(nwbfile)


def overwrite_nwbfile(save_path, config):
    tmp_save_path = os.path.join(
        os.path.dirname(save_path),
        "tmp_" + os.path.basename(save_path),
    )
    with NWBHDF5IO(save_path, "r") as src_io:
        old_nwbfile = src_io.read()
        nwbfile = set_nwbconfig(old_nwbfile, config)
        nwbfile.set_modified()
        with NWBHDF5IO(tmp_save_path, mode="w") as io:
            io.export(src_io=src_io, nwbfile=nwbfile)
    shutil.copyfile(tmp_save_path, save_path)
    os.remove(tmp_save_path)


def overwrite_nwb(config, save_path, nwb_file_name):
    # バックアップファイルを作成
    nwb_path = os.path.join(save_path, nwb_file_name)
    tmp_nwb_path = os.path.join(save_path, "tmp_" + nwb_file_name)

    # NWBファイルの読み込み
    with NWBHDF5IO(nwb_path, "r") as io:
        nwbfile = io.read()
        # acquisition を元ファイルから作成する
        new_nwbfile = NWBCreater.reaqcuisition(nwbfile)
        new_nwbfile = set_nwbconfig(new_nwbfile, config)

        with NWBHDF5IO(tmp_nwb_path, "w") as io:
            io.write(new_nwbfile)
    shutil.copyfile(tmp_nwb_path, nwb_path)
    os.remove(tmp_nwb_path)


def merge_nwbfile(old_nwbfile, new_nwbfile):
    for pattern in [
        NWBDATASET.POSTPROCESS,
        NWBDATASET.TIMESERIES,
        NWBDATASET.MOTION_CORRECTION,
        NWBDATASET.ROI,
        NWBDATASET.COLUMN,
        NWBDATASET.FLUORESCENCE,
        NWBDATASET.BEHAVIOR,
        NWBDATASET.IMAGE_SERIES,
    ]:
        if pattern in old_nwbfile and pattern in new_nwbfile:
            for function_id in new_nwbfile[pattern]:
                if function_id in old_nwbfile[pattern]:
                    old_nwbfile[pattern][function_id].update(
                        new_nwbfile[pattern][function_id]
                    )
                else:
                    old_nwbfile[pattern][function_id] = new_nwbfile[pattern][
                        function_id
                    ]
        elif pattern in new_nwbfile:
            old_nwbfile[pattern] = new_nwbfile[pattern]

    return old_nwbfile
