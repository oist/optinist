from pynwb import NWBFile, NWBHDF5IO
from pynwb.ophys import (
    OpticalChannel, TwoPhotonSeries, ImageSegmentation,
    RoiResponseSeries, Fluorescence, ImageSeries, TimeSeries,
    CorrectedImageStack, MotionCorrection,
)

from datetime import datetime
from dateutil.tz import tzlocal

from optinist.api.nwb.optinist_data import PostProcess
from optinist.api.nwb.nwb import NWBDATASET


class NWBCreater:
    @classmethod
    def acquisition(cls, config):
        nwbfile = NWBFile(
            session_description=config['session_description'],
            identifier=config['identifier'],
            experiment_description=config['experiment_description'],
            session_start_time=datetime.now(tzlocal()),
        )

        # 顕微鏡情報を登録
        device = nwbfile.create_device(
            name=config['device']['name'], 
            description=config['device']['description'],
            manufacturer=config['device']['manufacturer']
        )

        # 光チャネルを登録
        optical_channel = OpticalChannel(
            name=config['optical_channel']['name'], 
            description=config['optical_channel']['description'], 
            emission_lambda=float(config['optical_channel']['emission_lambda'])
        )

        # imaging planeを追加
        imaging_plane = nwbfile.create_imaging_plane(
            name=config['imaging_plane']['name'],
            description=config['imaging_plane']['description'],
            optical_channel=optical_channel,   # 光チャネル
            device=device,   # 電極デバイス
            imaging_rate=float(config['imaging_plane']['imaging_rate']),   # 画像の比率Hz
            excitation_lambda=float(config['imaging_plane']['excitation_lambda']), # 励起（れいき）波長
            indicator=config['imaging_plane']['indicator'],   # カルシウムインディケーター
            location=config['imaging_plane']['location'],
        )

        # using internal data. this data will be stored inside the NWB file
        if NWBDATASET.IMAGE_SERIES in config and 'external_file' in config[NWBDATASET.IMAGE_SERIES]:
            # image_data = config[NWBDATASET.IMAGE_SERIES]['external_file'].data
            image_path = config[NWBDATASET.IMAGE_SERIES]['external_file'].path
            starting_frames = config[NWBDATASET.IMAGE_SERIES]['starting_frame'] if 'starting_frame' in config[NWBDATASET.IMAGE_SERIES] else None

            if isinstance(image_path, list) and len(image_path) > 1:
                if starting_frames == 0 or starting_frames == [0]:
                    starting_frames = [0 for _ in range(len(image_path))]
                elif isinstance(starting_frames, str):
                    starting_frames = starting_frames.split(',')
                    starting_frames =  list(map(int, starting_frames))

            image_series = TwoPhotonSeries(
                name='TwoPhotonSeries',
                starting_frame = starting_frames,
                external_file=image_path,
                imaging_plane=imaging_plane,
                rate=1.0,
                unit='normalized amplitude'
            )
            nwbfile.add_acquisition(image_series)

        return nwbfile

    @classmethod
    def ophys(cls, nwbfile):
        img_seg = ImageSegmentation()
        nwbfile.processing['ophys'].add(img_seg)

        if 'TwoPhotonSeries' in nwbfile.acquisition:
            reference_images = nwbfile.acquisition['TwoPhotonSeries']

            img_seg.create_plane_segmentation(
                name='PlaneSegmentation',
                description='output',
                imaging_plane=nwbfile.imaging_planes['ImagingPlane'],
                reference_images=reference_images,
            )

        return nwbfile

    @classmethod
    def motion_correction(cls, nwbfile, mc_data, xy_trans_data):
        # image_data = mc_data.data
        image_path = mc_data.path
        corrected = ImageSeries(
            name='corrected',  # this must be named "corrected"
            # data=image_data,
            external_file=image_path,
            unit='na',
            format='external',
            starting_time=0.0,
            rate=1.0
        )

        xy_translation = TimeSeries(
            name='xy_translation',
            data=xy_trans_data,
            unit='pixels',
            starting_time=0.0,
            rate=1.0,
        )

        corrected_image_stack = CorrectedImageStack(
            corrected=corrected,
            original=nwbfile.acquisition['TwoPhotonSeries'],
            xy_translation=xy_translation,
        )

        motion_correction = MotionCorrection(
            corrected_image_stacks=corrected_image_stack
        )
        nwbfile.processing['ophys'].add(motion_correction)

        return nwbfile

    @classmethod
    def column(cls, nwbfile, name, discription, data):
        data_interfaces = nwbfile.processing['ophys'].data_interfaces
        plane_seg = data_interfaces['ImageSegmentation'].plane_segmentations['PlaneSegmentation']
        plane_seg.add_column(name, discription, data)

        return nwbfile

    @classmethod
    def roi(cls, nwbfile, roi_list):
        data_interfaces = nwbfile.processing['ophys'].data_interfaces
        plane_seg = data_interfaces['ImageSegmentation'].plane_segmentations['PlaneSegmentation']

        for col in roi_list[0]:
            if col != 'pixel_mask' and col not in plane_seg.colnames:
                plane_seg.add_column(col, f'{col} list')

        for col in roi_list:
            plane_seg.add_roi(**col)

        return nwbfile

    @classmethod
    def fluorescence(
        cls, nwbfile, table_name, region, name,
        data, unit, timestamps=None, rate=0.0
    ):

        data_interfaces = nwbfile.processing['ophys'].data_interfaces
        plane_seg = data_interfaces['ImageSegmentation'].plane_segmentations['PlaneSegmentation']

        region_roi = plane_seg.create_roi_table_region(
            table_name, region=region)

        roi_resp_series = RoiResponseSeries(
            name=name,
            data=data,
            rois=region_roi,
            unit=unit,
            timestamps=timestamps,
            rate=float(rate)
        )

        fluo = Fluorescence(
            name=name,
            roi_response_series=roi_resp_series
        )

        nwbfile.processing['ophys'].add(fluo)

        return nwbfile

    @classmethod
    def timeseries(cls, nwbfile, key, value):
        timeseries_data = TimeSeries(
            name=key,
            data=value.data,
            unit='second',
            starting_time=0.0,
            rate=1.0,
        )

        nwbfile.processing['ophys'].add(timeseries_data)

        return nwbfile

    @classmethod
    def behavior(cls, nwbfile, key, value):
        timeseries_data = TimeSeries(
            name=key,
            data=value.data,
            unit='second',
            starting_time=0.0,
            rate=1.0,
        )

        nwbfile.processing['optinist'].add(timeseries_data)

        return nwbfile

    @classmethod
    def postprocess(cls, nwbfile, key, value):
        postprocess = PostProcess(
            name=key,
            data=value,
        )

        nwbfile.processing['optinist'].add_container(postprocess)

        return nwbfile


def save_nwb(save_path, input_config, config):
    nwbfile = NWBCreater.acquisition(input_config)
    nwbfile.create_processing_module(
        name='ophys',
        description='optical physiology processed data'
    )
    NWBCreater.ophys(nwbfile)

    nwbfile.create_processing_module(
        name='optinist', 
        description='description'
    )

    if NWBDATASET.POSTPROCESS in config:
        for key, value in config[NWBDATASET.POSTPROCESS].items():
            NWBCreater.postprocess(nwbfile, key, value)

    if NWBDATASET.TIMESERIES in config:
        for key, value in config[NWBDATASET.TIMESERIES].items():
            NWBCreater.timeseries(nwbfile, key, value)

    if NWBDATASET.BEHAVIOR in config:
        for key, value in config[NWBDATASET.BEHAVIOR].items():
            NWBCreater.behavior(nwbfile, key, value)

    if NWBDATASET.MOTION_CORRECTION in config:
        for mc in config[NWBDATASET.MOTION_CORRECTION].values():
            nwbfile = NWBCreater.motion_correction(nwbfile, **mc)

    if NWBDATASET.ROI in config:
        for roi_list in config[NWBDATASET.ROI].values():
            NWBCreater.roi(nwbfile, roi_list)

    if NWBDATASET.COLUMN in config:
        for value in config[NWBDATASET.COLUMN].values():
            nwbfile = NWBCreater.column(nwbfile, **value)

    if NWBDATASET.FLUORESCENCE in config:
        for value in config[NWBDATASET.FLUORESCENCE].values():
            nwbfile = NWBCreater.fluorescence(nwbfile, **value)

    with NWBHDF5IO(save_path, 'w') as f:
        f.write(nwbfile)


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
            old_nwbfile[pattern].update(new_nwbfile[pattern])
        elif pattern in new_nwbfile:
            old_nwbfile[pattern] = new_nwbfile[pattern]

    return old_nwbfile