import os
from pynwb import NWBFile, NWBHDF5IO
from pynwb.ophys import (
    OpticalChannel, TwoPhotonSeries, ImageSegmentation,
    RoiResponseSeries, Fluorescence, ImageSeries, TimeSeries,
    CorrectedImageStack, MotionCorrection,
)

from pynwb.core import NWBDataInterface

from datetime import datetime
from dateutil.tz import tzlocal

from .optinist_data import PostProcess
from wrappers.nwb_wrapper.const import NWBDATASET


def nwb_add_acquisition(nwb_dict):
    nwbfile = NWBFile(
        session_description=nwb_dict['session_description'],
        identifier=nwb_dict['identifier'],
        experiment_description=nwb_dict['experiment_description'],
        session_start_time=datetime.now(tzlocal()),
    )

    # 顕微鏡情報を登録
    device = nwbfile.create_device(
        name=nwb_dict['device']['name'], 
        description=nwb_dict['device']['description'],
        manufacturer=nwb_dict['device']['manufacturer']
    )

    # 光チャネルを登録
    optical_channel = OpticalChannel(
        name=nwb_dict['optical_channel']['name'], 
        description=nwb_dict['optical_channel']['description'], 
        emission_lambda=nwb_dict['optical_channel']['emission_lambda']
    )

    # imaging planeを追加
    imaging_plane = nwbfile.create_imaging_plane(
        name=nwb_dict['imaging_plane']['name'],
        description=nwb_dict['imaging_plane']['description'],
        optical_channel=optical_channel,   # 光チャネル
        device=device,   # 電極デバイス
        imaging_rate=nwb_dict['imaging_plane']['imaging_rate'],   # 画像の比率Hz
        excitation_lambda=nwb_dict['imaging_plane']['excitation_lambda'], # 励起（れいき）波長
        indicator=nwb_dict['imaging_plane']['indicator'],   # カルシウムインディケーター
        location=nwb_dict['imaging_plane']['location'],
    )

    # using internal data. this data will be stored inside the NWB file
    if 'external_file' in nwb_dict['image_series'].keys():
        image_data = nwb_dict['image_series']['external_file'].data
        image_series = TwoPhotonSeries(
            name='TwoPhotonSeries',
            data=image_data,
            imaging_plane=imaging_plane,
            rate=1.0,
            unit='normalized amplitude'
        )
        nwbfile.add_acquisition(image_series)

    return nwbfile


def nwb_add_ophys(nwbfile):
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


def nwb_motion_correction(nwbfile, mc_data, xy_trans_data):
    corrected = ImageSeries(
        name='corrected',  # this must be named "corrected"
        data=mc_data.data,
        unit='na',
        format='raw',
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


def nwb_add_column(nwbfile, name, discription, data):
    data_interfaces = nwbfile.processing['ophys'].data_interfaces
    plane_seg = data_interfaces['ImageSegmentation'].plane_segmentations['PlaneSegmentation']
    plane_seg.add_column(name, discription, data)

    return nwbfile


def nwb_add_roi(nwbfile, roi_list):
    data_interfaces = nwbfile.processing['ophys'].data_interfaces
    plane_seg = data_interfaces['ImageSegmentation'].plane_segmentations['PlaneSegmentation']

    for col in roi_list[0].keys():
        if col != 'pixel_mask' and col not in plane_seg.colnames:
            plane_seg.add_column(col, f'{col} list')

    for col in roi_list:
        plane_seg.add_roi(**col)

    return nwbfile


def nwb_add_fluorescence(
        nwbfile, table_name, region, name, data, unit, 
        timestamps=None, rate=0.0):

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


def nwb_add_timeseries(nwbfile, key, value):

    timeseries_data = TimeSeries(
        name=key,
        data=value.data,
        unit='second',
        starting_time=0.0,
        rate=1.0,
    )

    nwbfile.processing['ophys'].add(timeseries_data)

    return nwbfile


def nwb_add_postprocess(nwbfile, key, value):
    postprocess = PostProcess(
        name=key,
        data=value,
    )

    nwbfile.processing['optinist'].add_container(postprocess)

    return nwbfile


def save_nwb(nwb_dict, save_path):
    nwbfile = nwb_add_acquisition(nwb_dict)
    nwbfile.create_processing_module(
        name='ophys',
        description='optical physiology processed data'
    )
    nwb_add_ophys(nwbfile)

    nwbfile.create_processing_module(
        name='optinist', 
        description='description'
    )

    if NWBDATASET.POSTPROCESS in nwb_dict:
        for key, value in nwb_dict[NWBDATASET.POSTPROCESS].items():
            nwb_add_postprocess(nwbfile, key, value)

    if NWBDATASET.TIMESERIES in nwb_dict:
        for key, value in nwb_dict[NWBDATASET.TIMESERIES].items():
            nwb_add_timeseries(nwbfile, key, value)

    if NWBDATASET.MOTION_CORRECTION in nwb_dict:
        for mc in nwb_dict[NWBDATASET.MOTION_CORRECTION].values():
            nwbfile = nwb_motion_correction(
                nwbfile, **mc)

    if NWBDATASET.ROI in nwb_dict:
        for roi_list in nwb_dict[NWBDATASET.ROI].values():
            nwb_add_roi(nwbfile, roi_list)

    if NWBDATASET.COLUMN in nwb_dict:
        for value in nwb_dict[NWBDATASET.COLUMN].values():
            nwbfile = nwb_add_column(nwbfile, **value)

    if NWBDATASET.FLUORESCENCE in nwb_dict:
        for value in nwb_dict[NWBDATASET.FLUORESCENCE].values():
            nwbfile = nwb_add_fluorescence(nwbfile, **value)

    with NWBHDF5IO(f'{save_path}.nwb', 'w') as f:
        f.write(nwbfile)
