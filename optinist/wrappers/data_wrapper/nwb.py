from pynwb import NWBFile
from pynwb.ophys import (
    OpticalChannel, TwoPhotonSeries, ImageSegmentation,
    RoiResponseSeries, Fluorescence, ImageSeries, TimeSeries,
    CorrectedImageStack, MotionCorrection
)
from datetime import datetime
from dateutil.tz import tzlocal


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
    image_series = TwoPhotonSeries(
        name='TwoPhotonSeries',
        data=nwb_dict['image_series']['external_file'],
        imaging_plane=imaging_plane,
        rate=1.0,
        unit='normalized amplitude'
    )

    nwbfile.add_acquisition(image_series)

    return nwbfile


def nwb_add_ophys(nwbfile):
    img_seg = ImageSegmentation()
    nwbfile.processing['ophys'].add(img_seg)

    img_seg.create_plane_segmentation(
        name='PlaneSegmentation',
        description='suite2p output',
        imaging_plane=nwbfile.imaging_planes['ImagingPlane'],
        reference_images=nwbfile.acquisition['TwoPhotonSeries']
    )

    return nwbfile


def nwb_motion_correction(nwbfile, mc_data, xy_trans_data):
    corrected = ImageSeries(
        name='corrected',  # this must be named "corrected"
        data=mc_data,
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


def nwb_add_ps_column(nwbfile, roi_list):
    ohys = nwbfile.processing['ophys']
    image_seg = ohys.data_interfaces['ImageSegmentation']
    plane_seg = image_seg.plane_segmentations['PlaneSegmentation']

    for col in roi_list[0].keys():
        plane_seg.add_column(col, f'in {col} list')

    for col in roi_list:
        plane_seg.add_roi(**col)

    return nwbfile
