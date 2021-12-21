

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
    nwbfile.create_imaging_plane(
        name=nwb_dict['imaging_plane']['name'],
        description=nwb_dict['imaging_plane']['description'],
        optical_channel=optical_channel,   # 光チャネル
        device=device,   # 電極デバイス
        imaging_rate=nwb_dict['imaging_plane']['imaging_rate'],   # 画像の比率Hz
        excitation_lambda=nwb_dict['imaging_plane']['excitation_lambda'], # 励起（れいき）波長
        indicator=nwb_dict['imaging_plane']['indicator'],   # カルシウムインディケーター
        location=nwb_dict['imaging_plane']['location'],
    )

    image_series = ImageSeries(
        name=nwb_dict['image_series']['name'],
        external_file=nwb_dict['image_series']['external_file'],
        format='external',
        rate=nwb_dict['imaging_plane']['imaging_rate'],
        starting_frame=[nwb_dict['image_series']['starting_frame']]
    )

    nwbfile.add_acquisition(image_series)

    return nwbfile


def nwb_motion_correction(nwb, mc_data, xy_trans_data):
	# motion correction
	original = ImageSeries(
		name='original',  # this must be named "corrected"
		data=list(nwb['nwbfile'].acquisition.values())[0],
		unit='na',
		format='raw',
		starting_time=0.0,
		rate=1.0
	)

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
		original=original,
		xy_translation=xy_translation,
	)

	motion_correction = MotionCorrection(
		corrected_image_stacks=corrected_image_stack
	)
	nwb['ophys'].add(motion_correction)

    return nwb

