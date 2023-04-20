# Dataclass

## BaseData クラス

BaseData クラスは、他のデータクラスの継承元となる基本的なデータクラスです。このクラスには以下の機能があります。

- **`__init__` メソッド**: `file_name` を引数に取り、インスタンス変数に保存します。
- **`save_json` メソッド**: このメソッドは基本的には何も行いませんが、継承したクラスでオーバーライドされることを想定しています。オーバーライドされた場合、指定された `json_dir` に JSON ファイルを保存する機能が実装されることが期待されます。
- **`__del__` メソッド**: インスタンスが削除される際に、ガーベジコレクションを実行してリソースを解放します。

## ImageData クラス

- **ファイル名**: `{file_name}.json` （`file_name` は `__init__` メソッドで指定）
- **内容**: `ImageData` インスタンスの `data` プロパティから生成された画像リスト

## RoiData クラス

- **ファイル名**: `{file_name}.json` （`file_name` は `__init__` メソッドで指定）
- **内容**: `RoiData` インスタンスの `data` プロパティから生成された画像リスト

> **注意**: 両クラスともに `create_images_list` 関数を使用して `data` プロパティから画像リストを生成。実際の画像データは、JSONファイルには直接含まれず、画像リストの情報が格納されている。

## Suite2pData クラス

- **ファイル名**: `{file_name}.npy`（`file_name` は `__init__` メソッドで指定）
- **内容**: `Suite2pData` インスタンスの `data` プロパティの内容を numpy ファイルとして保存

## TimeSeriesData クラス

- **ディレクトリ名**: `{file_name}`（`file_name` は `__init__` メソッドで指定）
- **内容**: `TimeSeriesData` インスタンスの `data` プロパティを元に、セルごとに JSON ファイルを生成
  - **各 JSON ファイル名**: `{cell_number}.json`（`cell_number` はセル番号）
  - **各 JSON ファイル内容**: セルのデータと標準偏差（`std`）が含まれるデータフレーム

## FluoData クラス

- FluoData クラスは TimeSeriesData クラスを継承し、`save_json` メソッドは TimeSeriesData クラスと同様のファイルを出力。

## BehaviorData クラス

- BehaviorData クラスも TimeSeriesData クラスを継承し、`save_json` メソッドは TimeSeriesData クラスと同様のファイルを出力。

## IscellData クラス

IscellData クラスは BaseData クラスを継承していますが、`save_json` メソッドが実装されていないため、このクラスのインスタンスによってファイル出力が行われることはありません。ただし、基本的な構造は以下のようになります。

- **ファイル名**: このクラスでは、`save_json` メソッドが実装されていないため、ファイル名はありません。
- **内容**: `IscellData` インスタンスの `data` プロパティには、`__init__` メソッドで指定されたデータが格納されます。

## HeatMapData クラス

HeatMapData クラスは BaseData クラスを継承しており、`save_json` メソッドを実装しています。このクラスのインスタンスによって生成されるファイル出力は以下のようになります。

- **ファイル名**: `{file_name}.json`（`file_name` は `__init__` メソッドで指定されたもの）
- **内容**: `HeatMapData` インスタンスの `data` プロパティを元に、データフレームが作成され、そのデータフレームが JSON ファイルに書き込まれます。データフレームのカラムは、`columns` プロパティに格納されている値を使用します。

## ScatterData クラス

ScatterData クラスは BaseData クラスを継承しており、`save_json` メソッドを実装しています。このクラスのインスタンスによって生成されるファイル出力は以下のようになります。

- **ファイル名**: `{file_name}.json`（`file_name` は `__init__` メソッドで指定されたもの）
- **内容**: `ScatterData` インスタンスの `data` プロパティを JSON ファイルに書き込みます。このプロパティは、`__init__` メソッドで指定されたデータの転置（`.T`）が格納されています。

## BarData クラス

BarData クラスは BaseData クラスを継承しており、`save_json` メソッドを実装しています。このクラスのインスタンスによって生成されるファイル出力は以下のようになります。

- **ファイル名**: `{file_name}.json`（`file_name` は `__init__` メソッドで指定されたもの）
- **内容**: `BarData` インスタンスの `data` プロパティを JSON ファイルに書き込みます。このプロパティは、`__init__` メソッドで指定されたデータが格納されています。また、`index` プロパティが指定されている場合、データフレームのインデックスとして設定されます。指定されていない場合、デフォルトで連番のインデックスが設定されます。

## HTMLData クラス

HTMLData クラスは BaseData クラスを継承しており、`save_json` メソッドを実装しています。このクラスのインスタンスによって生成されるファイル出力は以下のようになります。

- **ファイル名**: `{file_name}.html`（`file_name` は `__init__` メソッドで指定されたもの）
- **内容**: `HTMLData` インスタンスの `data` プロパティを HTML ファイルに書き込みます。このプロパティは、`__init__` メソッドで指定されたデータが格納されています。

# NWB File I/O

## Initialize

- NWB File Member
  - [Acquisition](https://pynwb.readthedocs.io/en/stable/pynwb.file.html#pynwb.file.NWBFile.add_acquisition)
    - [TwoPhotonSeries](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.TwoPhotonSeries)
  - [Device object](https://pynwb.readthedocs.io/en/stable/pynwb.file.html#pynwb.file.NWBFile.create_device)
  - [Imaging plane](https://pynwb.readthedocs.io/en/stable/pynwb.file.html#pynwb.file.NWBFile.create_imaging_plane)
    - [Optical Channel](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.OpticalChannel)
    - [Device object](https://pynwb.readthedocs.io/en/stable/pynwb.file.html#pynwb.file.NWBFile.create_device)

`data.py`の`file_write.py`ファイルに含まれる`FileWrite`クラスでNWBファイルに記載している。

## suite2p

### suite2p_file_convert

#### NWB output

No access to NWB file

#### Function output

info = {
    'meanImg': ImageData(ops['meanImg'], file_name='meanImg'),
    'ops': Suite2pData(ops, file_name='ops')
}

### suite2p_registration

#### NWB output

No access to NWB file

#### Function output

info = {
    'refImg': ImageData(ops['refImg'], file_name='refImg'),
    'meanImgE': ImageData(ops['meanImgE'], file_name='meanImgE'),
    'ops': Suite2pData(ops, file_name='ops'),
}

### suite2p_roi

#### NWB output

- [ophys](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [ImageSegmentation](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.data_interfaces)
    - [PlaneSegmentation](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation)
      - [pixel_mask](https://suite2p.readthedocs.io/en/latest/api/suite2p.extraction.html#suite2p.extraction.extract.create_masks_and_extract)
      - [columns](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation): 'iscell'
        - [description]: 'two columns - iscell & probcell'
        - [data](https://suite2p.readthedocs.io/en/latest/api/suite2p.classification.html#suite2p.classification.classify.classify)
  - [Fluorescence](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.Fluorescence)
    - [Fluorescence](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.RoiResponseSeries)
      - [data](https://suite2p.readthedocs.io/en/latest/api/suite2p.extraction.html#suite2p.extraction.extract.create_masks_and_extract):F
      - [rois](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation.create_roi_table_region)
  - [Neuropil](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.Fluorescence)
    - [Neuropil](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.RoiResponseSeries)
      - [data](https://suite2p.readthedocs.io/en/latest/api/suite2p.extraction.html#suite2p.extraction.extract.create_masks_and_extract):F_neu
      - [rois](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation.create_roi_table_region)

#### Function output

info = {
    'ops': Suite2pData(ops),
    'max_proj': ImageData(ops['max_proj'], file_name='max_proj'),
    'Vcorr': ImageData(ops['Vcorr'], file_name='Vcorr'),
    'fluorescence': FluoData(F, file_name='fluorescence'),
    'iscell': IscellData(iscell, file_name='iscell'),
    'all_roi': RoiData(np.nanmax(im, axis=0), file_name='all_roi'),
    'non_cell_roi': RoiData(np.nanmax(im[~iscell], axis=0), file_name='noncell_roi'),
    'cell_roi': RoiData(np.nanmax(im[iscell], axis=0), file_name='cell_roi'),
    'nwbfile': nwbfile,
}

### suite2p_spike_deconv

#### NWB output

- [ophys](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [ImageSegmentation](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.data_interfaces)
    - [PlaneSegmentation](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation)
      - [pixel_mask](https://suite2p.readthedocs.io/en/latest/api/suite2p.extraction.h
  - [Deconvolved](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.Fluorescence)
    - [Deconvolved](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.RoiResponseSeries)
      - [data](https://suite2p.readthedocs.io/en/latest/api/suite2p.extraction.html#suite2p.extraction.dcnv.oasis)

#### Function output

info = {
    'ops': Suite2pData(ops),
    'spks': FluoData(spks, file_name='spks'),
    'nwbfile': nwbfile
}

## CaImAn

### caiman_mc

#### NWB output

- [CorrectedImageStack](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.CorrectedImageStack)
  - [correncted](https://pynwb.readthedocs.io/en/stable/pynwb.image.html#pynwb.image.ImageSeries)
    - [external_file](https://pynwb.readthedocs.io/en/stable/pynwb.image.html#pynwb.image.ImageSeries.external_file)
    - [format](https://pynwb.readthedocs.io/en/stable/pynwb.image.html#pynwb.image.ImageSeries.format): 'external'
  - [xy_translation](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.TimeSeries)
    - [data](https://caiman.readthedocs.io/en/latest/core_functions.html#motion-correction)
  - [original] : [Acquisition](https://pynwb.readthedocs.io/en/stable/pynwb.file.html#pynwb.file.NWBFile.add_acquisition)

#### Function output

info = {
    'mc_images': mc_images,
    'meanImg': ImageData(meanImg, file_name='meanImg'),
    'rois': RoiData(rois, file_name='rois'),
    'nwbfile': nwbfile,
}

### caiman_cnmf

#### NWB output

- [ophys](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [ImageSegmentation](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.data_interfaces)
    - [PlaneSegmentation](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation)
      - [image_mask](https://caiman.readthedocs.io/en/latest/core_functions.html#caiman.source_extraction.cnmf.estimates.Estimates)
      - [columns](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation): 'iscell'
        - [description]: 'two columns - iscell & probcell'
        - [data](https://caiman.readthedocs.io/en/latest/core_functions.html#caiman.source_extraction.cnmf.cnmf.CNMF.fit)
  - [RoiResponseSeries](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.Fluorescence)
    - [RoiResponseSeries](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.RoiResponseSeries)
      - [data](https://caiman.readthedocs.io/en/latest/core_functions.html#caiman.source_extraction.cnmf.estimates.Estimates): estimates.C
      - [rois](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation.create_roi_table_region)
  - [Background_Fluorescence_Response](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.Fluorescence)
    - [Background_Fluorescence_Response](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.RoiResponseSeries)
      - [data](https://caiman.readthedocs.io/en/latest/core_functions.html#caiman.source_extraction.cnmf.estimates.Estimates): estimates.f
      - [rois](https://pynwb.readthedocs.io/en/stable/pynwb.ophys.html#pynwb.ophys.PlaneSegmentation.create_roi_table_region)

#### Function output

info = {
    'images': ImageData(np.array(Cn * 255, dtype=np.uint8), file_name='images'),
    'fluorescence': FluoData(fluorescence, file_name='fluorescence'),
    'iscell': IscellData(iscell, file_name='iscell'),
    'all_roi': RoiData(all_roi, file_name='all_roi'),
    'cell_roi': RoiData(cell_roi, file_name='cell_roi'),
    'non_cell_roi': RoiData(non_cell_roi, file_name='non_cell_roi'),
    'nwbfile': nwbfile
}

## Lccd

### lccd_detect

#### NWB output

No access to NWB file

#### Function output

info = {
    'rois': RoiData(np.nanmax(roi_list, axis=0), file_name='cell_roi'),
    'fluorescence': FluoData(timeseries, file_name='fluorescence'),
    'dff': FluoData(timeseries_dff, file_name='dff'),
}

## Edit ROI

#### NWB output

No access to NWB file

#### Function output

save_json_data(ops, im, save_path=node_dirpath,
save_data=['ops', 'fluorescence', 'all_roi', 'non_cell_roi', 'cell_roi']
)

### Add ROI

#### NWB output

No access to NWB file

#### Function output

save_json_data(ops, im, save_path=node_dirpath,
save_data=['ops', 'fluorescence', 'all_roi', 'non_cell_roi', 'cell_roi'])

### Delete ROI

#### NWB output

No access to NWB file

#### Function output

save_json_data(ops, im, save_path=node_dirpath,
save_data=['ops', 'fluorescence', 'all_roi', 'non_cell_roi', 'cell_roi'])

### Merge ROI

#### NWB output

No access to NWB file

#### Function output

save_json_data(ops, im, save_path=node_dirpath,
save_data=['ops', 'fluorescence', 'all_roi', 'non_cell_roi', 'cell_roi'],
)

## Optinistモジュール

### Basic neural analysis

#### eta

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [mean](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [sem](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [num_sample](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info['mean'] = TimeSeriesData(
    mean,
    std=sem,
    index=list(np.arange(params['start_time'], params['end_time'])),
    cell_numbers=cell_numbers if iscell is not None else None,
    file_name='mean'
)
info['mean_heatmap'] = HeatMapData(
    norm_mean,
    columns=list(np.arange(params['start_time'], params['end_time'])),
    file_name='mean_heatmap'
)
info['nwbfile'] = nwbfile

### Dimension reduction

#### cca

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [projectedNd](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [x_weights](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [y_weights](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [x_loadings_](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [y_loadings_](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [coef](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [n_iter_](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'projectedNd': ScatterData(proj, file_name='projectedNd'),
    'coef': BarData(cca.coef_.flatten(), file_name='coef'),
    'nwbfile': nwbfile,
}

#### pca

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [pca_projectedNd](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [components](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [explained_variance](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [explained_variance_ratio](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [singular_values](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [mean](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [n_components](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [n_samples](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [noise_variance](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [n_features_in](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'explained_variance': BarData(pca.explained_variance_ratio_, file_name='evr'),
    'projectedNd': ScatterData(proj_X, file_name='projectedNd'),
    'contribution': BarData(
        pca.components_,
        index=[f'pca: {i}' for i in range(len(pca.components_))],
        file_name='contribution'
    ),
    'cumsum_contribution': BarData(
        np.cumsum(pca.components_, axis=0),
        index=[f'pca: {i}' for i in range(len(pca.components_))],
        file_name='cumsum_contribution'
    ),
    'nwbfile': nwbfile,
}

#### tsne

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [projectedNd](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'projectedNd': ScatterData(proj_X, file_name='projectedNd'),
    'nwbfile': nwbfile,
}

### Neural decoding

#### glm

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [actual_predicted](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [params](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [pvalues](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [tvalues](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [aic](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [bic_llf](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [llf](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [pearson_chi2](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [df_model](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [df_resid](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [scale](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [mu](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [endog](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
        'actual_predicted': ScatterData(
            np.array([Res._endog, Res.mu]).transpose(),
            file_name='actual_predicted'
        ),
        'params': BarData(Res.params.values, file_name='params'),
        'textout': HTMLData(Res.summary().as_html(), file_name='textout'),
        'nwbfile': nwbfile,
}

#### lda

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [score](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'score': BarData(score, file_name='score'),
    'nwbfile': nwbfile
}

#### svm

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [score](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'score': BarData(score, file_name='score'),
    'nwbfile': nwbfile
}

### Neural population analysis

#### correlation

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [corr](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'corr': HeatMapData(corr, file_name='corr'),
    'nwbfile': nwbfile,
}

#### cross_correlation

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [mat](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [baseline](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [base_confint](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)

##### Function output

info = {
    'nwbfile': nwbfile
}
info[name] = TimeSeriesData(arr1.T, file_name=name)
info[name] = TimeSeriesData(arr2.T, file_name=name)


#### granger

##### NWB output

- [optinist](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule)
  - [Granger_fval_mat](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [gc_combinations](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [gc_ssr_ftest](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [gc_ssr_chi2test](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [gc_lrtest](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [gc_params_ftest](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [cit_pvalue](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)
  - [adf_pvalue](https://pynwb.readthedocs.io/en/stable/pynwb.base.html#pynwb.base.ProcessingModule.add_container)


##### Function output

info['Granger_fval_mat_heatmap'] = HeatMapData(
    GC['Granger_fval_mat'][0],
    file_name='gfm_heatmap'
)
info['Granger_fval_mat_scatter'] = ScatterData(
    GC['Granger_fval_mat'][0],
    file_name='gfm'
)
info['nwbfile'] = nwbfile
