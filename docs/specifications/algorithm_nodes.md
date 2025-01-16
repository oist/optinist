(algorithm_nodes)=
Algorithm Nodes
=================

OptiNiSt includes a variety of third-party calcium (Ca<sup>2+</sup>) imaging software packages for processing video data and extracting cell ROI information. OptiNiSt also includes a selection of data analysis tools. Each step of data processing and analysis are contained within a Node, which can be mixed-and-matched. Here we summarise the usage and requirements of each node.


#### Data Processing
##### [CaImAn](https://caiman.readthedocs.io/en/latest/)
###### caiman_mc
  - **Description:** Motion Correction ([NoRMCorre](https://caiman.readthedocs.io/en/latest/CaImAn_features_and_references.html#normcorre)). The function will perform either rigid or piecewise rigid motion correction depending on the attribute self.pw_rigid and will perform high pass spatial filtering for determining the motion (used in 1p data) if the attribute self.gSig_filt is not None.
  - **Input:** ImageData (Tiff), HDF5
  - **Output:** ImageData
  - **Parameters:**
    - **border_nan** [true, false, copy, default: copy]: flag for allowing NaN in the boundaries. True allows NaN, whereas ‘copy’ copies the value of the nearest data point.
    - **g_Sig_filt**　[int or None, default: None]:　Size of kernel for high-pass spatial filtering in 1P data. If None no spatial filtering is performed.
    - **is3D** [bool, default: False]: Flag for 3D recordings for motion correction
    - **max_deviation_rigid** [int, default: 3]: Maximum deviation in pixels between rigid shifts and shifts of individual patches
    - **max_shifts** [(int, int), default: (6,6)]: Maximum shifts per dimension in pixels.
    - **min_mov** [float or None, default: None]: Minimum value of movie, used to adjust the brightness or offset. If None it gets computed.
    - **niter_rig** [int, default: 1]: Number of iterations rigid motion correction.
    - **nonneg_movie** [bool, default: True]: Ensures that the minimum value in the motion-corrected movie is non-negative by subtracting the minimum value if it's negative.
    - **num_frames_split** [int, default: 80]: Splits movie every x frames for parallel processing
    - **overlaps** [(int, int), default: (32, 32)]: Overlap between patches in pixels in piecewise-rigid motion correction.
    - **pw_rigid** [ bool, default: False]: flag for performing piecewise-rigid motion correction.
    - **splits_els** [int, default: 14]: Number of splits across time for pw-rigid registration.
    - **splits_rig** [int, default: 14]: Number of splits across time for rigid registration
    - **shifts_opencv** [bool, default: True]: Flag for applying shifts using cubic interpolation (otherwise FFT).
    - **strides** [(int, int), default: (96, 96)]: How often to start a new patch in pw-rigid registration. Size of each patch will be strides + overlaps
    - **upsample_factor_grid** [int, default: 4]:, Motion field upsampling factor during FFT shifts.
    - **use_cuda** [bool, default: False]: Flag for using a GPU.
    - **advanced** (See [CaImAn documentation](https://caiman.readthedocs.io/en/latest/core_functions.html#))
###### caiman_cnmf
  - **Description:** Constrained nonnegative matrix factorization, find ROIs 2P data
  - **Input:**　Image (Tiff), HDF5
  - **Output:** FluoData, IsCellData
  - **Parameters:**
    - **decay_time** [float, default: 0.4] Length of typical transient in seconds. Dependent on fluorescent indicator dynamics.
    - **use_online** Provides faster processing and reduced memory requirements, slight reduction in accuracy. Good for large data.
    - **do_refit** Refit reprocesses data after an initial fit, refining results by using the centres of components from the initial fit as seeds for another iteration.
    - **K [int, default: 10]: Number of components to be found. This determines how many cells can be found** (per patch or whole FOV depending on whether rf=None)
    - **gSig** [default: (4,4)]: Standard deviation of neuron size along x and y.
    - **ssub** [float, default: 2]: Spatial downsampling factor.
    - **tsub** [float, default 1, no downsampling]: temporal downsampling factor recommended for long datasets.
    - **nb** [int, default: 2]: Number of global background components (noise complexity).
    - **method_init** [‘greedy_roi’,’corr_pnr’, default: ‘greedy_roi’]: Use ‘greedy_roi’ for standard 2P, ‘corr_pnr’ for 1p (CNMF-E).
    - **roi_thr** [scalar between 0 and 1, default: 0.9]: Energy threshold for computing contours
    - **p** [int, default: 2]: Order of the autoregression model.
    - **rf** [int or None, default: None]: Half-size of patch in pixels. If None, no patches are constructed and the whole FOV is processed jointly.
    - **stride** [int or 0, default: 6]: Overlap between neighboring patches in pixels. Used to optimize memory consumption and parallelizing computations, as it allows for the processing of data in smaller segments while maintaining some continuity between them.
    - **merge_thr** [float, default: 0.8]: Trace correlation threshold for merging two components. A higher threshold value means that only components with very similar activity traces will be merged.
    - **advanced** (See [CaImAn documentation](https://caiman.readthedocs.io/en/latest/core_functions.html#))
###### caiman_cnmfe
  - **Description:** Constrained nonnegative matrix factorization-endoscope, find ROIs for 1P data.
  - **Input:**　Image (Tiff), HDF5
  - **Output:** FluoData, IsCellData
  - **Parameters:**
    - Additional information on CNMF-E can be found on [Inscopix github](https://github.com/inscopix/inscopix-cnmfe?tab=readme-ov-file#algorithm-parameters)
    - **decay_time** [float, default: 0.4] Length of typical transient in seconds. Dependent on fluorescent indicator dynamics
    - **do_refit** Refit reprocesses data after an initial fit, refining results by using the centres of components from the initial fit as seeds for another iteration.
    - **K [int, default: 10]: Number of components to be found. This determines how many cells can be found** (per patch or whole FOV depending on whether rf=None)
    - **gSig** [default: (4,4)]: Standard deviation of neuron size along x and y.
    - **ssub** [float, default: 2]: Spatial downsampling factor.
    - **tsub** [float, default 1, no downsampling]: temporal downsampling factor recommended for long datasets.
    - **nb** [int, default: 2]: Number of global background components (noise complexity).
    - **roi_thr** [scalar between 0 and 1, default: 0.9]: Energy threshold for computing contours
    - **p** [int, default: 2]: Order of the autoregression model.
    - **rf** [int or None, default: None]: Half-size of patch in pixels. If None, no patches are constructed and the whole FOV is processed jointly.
    - **stride** [int or 0, default: 6]: Overlap between neighboring patches in pixels. Used to optimize memory consumption and parallelizing computations, as it allows for the processing of data in smaller segments while maintaining some continuity between them.
    - **merge_thr** [float, default: 0.8]: Trace correlation threshold for merging two components. A higher threshold value means that only components with very similar activity traces will be merged.

###### caiman_cnmf_multisession
  - **Description:** Keep union of components found in multiple sessions
  - **Input:**　Image (Tiff), HDF5
  - **Output:** FluoData, IsCellData
  - **Parameters:** Same as caiman_cnmf with the addition of:
    - n_reg_files
    - reg_file_rate
    - **align_flag** [bool, default: true]: Align the templates before matching
    - **max_thr** [scalar, default 0]: Max threshold parameter before binarization.
    - **use_opt_flow** [bool, default: true]: Use dense optical flow to align templates
    - **thresh_cost** [scalar, default: 0.7]: Maximum distance considered.
    - **max_dist** [scalar, default: 10]: Maximum distance between centroids.
    - **enclosed_thr** [float or none, default: None]: If not None, set distance to at most the specified value when ground truth is a subset of inferred.
##### [Suite2p](https://suite2p.readthedocs.io/en/latest/)

###### suite2p_file_convert
  - **Description:** Convert to int16
  - **Input:** Image (Tiff), HDF5
  - **Output:** Suite2pData
  - **Parameters:**
    - force sktiff [boolean, default: False] Specifies whether or not to use scikit-image for reading in tiffs
    - batch_size [int, default: 500]: How many frames to register simultaneously in each batch. This depends on memory constraints it will be faster to run if the batch is larger, but it will require more RAM.

###### suite2p_registration
  - **Description:** Motion correction
  - **Input:** Suite2pData
  - **Output:** Suite2pData, ImageData
  - **Parameters:**
    - **I/O:**
        - **nplanes** [int, default: 1]: each tiff has this many planes in sequence
        - **nchannels** [int, default: 1]: each tiff has this many channels per plane
        - **functional_chan** [int, default: 1]: this channel is used to extract functional ROIs (e.g. 1 = first channel, and 2 = second channel)
        **frames_include** [int, default: -1]: If greater than zero, only frames_include frames are processed. useful for testing parameters on a subset of data.

    - **Registration**
        - **smooth_sigma** [float, default: 1.15]: Standard deviation in pixels of the gaussian used to smooth the phase correlation between the reference image and the frame which is being registered. A value of >4 is recommended for one-photon recordings (with a 512x512 pixel FOV).
        - **smooth_sigma_time**: [float, default: 0]: Standard deviation in time frames of the gaussian used to smooth the data before phase correlation is computed. Might need this to be set to 1 or 2 for low SNR data.
        - **maxregshift** [float, default: 0.1]: Maximum shift as a fraction of the frame size. If the frame is Ly pixels x Lx pixels, then the maximum pixel shift in pixels will be max(Ly,Lx) * ops[‘maxregshift’].
        - **align_by_chan** [int, default: 1]: Which channel to use for alignment (1-based, so 1 means 1st channel and 2 means 2nd channel). If you have a non-functional channel with something like td-Tomato expression, you may want to use this channel for alignment rather than the functional channel.
        - **reg_tif** [bool, default: False]: Whether or not to write the registered binary to tiff files
        - **th_badframes** [float, default: 1.0]: Involved with setting threshold for excluding frames for cropping. Set this smaller to exclude more frames.

    - **1P:**
        - **1Preg** [bool, default: False]: Whether to perform high-pass spatial filtering and tapering (parameters set below), which help with 1P registration.
        - **spatial_hp_reg** [int, default: 42]: Window in pixels for spatial high-pass filtering before registration.
        - **pre_smooth** [float, default: 0]: If > 0, defines stddev of Gaussian smoothing, which is applied before spatial high-pass filtering.
        - **spatial_taper** [float, default: 40]: How many pixels to ignore on edges - they are set to zero (important for vignetted windows, for FFT padding do not set BELOW 3*ops[‘smooth_sigma’]).

    - **Non-rigid:**
        - **nonrigid** [bool, default: True]: Whether to use nonrigid registration, which splits the field of view into blocks and computes registration offsets in each block separately.
        - **block_size** [list of int, default: [128, 128]]: Size of blocks for non-rigid registration, in pixels. HIGHLY recommend keeping this a power of 2 and/or 3 (e.g. 128, 256, 384, etc) for efficient fft.
        - **snr_thresh** [float, default: 1.2]: How big the phase correlation peak has to be relative to the noise in the phase correlation map for the block shift to be accepted. In low SNR recordings like one-photon, I’d recommend a larger value like 1.5, so that block shifts are only accepted if there is significant SNR in the phase correlation.
        - **maxregshiftNR** [int, default: 5]: Maximum shift in pixels of a block relative to the rigid shift.

###### suite2p_roi
  - **Description:** ROI detection 2P or 1P
  - **Input:** Suite2pData
  - **Output:** Suite2pData, FluoData, IsCellData
  - **Parameters:**
    - **tau** [float, default 1.0]: The timescale of the sensor (in seconds), used for deconvolution kernel. The kernel is fixed to have this decay and is not fit to the data.
      - 0.7 for GCaMP6f
      - 1.0 for GCaMP6m
      - 1.25-1.5 for GCaMP6s

    - **Classification:**
        - **soma_crop** [bool, default,: True] Crop dendrites for cell classification stats like compactness

    - **Cell_detection:**
      - **high_pass** [int, default: 100]: Running mean subtraction across time with window of size ‘high_pass’. Values of less than 10 are recommended for 1P data where there are often large full-field changes in brightness.
      - **sparse_mode** [bool, default: True]: Whether or not to use sparse_mode cell detection.
      - **max_overlap** [float, default: 0.75] Overlapping ROIs are allowed during cell detection. After detection, ROIs with more than ops[‘max_overlap’] fraction of their pixels overlapping with other ROIs will be discarded. Therefore, to throw out NO ROIs, set this to 1.0.
      - **nbinned** [int, default: 5000]: Maximum number of binned frames to use for ROI detection.
      - **spatial_scale** [int, default: 0]: If set to 0, then the algorithm determines the optimal scale of the recording is in pixels automatically. If it seems off, set it yourself to the following values: 1 (=6 pixels), 2 (=12 pixels), 3 (=24 pixels), or 4 (=48 pixels).
      - **threshold_scaling** [float, default: 1.0]: Controls the threshold at which to detect ROIs (how much the ROIs have to stand out from the noise to be detected). If higher, then fewer ROIs will be detected, and if lower, more ROIs will be detected.
      - **max_iterations** [int, default: 20] How many iterations over which to extract cells - at most ops[‘max_iterations’], but usually stops before due to ops[‘threshold_scaling’] criterion.
      - **spatial_hp_detect** [int, default: 25]: Window for spatial high-pass filtering for neuropil subtracation before ROI detection takes place.

    - **Output:**
      - **preclassify** [float, default: 0.0] Apply classifier before signal extraction with probability threshold of “preclassify”. If this is set to 0.0, then all detected ROIs are kept and signals are computed.

    - **ROI_extraction:**
      - **allow_overlap** [bool, default: False]: Whether or not to extract signals from pixels which belong to two ROIs. By default, any pixels which belong to two ROIs (overlapping pixels) are excluded from the computation of the ROI trace.
      - **min_neuropil_pixels** [int, default: 350]: Minimum number of pixels used to compute neuropil for each cell.

###### suite2p_spike_deconv
  - **Description:**  Estimates voltage spikes from calcium imaging data
  - **Input:** FluoData
  - **Output:** FluoData
  - **Parameters:**
    - **neucoeff** [float, default: 0.7]: Neuropil coefficient for all ROIs. Used to subtract neuropil contamination from cell fluorescence.
    - **baseline** [string,'maximin',‘constant’,'prctile', default: maximin]: Baseline of each trace to be subtracted from each cell.
      - maximin: computes a moving baseline by filtering the data with a Gaussian of width ops['sig_baseline'] * ops['fs'], and then minimum filtering with a window of ops['win_baseline'] * ops['fs'], and then maximum filtering with the same window.
      - constant: computes a constant baseline by filtering with a Gaussian of width ops['sig_baseline'] * ops['fs'] and then taking the minimum value of this filtered trace.
      - prctile: computes a constant baseline by taking the ops['prctile_baseline'] percentile of the trace.
        <!--Note to check these are actually possible-->
    - **win_baseline** [float, default: 60.0]: Window for maximin filter in seconds
    - **sig_baseline** [float, default: 10.0] Gaussian filter width in seconds, used before maximin filtering or taking the minimum value of the trace, ops['baseline'] = 'maximin' or 'constant'.
    - **prctile_baseline** [float, optional, default: 8]: percentile of trace to use as baseline if ops['baseline'] = 'constant_percentile'.

##### [LCCD](https://doi.org/10.1016/j.neures.2022.02.008)
###### lccd_detect
  - **Description:** Low Computational-cost Cell Detection
  - **Input:** ImageData
  - **Output:** FluoData, IsCellData
  - **Parameters:**
    - **blob_detector:** Initial step that identifies potential cells or "blobs" in each frame of the video.
      - **filtersize1** [int (filtersize1 <= Time / 2), default: 100]: Size of the larger temporal filter for smoothing.
      - **filtersize2** [int (filtersize2 < filtersize1), default: 4]: Size of the smaller temporal filter for edge detection.
      - **sigma** [float, smaller than T & filtersize2 < filtersize1, default: 1.25]: Standard deviation for the Gaussian filter.
      - **fsize** [int (odd number recommended), default: 9]: Size of the 2D Gaussian filter for spatial filtering.
      - **min_area** [int (min_area < max_area < X * Y), default: 16]: Define the acceptable size range for detected blobs.
      - **max_area** [int, (max_area < X * Y), default: 100]: Define the acceptable size range for detected blobs.
      - **sparse** [bool, default: False]: Determines whether to use sparse matrix operations.

    - **roi_integration:** Takes blobs detected across all frames and integrates them into coherent ROI that represent consistent cell locations over time.
      - **overlap_threshold** [float (0 < overlap_threshold < 1), default: 0.25]: Threshold for determining significant overlap between ROIs
      - **min_area** [int (min_area < max_area < X * Y), default: 16]: Define the acceptable size range for integrated ROIs.
      - **max_area** [int (max_area < X * Y), default: 100]: Maximum acceptable size for integrated ROIs.
      - **sparse** [bool, default: False]: Determines whether to use sparse matrix operations.

    - **lccd:**
      - **frame_divider** [int (Time / frame_divider >= 2), default: 100]: Divides the time dimension into segments for faster processing.

    - **dff:**
      - **f0_frames** [int (0, f0_frames, Time), default: 100]: Number of frames used to calculate the baseline fluorescence.
      - **f0_percentile** [float (0 <= f0_percentile <= 100), default: 8]: Percentile used for calculating the baseline fluorescence.

#### OptiNiSt (Analysis)
##### Basic neural analysis
###### ETA (Event-Triggered Average)
  - **Description:** Calculates the average neural response around specific events in behavioral data.
  - **Input:** FluoData, BehaviorData, IsCellData (optional)
    - **Neural data (X) and behavior data (Y) must have the same number of time points: X.shape[0] == Y.shape[0].**
  - **Output:** mean (TimeSeriesData), mean_heatmap (HeatMapData), nwbfile
  - **Parameters:**
    - **transpose_x** [bool, default: true]: Whether to transpose the neural data.
    - **transpose_y** [bool, default: false]: Whether to transpose the behaviour data.
    - **event_col_index** [int, default 1]: Index of column in behavioral data to use for event detection.
    - **trigger_type** ['up', 'down', 'cross', default: 'up']:
      - 'up' detects transitions 0 to trigger_threshold
      - 'down' detects transitions trigger_threshold to 0
      - 'cross' detects either up or down transitions.
    - **trigger_threshold** [float, default 0.5]: Threshold value for trigger detection
    - **pre_event** [int (0 < pre_event < T/2), default: -10]: Number of time points before the trigger to include.
    - **post_event** [int (T/2 > post_event < 0), default: 10]: Number of time points after the trigger to include.

##### Dimensionality Reduction
###### [CCA](https://scikit-learn.org/stable/modules/generated/sklearn.cross_decomposition.CCA.html) (Canonical Correlation Analysis)
* **Description:** Finds linear combinations of variables in two datasets that have maximum correlation with each other
* **Input:** FluoData (neural activity), BehaviorData (behavioral variables), IsCellData (optional)
    - **Neural data (X) and behavior data (Y) must have the same number of time points: X.shape[0] == Y.shape[0].**
    - **Y must be a single column after selecting the target_index.**
* **Output:** ScatterData (projected data), BarData (correlation coefficients).
* **Parameters:**
    - **I/O:**
      - **target_index** [int (columns start from 0), default 1]: Index of the target behavioral variable.
      - **transpose_x** [bool, default: False]: Whether to transpose the neural data matrix.
      - **transpose_y** [bool, default: False]: Whether to transpose the behavior data matrix.
      - **standard_x_mean** [bool, default: True]: Whether to standardize X by subtracting mean.
      - **standard_x_std** [bool, default: True]: Whether to standardize X by dividing by std.
      - **standard_y_mea** [bool, default: True]: Whether to standardize Y by subtracting mean.
      - **standard_y_std** [bool, default: True]: Whether to standardize Y by dividing by std.

    - **CCA:**
      - **n_components** [int, 1 <= n_components <= min(n_features, n_samples), default: 2]: Number of components to keep.
      - **scale** [bool, default: True]: Whether to scale the data to unit variance before performing CCA. Ensures that all variables contribute equally to the analysis.
      - **max_iter** [int, default: 500]: Maximum number of iterations of the power method.
      - **tol** [float, default: 0.000001]: Convergence tolerance for the power method.  Will stop iterating when the change between iterations is less than this value.
      - **copy** [bool, default: True]: Whether to copy X and Y before fitting and transforming. Ensures that the original data is not modified during the CCA process.


###### [PCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html)
  - **Description:** Dimensionality-reduction, linear reduction in components while retaining max variance.
  - **Input:** FluoData (neural activity) , IsCellData (optional)
    - **Neural data shape should be (timepoints, ROI)**
  - **Output:** ScatterData (projected data), BarData (explained variance ratios)
  - **Parameters:**
    - **I/O:**
      - **transpose** [bool, default: False]: Whether to transpose the neural data matrix.
      - **standard_mean** [bool, default: True]: Whether to standardize X by subtracting mean.
      - **standard_std** [bool, default: True]: Whether to standardize X by dividing by std.

    - **PCA:**
      - **n_components** [int, float, default: None]: Number of components to keep. If int, it specifies the exact number of components. If float between 0 and 1, it specifies the proportion of variance to be retained.

    - **advanced:**
      - See [SciKit link](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html) for more information.
      - copy=True, whiten=False, svd_solver='auto', tol=0.0, iterated_power='auto'
      <!-- check which other scikit can be made simpler by splitting to advanced-->

###### [dPCA](https://github.com/machenslab/dPCA) (demixed Principal Component Analysis)
  - **Description:** Decompose neural activity into components related to task-related components.
   - **Input:** FluoData (neural activity), BehaviorData (behavioral variables) , IsCellData (optional)
     - **Neural data should have shape (n_neurons, n_timepoints).**
     - **Behavior data should have shape (n_timepoints, n_behavioral_variables).**
     - **Behavior data should include columns for trigger_column and feature_columns**
     - **Set trigger_column and feature_column indices within total number of behavioural data columns**
  - **Output:** HeatMapData (component activations for different task conditions).
  - **Parameters:**
    - **I/O:**
      - **transpose** [bool, default: False]: Whether to transpose the neural data matrix.
      - **trigger_column** [int, default: 1]: Column index in behavior data for event triggers.
      - **trigger_type** ['up', 'down', 'cross', default: 'up']:
        - 'up' detects transitions 0 to trigger_threshold
        - 'down' detects transitions trigger_threshold to 0
        - 'cross' detects either up or down transitions.
      - **trigger_threshold** [float, default 0.5]: Threshold value for trigger detection
      - **trigger_duration** [list of 2 ints, default: [10, 10]]: Frames before and after trigger to include.
      - **feature_columns** [list of ints, default: [3, 4]]: Columns in behaviors_data to use as features/conditions for the dPCA analysis.

      - **standard_mean** [bool, default: True]: Whether to standardize by subtracting mean.
      - **standard_std** [bool, default: True]: Whether to standardize by dividing by std.

    - **dPCA:**
    <!-- check 't', 'b', 'c', 'tbc', 's', 'ts'. Are these all used correctly?-->
      - **labels** [str, default: 'tbc']: list of characters with which to describe the parameter axes, e.g. 'tsd' to denote time, stimulus and decision axis. All marginalizations (e.g. time-stimulus) are referred to by subsets of those characters (e.g. 'ts'). Must match the number of dimensions in the data.
      - **regularizer** [float, default: 0.001]: If > 0, the regularization weight is regularizer*var(data). Helps prevent overfitting by adding a penalty term to the optimization objective.
      <!-- check if can use 'auto' regulariser-->
      - **n_components** [int or dict, default: 10]: Number of components to keep. If int, same number of components are kept in each marginalization (e.g. {'t' : 10, 'ts' : 5}).
      - **join** [dict or None, default: None]: How to join task parameters. If a data set has parametrized by time t and stimulus s, then dPCA will split the data into marginalizations corresponding to 't', 's' and 'ts'. At times, we want to join different marginalizations (like 's' and 'ts'), e.g. if we are only interested in the time-modulated stimulus components. In this case,we would pass {'ts' : ['s','ts']}.
      - **n_iter:** [int, default: 0]: Number of iterations for randomized SVD solver (sklearn).
      - **copy** [bool, default: True]: Whether to copy X and Y before fitting and transforming. Ensures that the original data is not modified.

    - **Plot:**
      - **figure_components** [list of ints, default: [0, 1]]: Specifies which dPCA components to plot. Default will create figures for the first two components (0 and 1).
      - **figure_features**: ['t', 'b', 'c', 'tbc']: Defines which marginalizations or feature combinations to plot. 't' might represent time, 'b' and 'c' the two behavioral features, and 'tbc' their interaction.

###### [TSNE](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html) (t-distributed stochastic neighbor embedding)
  - **Description:** Reduces dimensionality while preserving local structure of high-dimensional data.
  - **Input:** FluoData (neural activity), IsCellData (optional)
    - **Neural data should have shape (n_neurons, n_timepoints).**
  - **Output:** ScatterData (low-dimensional embedding)
  - **Parameters:**
    - **I/O:**
      - **transpose** [bool, default: False]: Whether to transpose the neural data matrix.
      - **standard_mean** [bool, default: True]: Whether to standardize X by subtracting mean.
      - **standard_std** [bool, default: True]: Whether to standardize X by dividing by std.

    - **TSNE**
      - **n_components** [int, default: 2]: Dimension of the embedded space, i.e., the number of dimensions in the output. Most commonly set to 2 or 3 for visualization purposes.
      - **perplexity** [float, 5 <= perplexity <= n_samples - 1, default: 30]: Related to the number of nearest neighbors used. It can be thought of as a guess about the number of close neighbors each point has. Perplexity balances attention between local and global aspects of the data. Typical values range between 5 and 50
      - **early_exaggeration** [float, default: 12]: Controls how tight natural clusters are in the embedded space
      - **learning_rate** [float or 'auto', default: 'auto']: The learning rate for t-SNE

    - **advanced:**
      - See [SciKit link](https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html) for more information.
      -   n_iter: 1000, n_iter_without_progress: 300, min_grad_norm: 0.0000001, metric: 'euclidean', init: 'pca', random_state: 0, method: 'barnes_hut', angle: 0.5



##### Neural population analysis
###### [GLM](https://www.statsmodels.org/stable/glm.html) (Generalized Linear Model)
  - **Description:** Generalized Linear Model. Extends linear regression to allow for non-normal distributions of the response variable. Can model relationships between neural activity and behavioral variables.
  - **Input:** FluoData (neural activity), BehaviorData (behavioral variables), IsCellData (optional)
  - **Output:** BarData (model parameters), ScatterData (actual vs predicted values), HTMLData (model summary statistics)
  - **Parameters:**
    - **I/O:**
      - **target_index** [int (columns start from 0), default 1]: Column index of behavioral variable to model.
      - **transpose_x** [bool, default: True]: Whether to transpose the neural data matrix.
      - **transpose_y** [bool, default: True]: Whether to transpose the behavior data matrix.
      - **standard_x_mean** [bool, default: True]: Standardize X by subtracting mean.
      - **standard_x_std*** [bool, default: True]: Standardize X by dividing by std.
      - **standard_y_mean** [bool, default: True]: Standardize Y by subtracting mean.
      - **standard_y_std** [bool, default: False]: Standardize Y by dividing by std.

    - **GLM:**
      - **family** [str ('Binomial', 'Gamma', 'Gaussian', 'InverseGaussian', 'NegativeBinomial', 'Poisson', 'Tweedie'), default: 'Gaussian']: Distribution family
      - **link** [str ('CDFLink', 'CLogLog', 'LogLog', 'LogC', 'Log', 'Logit', 'NegativeBinomial', 'Power', 'Cauchy', 'Identity', 'InversePower', 'InverseSquared', 'Probit'), default: 'log']: Link function
      - **add_constant** [bool, default: False]: Add a constant term to the model.
      - **offset** [Array or None]: An offset to be included in the model. If provided, must be an array whose length is the number of rows in exog.
      - **exposure** [Array or None]: Log(exposure) will be added to the linear prediction in the model. Exposure is only valid if the log link is used. If provided, it must be an array with the same length as endog.
      - **missing** [str (‘none', ‘drop’, ‘raise’), default: 'raise']: If ‘none’, no nan checking is done. If ‘drop’, any observations with nans are dropped. If ‘raise’, an error is raised. Default is ‘none’.

###### [LDA](https://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html) (Linear Discriminant Analysis)
  - **Description:** Finds linear combinations of features that best separate classes. Can be used for dimensionality reduction and classification of neural data,  by projecting input to the most discriminative directions. The model fits a Gaussian density to each class, assuming that all classes share the same covariance matrix.
  - **Input:** FluoData (neural activity), BehaviorData (class labels), IsCellData (optional)
    - **Neural data (X) and behavior data (Y) must have the same number of time points: X.shape[0] == Y.shape[0].**
    - **To predict 2 classes, use binary, for 3+ classes use int, e.g. 0,1,2,3**
  - **Output:** BarData (cross-validation scores)
  - **Parameters:**
    - **I/O:**
      - **target_index** [int (columns start from 0), default 1]: Index of target behavioral variable to use as class labels.
      - **transpose_x** [bool, default: True]: Whether to transpose the neural data matrix.
      - **transpose_y** [bool, default: True]: Whether to transpose the behavior data matrix.
      - **standard_x_mean** [bool, default: True]: Standardize X by subtracting mean.
      - **standard_x_std*** [bool, default: True]: Standardize X by dividing by std.

    - **cross_validation:**
      - **n_splits** [int, default: 5]: Number of folds for cross-validation.
      - **shuffle** [bool, default: False]: Whether to shuffle the data before splitting into batches.

    - **LDA:**
      - **solver** ['svd', 'lsqr', 'eigen', default: 'svd']: Solver to use.
        - 'svd': Singular value decomposition, doesn't compute the covariance matrix, therefore this solver is recommended for data with a large number of features.
        - 'lsqr': Least squares solution. Can be combined with shrinkage.
        - 'eigen': Eigenvalue decomposition. Can be combined with shrinkage.
      - **shrinkage** shrinkage [float or 'auto', default: None]: Shrinkage parameter, between 0 and 1. Only available if solver is 'lsqr' or 'eigen'.
      - **priors** [array-like of shape (n_classes,), default: None]: Class priors. By default, the class proportions are inferred from the training data.
      - **n_components** :
      - **store_covariance**:
      - **tol**  : 0.0001*
      - **covariance_estimator**:

###### [SVM](https://scikit-learn.org/stable/modules/svm.html#svm) (Support Vector Machine)
  - **Description:**  Finds a hyperplane that best separates classes in high-dimensional space. Can be used for classification of neural states.  It can handle non-linear relationships through kernel functions and is effective for both linear and non-linear classification tasks.
  - **Input:** FluoData (neural activity), BehaviorData (class labels), IsCellData (optional)
  - **Output:** BarData (cross-validation scores)
  - **Parameters:**
      - **I/O:**
      - **target_index** [int (columns start from 0), default 1]: Index of target behavioral variable to use as class labels.
      - **transpose_x** [bool, default: True]: Whether to transpose the neural data matrix.
      - **transpose_y** [bool, default: True]: Whether to transpose the behavior data matrix.
      - **standard_x_mean** [bool, default: True]: Standardize X by subtracting mean.
      - **standard_x_std*** [bool, default: True]: Standardize X by dividing by std.

    - **Grid Search**: Find optimal parameters by comparing models within specified range
        (Otherwise specify yourself in SVM parameters)
      - **use_grid_search** [bool, default: True]: Whether to perform grid search for hyperparameter optimization.
      - **params_to_grid_search**: Parameters and their ranges to search during grid search:
        - **C** [list of float, default: [0.001, 0.01, 0.1]]: Range of regularization parameter (inverse of regularization strength).
        - **kernel** [list of str (‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’), default: ['linear', 'rbf']]: Kernel types to test; (See SVM below for more info).
        - **degree** [list of int, default: [2,3]]: Polynomial degree for 'poly' kernel. Degree of the polynomial kernel function. Must be non-negative. Ignored by all other kernels.
        - **gamma** [list of str (‘scale’, ‘auto’) or float, default: ['scale', 'auto']]: Kernel coefficient for 'rbf', 'poly' and 'sigmoid'; (See SVM below for more info).
        - **coef0** [list of float, default: [0.0]]: Independent term in kernel function. It is only significant in ‘poly’ and ‘sigmoid’.
        - **shrinking** [list of bool, default: [True]]: Whether to use the shrinking heuristic; (See SVM below for more info).
        - **tol** [list of float, default: [0.001]]: Tolerance for stopping criterion.
        - **decision_function_shape** [list of str (‘ovo’, ‘ovr’), default: ['ovr']]: Decision function shape ('ovr' for one-vs-rest, 'ovo' for one-vs-one); (See SVM below for more info).

    - **Cross Validation:**
      - **n_splits** [int, default: 5]: Number of folds in k-fold cross-validation.
      - **shuffle** [bool, default: True]: Whether to shuffle the data before splitting.

    - **SVM**:
      - **C** [float, default: 1.0]: Regularization parameter (inverse of regularization strength). Must be strictly positive. The penalty is a squared L2 penalty.
      - **kernel** [str (‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’), default: 'rbf']: Kernel type to be used in the algorithm.  For an intuitive visualization of different kernel types see [Plot classification boundaries with different SVM Kernels](https://scikit-learn.org/stable/auto_examples/svm/plot_svm_kernels.html#sphx-glr-auto-examples-svm-plot-svm-kernels-py).
        - 'linear': Linear kernel, equivalent to no transformation.
        - 'poly': Polynomial kernel, allows for curved decision boundaries.
        - 'rbf': Radial Basis Function kernel, creates circular/elliptical decision boundaries.
        - 'sigmoid': Sigmoid kernel, similar to a neural network.
      - **degree** [int, default: 3]: Polynomial degree for 'poly' kernel. Degree of the polynomial kernel function. Must be non-negative. Ignored by all other kernels.
      - **gamma** [str (‘scale’, ‘auto’) or float, default: 'scale']: Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.
         - 'scale' (default) is passed then it uses 1 / (n_features * X.var()) as value of gamma.
         - ‘auto’, uses 1 / n_features.
         - float, must be non-negative.
      - **coef0** [float, default: [0.0]]: Independent term in kernel function. It is only significant in ‘poly’ and ‘sigmoid’.
      - **shrinking** [bool, default: True]: Whether to use the shrinking heuristic. If the number of iterations is large, then shrinking can shorten the training time. However, if you loosely solve the optimization problem (e.g., by using a large stopping tolerance), the code without using shrinking may be much faster.
      - **probability** [bool, default: False]: Whether to enable probability estimates. This must be enabled prior to calling fit, will slow down that method as it internally uses 5-fold cross-validation, and predict_proba may be inconsistent with predict.
      - **tol** [float, default: 0.001]: Tolerance for stopping criterion.
      - **class_weight** [dict or 'balanced', default: None]: Weights associated with classes. Set the parameter C of class i to class_weight[i]*C for SVC. If not given, all classes are supposed to have weight one. The “balanced” mode uses the values of y to automatically adjust weights inversely proportional to class frequencies in the input data as n_samples / (n_classes * np.bincount(y)).
      - **max_iter** [int, default: -1]: Hard limit on iterations within solver, or -1 for no limit.
      - **decision_function_shape** [str (‘ovo’, ‘ovr’), default: 'ovr']: Decision function shape ('ovr' for one-vs-rest, 'ovo' for one-vs-one). Defines the shape of the decision function for multi-class classification. It can be 'ovr' (one-vs-rest) for a (n_samples, n_classes) shape or 'ovo' (one-vs-one) for the original (n_samples, n_classes * (n_classes - 1) / 2) shape. Internally, models are always trained using 'ovo', and 'ovr' is derived from it. This parameter is not applicable to binary classification.
      - **break_ties** [bool, default: False]: If true, decision_function_shape='ovr', and number of classes > 2, predict will break ties according to the confidence values of decision_function; otherwise the first class among the tied classes is returned. Please note that breaking ties comes at a relatively high computational cost compared to a simple predict.

##### Neural Decoding
###### [correlation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.pearsonr.html)
  - **Description:**  Calculates pairwise correlations between neural signals.
  - **Input:** FluoData (neural activity), IsCellData (optional)
   - **Neural data shape should be (ROI, timepoints)**
  - **Output:** HeatMapData (correlation matrix)
  - **Parameters:**
   - **transpose** [bool, default: False]: Whether to transpose the neural data matrix.

###### [cross_correlation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.correlate.html)
  - **Description:** Estimate the similarity between pairs of neural signals as a function of a specified time-lag.
  - **Input:** FluoData (neural activity), IsCellData (optional)
    - **Neural data shape should be (ROI, timepoints)**
  - **Output:** TimeSeriesData (cross-correlation values for each pair)
  - **Parameters:**
   - **transpose** [bool, default: False]: Whether to transpose the neural data matrix.
   - **lags** [int, default: 10]: Maximum number of time lags to compute correlations for.
   - **method** [str (‘auto’, ‘direct’, ‘fft’), default 'auto']: Method for correlation calculation.
     - direct: The correlation is determined directly from sums, the definition of correlation.
     - fft: The Fast Fourier Transform is used to perform the correlation more quickly.
     - auto: Automatically chooses direct or Fourier method based on an estimate of which is faster.
   - **shuffle_sample_number** [int, default: 100]: Determines how many times the data is shuffled to create the baseline distribution for testing.
   - **shuffle_confidence_interval** [float, default: 0.95]: determines the confidence level used when calculating the interval of this baseline distribution.

###### [Granger](https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.grangercausalitytests.html) (Granger causality test)
  - **Description:** Performs Granger causality tests to assess potential causal relationships between neural signals.
  - **Input:** FluoData (neural activity), IsCellData (optional)
    - **Neural data shape should be (ROI, timepoints)**
  - **Output:** HeatMapData (Granger F-values), ScatterData (Granger F-values)
  - **Parameters:**
  - **I/O:**
    - **transpose** [bool, default: False]: Whether to transpose the neural data matrix.
    - **standard_mean** [bool, default: True]: Whether to standardize X by subtracting mean.
    - **standard_std** [bool, default: True]: Whether to standardize X by dividing by std.

  - **Granger:**
    - **Granger_maxlag** [int, default: 1 (no lag)]: If an integer, computes the test for all lags up to maxlag. If an iterable int e.g.[1,3,5], computes the tests only for the lags in maxlag.
    - **Granger_addconst** [bool]: Whether to add a constant term to the model.
    - **use_adfuller_test** [bool, default: ]: Whether to perform Augmented Dickey-Fuller test for stationarity. Determines if the time series is suitable for Granger causality analysis. See[Augmented Dickey-Fuller test documentation](https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.adfuller.html)
    - **use_coint_test** [bool]: Whether to perform cointegration test. Checks if non-stationary variables are cointegrated, which can affect the interpretation of Granger causality. See[cointegration test documentation](https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.coint.html)

    - **adfuller**: The Augmented Dickey-Fuller test can be used to test for a unit root in a univariate process in the presence of serial correlation.
        - **maxlag** [int or None]: Maximum lag order for the test,  default value of 12*(nobs/100)^{1/4} is used when None.
        - **regression** [str ('c', 'ct','ctt','n'), default: 'c']: Type of regression to use.
            - “c” : constant only (default).
            - “ct” : constant and trend.
            - “ctt” : constant, and linear and quadratic trend.
            - “n” : no constant, no trend.
        - **autolag** [str (“AIC”, “BIC”, “t-stat”, None), default: 'AIC']: Method to use for automatically determining the lag order.
          - If “AIC” (default) or “BIC”, then the number of lags is chosen to minimize the corresponding information criterion.
          - “t-stat” based choice of maxlag. Starts with maxlag and drops a lag until the t-statistic on the last lag length is significant using a 5%-sized test.
          - If None, then the number of included lags is set to maxlag.
        - regresults [bool, default: False]: If True, returns the full regression results.
        - store ?
        - regresults ?

    - **coint**: The null hypothesis is no cointegration. Variables in y0 and y1 are assumed to be integrated of order 1, I(1).Uses the augmented Engle-Granger two-step cointegration test. Constant or trend is included in 1st stage regression, i.e. in cointegrating equation.
      - **trend** [str ('c', 'ct','ctt','n'), default: 'c']: The trend term included in regression for cointegrating equation.
            - “c” : constant only (default).
            - “ct” : constant and trend.
            - “ctt” : constant, and linear and quadratic trend.
            - “n” : no constant, no trend.
      - **maxlag:** Argument for adfuller, largest or given number of lags; see above.
      - **autolag:** Argument for adfuller, lag selection criterion; see above.

#### OptiNiSt (Utils)
###### Microscope to Image
  - **Description:** Data type conversion: Converts microscope data to image data.
  - **Input:** MicroscopeData
  - **Output:** ImageData
  - **Parameters:**
   - **ch:** [int, default: 0] Channel to extract from the microscope data.

###### Fluo from HDF5
  - **Description:** Extracts fluorescence data from HDF5 file and transposes for visualization
  - **Input:** FluoData, HDF5Data
    - If using optinist nwb format:
        - processing/ophys/suite2p_roi_UNIQUE-ID/Fluorescence/data
        - processing/ophys/caiman_cnmf_UNIQUE-ID/Fluorescence/data
  - **Output:** FluoData

###### ROI from HDF5
  - **Description:** Extracts ROI data from HDF5 file and prepares it for visualization
  - **Input:** ImageData, HDF5Data, IscellData
    - If using optinist nwb format:
      cell_img:
        - processing/optinist/suite2p_roi_UNIQUE-ID_all_roi_img/data
        - processing/optinist/caiman_cnmf_UNIQUE-ID_all_roi_img/data
      iscell:
        - processing/ophys/ImageSegmentation/suite2p_roi_UNIQUE-ID/iscell
        - processing/ophys/ImageSegmentation/caiman_cnmf_UNIQUE-ID/iscell
  - **Output:** IscellData
