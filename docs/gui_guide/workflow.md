# Workflow

OptiNiSt makes it easy to create analysis pipelines using the GUI.

In the workflow field, you can:

- Select the data and the algorithms or analysis methods (node).
- Connect these nodes to define the order of processing (workflow).

The analysis pipeline can be parallel or bifurcating.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/whole.png" alt="Workflow window" />
<br/>

## Data Nodes

See [](data_nodes) for a description of which data types each node accepts.

- Data shape is displayed in the file select dialog. Please check data shape if you have unexpected results.
- If you replace the image file with the same file name, shape cannot be updated automatically. Please click reload icon besides the checkbox.

## Algorithm Nodes

See [](algorithm_nodes) for a description of each data processing and analysis node.

## Create a Workflow

You can create a new workflow by clicking the + button.

<br>
<p align="center">
  <img width="600px" src="../_static/workflow/create_new_workflow.png" alt="Create workflow" />
</p>
<br>

```{eval-rst}
.. caution::
   By creating a new workflow, current workflow will be deleted.

   The records are kept if you have already RUN the workflow. You can reproduce the workflow from RECORD tab. See details :ref:`here <ReproduceButton>`.
```

### Adding Nodes

Select or drag Data or Algorithm nodes from the treeview on the left.
Clicking the + button adds the analysis nodes to the Workflow field.

<br>
<p align="center">
  <img width="600px" src="../_static/workflow/add_algorithm.png" alt="Add Algorithm" />
</p>
<br>

The left side of the window displays all available analysis methods. ROI detection tools (currently Suite2P, CaImAn and LCCD) are in the "Algorithm" category, and all other pre-installed analyses are in the "optinist" category.

(AddData)=

### Input Data with Data Nodes

- Image (TIFF)

  By default, an Image node is displayed. This node defines the path to the data to use.

  <p align="center">
  <img width="250px" src="../_static/workflow/image_node.png" alt="Image node" />
  </p>
  <br/>

  To select the data, follow these steps:

  1. Click on the check list icon on the Image node.
  2. Select a file or a folder. Choosing a folder makes all the TIFF files in the shown sequence an input set of continuous frames.

  You can select files in `input` directory under `OPTINIST_DIR`.
  To add files to the directory, see [](DirectorySettings) section.

  <p align="center">
  <img width="400px" src="../_static/workflow/imageList.png" alt="Image file selection" />
  </p>
  <br/>

- CSV, FLUO, BEHAVIOR

  These nodes are for csv data.
  Once the file selected, you can preview and change settings for the data.

  <p align="center">
  <img width="300px" src="../_static/workflow/csv_settings.png" alt="CSV settings" />
  </p>

  - transpose

    You can use transposed the data as the input by checking this box.

    <p align="center">
    <img width="400px" src="../_static/workflow/csv_transpose.png" alt="CSV transpose" />
    </p>

    ETA, CCA, granger, GLM, LDA, and SVM assume the input neural data shape is frames x cells matrix.
    Because the output of CaImAn and Suite2P on the pipeline is cell x frames, the default setting for neural data for these analyses is set to transpose.

    correlation and cross_correlation take cells x frames directly, so their transpose option is off by default.

    PCA and TSNE can be done in either direction depending on your purpose.
    The function assumes their input to be samples x features.

  - set_header

    If your CSV data has header, check this box.
    Set the header index to your data's header row. (**first row is 0**)
    The data below the header row will be used as the data.

    <p align="center">
    <img width="400px" src="../_static/workflow/csv_set_header.png" alt="CSV Set Header" />
    </p>

  - set_index

    You can add index column to the data by checking this box.

    <p align="center">
    <img width="400px" src="../_static/workflow/csv_set_index.png" alt="CSV Set Index" />
    </p>

- HDF5

  This node accepts `.hdf5` and `.nwb` file types. NWB is compatible with the HDF5 data format. You can use you OptiNiSt data analysis results, which save in NWB format, as the input.

  <p align="center">
  <img width="300px" src="../_static/workflow/hdf5_node.png" alt="HDF5 node" />
  </p>

  Once a file is selected, you can search through the data structure and select the data to use.

  <p align="center">
  <img width="400px" src="../_static/workflow/hdf5_structure.png" alt="HDF5 node" />
  </p>

There are some additional utility nodes for processing HDF5 data.
For example, fluo_from_hdf5 extracts fluorescence data from HDF5 data or
roi_fluo_from_hdf5 that extracts roi and fluorescence data from HDF5 data.

For more details, refer to the [Algorithm Nodes Documentation](../specifications/algorithm_nodes.md#fluo-from-hdf5).

  <p align="center">
  <img width="650px" src="../_static/workflow/hdf5_to_fluo.png" alt="HDF5 to fluo" />
  </p>

- Matlab

  This node accepts `.mat` files.
  Once a file is selected, you can search through the data structure and select the data to use.

  ```{eval-rst}
  .. note::
    CSV, HDF5 and Matlab nodes have black output connectors.
    Black output connector can be connected to any input connector.
    Be careful; this means that it does not check the format correspondence between input and output.
  ```

- Microscope

  See [](data_nodes) for more information Microscope nodes. Currently, the Microscope node can accept following data formats.

  - Inscopix(.isxd)
  - NIKON(.nd2)
  - Olympus(.oir)

    <p align="center">
  <img width="650px" src="../_static/workflow/microscope.png" alt="Microscope" />
  </p>

  You can convert these data into ImageData so that you can use them in the pipeline.
  To convert, connect microscope_to_img Algorithm node to the microscope data node.

### Algorithm Node Settings

Each algorithm node has PARAM button and OUTPUT button.

<br>
<p align="center">
<img width="300px" src="../_static/workflow/algorithm_node_buttons.png" alt="Algo node buttons" />
</p>
<br/>

- **PARAM**

  You can see or edit parameters for the algorithm.
  Clicking `PARAM` button, the parameters are displayed in the right side of the window.
  More information on each parameter can be found by clicking the external link icon which links to our [Algorithm Specifications](algorithm_nodes) page.

  <br>
  <p align="center">
  <img width="600px" src="../_static/workflow/algorithm_param_form.png" alt="Algo node parameter form" />
  </p>
  <br/>

  The names, types, and default values of the parameters are the same as the original algorithms. Refer to the original documentation to confirm the meaning of the parameters. The link list is available at [Implemented Analysis](https://optinist.readthedocs.io/en/latest/utils/implemented_analysis.html).

- **OUTPUT**

  You can check the results of algorithm quickly after the successful execution of the pipeline. For details about the charts, see [](Visualize).

  <br>
  <p align="center">
  <img width="400px" src="../_static/workflow/algorithm_output.png" alt="Algo node output viewer" />
  </p>
  <br/>

(ConnectingNodes)=

### Connecting Nodes

Connect colored connectors of the nodes by dragging your cursor from the output connector(right side of the nodes) to the next input connector(left side of the nodes) to create connecting edges.

<br>
<p align="center">
<img width="500px" src="../_static/workflow/connect_edge.png" alt="Connect nodes" />
</p>
<br/>

The color of the connector indicates the data type of the input and the output.
You can only connect the input and output connectors of the same color.

**DataType List**

- <span style="color: red; ">ImageData</span>
- <span style="color: green; ">Suite2pData</span>
- <span style="color: orange; ">Fluorescence</span>
- <span style="color: #cfc22b; ">Behavior</span>
- <span style="color: blue; ">Iscell</span>
- <span style="color: black; ">Unspecified</span>

### Removing Nodes or Connects

Clicking on the x mark on a node or on an edge removes it from the workflow field.

(ImportWorkflowYaml)=

### Import existing workflow by yaml file

You can create same workflow by importing workflow config yaml format file. This feature is useful to share workflow template.

Workflow config file can be downloaded from RECORD tab's executed workflow.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/download_workflow_config.png" alt="Download workflow config file" />
</p>
<br/>

Clicking the import icon button and select shared workflow config yaml.
Then you can create same workflow as the file's record.
Input file selection will not reproduced because it may not be in your device.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/import_workflow.png" alt="Import workflow" />
</p>
<br/>

## Running Pipelines

Click the RUN button at the top right.

Note that while the workflow running, you **cannot**

- use [](Visualize) tab
- delete current workflow from [](Record) tab
- reproduce other workflow from [](Record) tab

There are 2 types of execution. You can select the type by clicking the dropdown button.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/run_buttons.png" alt="Run buttons" />
</p>
<br/>

### RUN ALL

- Runs the entire process.
- Assigns a new folder for saving the results. This folder name is only for the user’s convenience. The actual folder name is long digit random letter+number.

### RUN

- Available only when the pipeline has been executed.
- Useful to re-run the workflow with different parameters.
- By checking parameter changes and addition of new nodes, it would skip already executed processes.

  - e.g. Let's say you have executed following workflow.

    <br>
    <p align="left">
    <img src="../_static/workflow/run_before.png" alt="Run example before" />
    </p>
    <br/>

    Then you have changed the parameter "block_size" of suite2p_registration from 128,128 to 256,256.
    With "RUN", only following nodes would be executed again.

    - suite2p_registration: because you have changed parameter for this node
    - suite2p_roi: because its upstream data (from suite2p_registration) would be changed.

    Nodes which would be executed by "RUN" are highlighted in yellow.

    <br>
    <p align="left">
    <img src="../_static/workflow/run_after.png" alt="Run example after" />
    </p>
    <br/>

  - Following changes would affect whole workflow. So you cannot use "RUN" button after changing them.
    - Input data
    - NWB settings

```{eval-rst}
.. caution::
   With RUN, results will be overwritten. To avoid this, use RUN ALL.
```

```{eval-rst}
.. note::
   If you want to set up the Conda Environment for each node first, check :ref:`Setup Conda Environment <setup_conda_environment>`
```

### Filtering Data

#### Overview

Filtering data is only available after executing the pipeline. This feature is particularly useful when you need to examine a specific range of the Region of Interest (ROI) or fluorescence time series.

#### Step-by-Step Instructions

##### 1. Ensure Pipeline Execution

Before applying filters, make sure that you have executed the workflow. For example, if you followed **Tutorial 1**, your workflow might look like this:

<br>
<p align="left">
<img src="../_static/workflow/filter_before.png" alt="Workflow before filtering" />
</p>
<br/>

##### 2. Selecting the Filter Option

The filter placeholder values tell you the accepted range, based on the current data. Values can be added as comma-seaparated single integers or as a range, for example: 0, 1, 4:10.

If you want to filter ROI data within a specific range (e.g., 0 to 100), follow these steps:

1. Locate the `suite2p_roi` node in your workflow.
2. Click on the **Filter** button in the `suite2p_roi` node.

<br>
<p align="center">
<img width="300px" src="../_static/workflow/filter_button.png" alt="Filter button" />
</p>
<br/>

##### 3. Applying the Filter

1. In the filter settings, set the **ROI Data Range** from **0 to 100**.
2. Click the **Apply** button to confirm the filter settings.

<br>
<p align="center">
<img width="400px" src="../_static/workflow/filter_box.png" alt="Filter settings" />
</p>
<br/>

##### 4. Running the Workflow

Once you apply the filter you must `RUN` in order for the filter to be applied downstream.

- The `suite2p_roi` node will start **loading**.
- When loading completes, the `eta` node will be highlighted in **yellow**.
- Click the **RUN** button to execute the workflow and apply the filter changes.

##### 5. Filter output
Filtered data will be saved to the output NWB file once `RUN` has been performed. Filtered data will be saved in the same folders as the ROI detection algorithm, and additionally filter parameters will be saved in the optinist directory of the NWB. For example for CaImAn:
- ophys
  - ImageSegmentation
    - **filtered_caiman_cnmf_xxxx (id...iscell)**
    - caiman_cnmf_xxxx (id...iscell)
  - **filtered_caiman_cnmf_xxxx**
    - Fluorescence (data, rois, ...)
  - caiman_cnmf_xxxx
    - Fluorescence (data, rois, ...)
- optinist
  - **filtered_caiman_cnmf_xxxx_filter_roi_ind**
  - **filtered_caiman_cnmf_xxxx_filter_time_ind**

(DirectorySettings)=
## Directory Settings

OptiNiSt uses `OPTINIST_DIR` for retrieving data and saving results. OptiNiSt searches for input data in the 'input' directory within `OPTINIST_DIR`. The default `OPTINIST_DIR` is `/tmp/studio` on your computer.

Choosing a folder makes all the TIFF files in the shown sequence an input set of continuous frames.

You may not want to modify your original data folder, or you may want to make your data folder visible and accessible to OptiNiSt because imaging data can be large and take time to copy. You can take either strategy in assigning your data path:

1. **Upload from GUI**

   Click on the image icon on the node. The LOAD button copies the selected file to your `OPTINIST_DIR/input`.

   <p align="center">
   <img width="300px" src="../_static/workflow/image_node_upload_from_gui.png" alt="upload file" />
   </p>

   **By this method, you cannot upload multiple files or folder at once**.

   - If you want to upload multiple files or folder at once, use the next method.

2. **Copy Files to `OPTINIST_DIR`**

   Copy your raw data to `OPTINIST_DIR/input/1/` by your OS's file manager or command lines.

   ```{eval-rst}
   .. warning::
       Be sure to copy under ``input/1/``. ``1`` is the default workspace id for :ref:`standalone mode <about-multiuser-mode>`.
       If you copy under ``input/`` directly, the file cannot be found from GUI.
   ```

   You can copy folder into the input directory.

   - If you put folder, you can see the folder from GUI, SELECT IMAGE dialog like this.
     <br>
     <p align="center">
     <img width="400px" src="../_static/workflow/put_folder_to_input_dir.png" alt="Put folder to input dir" />
     </p>

3. **Get file via URL**

   Click on the link icon on the node.

   <p align="center">
   <img width="300px" src="../_static/workflow/image_node_from_cloud.png" alt="Get File from Cloud" />
   </p>

   Then input dialog will be shown. You can input the URL of the file you want to use.

   <p align="center">
   <img width="450px" src="../_static/workflow/image_node_from_cloud_input.png" alt="Input FILE URL" />
   </p>

   ```{eval-rst}
   .. note::
       The URL should be direct link to the one file.
       It should be

         - started with ``http://`` or ``https://``.
         - ended with the file name with extension.

       And more, download links require authentication are not supported.
   ```

4. **Change the Setting of `OPTINIST_DIR`**

   This requires modifying source codes. See [](native_platforms_developer) installation guide section.
   `OPTINIST_DIR` is defined in `optinist/studio/app/dir_path.py`. Change line for `OPTINIST_DIR`, `INPUT_DIR`, and `OUTPUT_DIR` according to your demand. Changing `dir_path.py` may also be necessary when running the pipeline on your cluster computers. Also, you can quickly change `OPTINIST_DIR` by changing the environment variable before launching. The change is effective after relaunching.

## SNAKEMAKE Settings

OptiNiSt's pipeline construction is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/), a pipeline controlling tool for Python scripts.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/snakemake_settings.png" alt="Snakemake settings" />
</p>
<br/>

The Snakemake parameter setting is following.

- use_conda: If this is on, snakemake uses conda environment.
- cores: Specifies the number of cores to use. If not specified, snakemake uses number of available cores in the machine.
- forceall: Flag to indicate the execution of the target regardless of already created output.
- forcetargets: Users may not want to change this.
- lock: Users may not want to change this.

For details about snakemake parameters please refer to [snakemake official page](https://snakemake.readthedocs.io/en/stable/executing/cli.html).

## NWB Settings

Defines the metadata associated with the input data.
By configuring this, the output NWB file includes the information set here.
You can leave this setting as the default.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/nwb_settings.png" alt="NWB settings" />
</p>
<br/>

The details of NWB setting in OptiNiSt are following.

- session_description: a description of the session where this data was generated
- identifier: a unique text identifier for the file
- experiment_description: general description of the experiment
- device: device used to acquire the data (information such as manufacturer, firmware version, model etc.)
- optical_channel: information about the optical channel used to acquire the data
- imaging_plane: information about imaging such as sampling rate, excitation wave length, calcium indicator.
- image_serises: information about imaing time
- ophys: general information about imaging

For general information about NWB, refer to [NWB official page](https://www.nwb.org/getting-started/).

```{eval-rst}
.. note::
   Basically, these parameters are description for your experiment.
   So, they are not used for analysis except for following parameters.

   - ``imaging_plane.imaging_rate``
      This will be used as parameter for frame rate.
      See details in :ref:`Switch time course plot units <SwitchTimeUnit>`.

   - ``image_series.save_raw_image_to_nwb``
      If True, raw image data will be saved to NWB file's acquisition.
      If not, only the path to the image data will be saved.
```

## Viewing Logs with Log Viewer

When running a workflow, logs are continuously generated to provide insights into the execution process. These logs can be monitored and analyzed using the **Log Viewer**.

To access the logs, simply click on the **Log Button** in the interface.

<br>
<p align="center">
<img width="300px" src="../_static/workflow/log_button.png" alt="Algo node buttons" />
</p>
<br/>

### Key Features of Log Viewer

The Log Viewer provides powerful capabilities for log analysis, including:

- **Real-time Log Monitoring** – View logs as they are generated while the workflow runs.
- **Highlight & Copy** – Easily highlight and copy log entries for further investigation or sharing.
- **Keyword Search** – Search for specific words or phrases to quickly locate relevant log entries.
- **Log Level Filtering** – Filter logs based on their severity level to focus on critical events.

### Supported Log Levels

Logs are categorized into different levels based on severity:

| Log Level    | Description                                                           |
| ------------ | --------------------------------------------------------------------- |
| **Debug**    | Detailed information useful for troubleshooting.                      |
| **Info**     | General information about the workflow execution.                     |
| **Warning**  | Potential issues that may not immediately impact execution.           |
| **Error**    | Issues that have affected the workflow but are not necessarily fatal. |
| **Critical** | Severe errors that require immediate attention.                       |

### Example Log Output

  <br>
  <p align="left">
  <img src="../_static/workflow/log_viewer_sample.png" alt="Log Viewer Sample" />
  </p>
  <br/>
