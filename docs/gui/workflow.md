Workflow
=================
OptiNiSt makes it easy to create analysis pipelines using the GUI.

In the workflow field, you can:
- Select the data and the algorithms or analysis methods (node).
- Connect these nodes to define the order of processing (workflow).

The analysis pipeline can be parallel or bifurcating.

<br>
<p align="center">
<img width="600px" src="../_static/workflow/whole.png" alt="Workflow window" />
<br/>


## Creating workflow
You can create a new workflow by clicking the + button.

<br>
<p align="center">
  <img width="600px" src="../_static/workflow/create_new_workflow.png" alt="Create workflow" />
</p>
<br>

```{eval-rst}
.. caution::
   By creating workflow, current workflow will be deleted.

   The records are kept if you have already run the workflow. You can reproduce the workflow from RECORD tab. See details :ref:`here <ReproduceButton>`.
```


### Setting Input data
By default, an Image node is displayed. This node defines the path to the data to use.

<br>
<p align="center">
<img width="250px" src="../_static/workflow/image_node.png" alt="Image node" />
</p>
<br/>

To select the data, follow these steps:
1. Click on the SELECT IMAGE button on the Image node.
2. Select a file or a folder. Choosing a folder makes all the TIFF files in the shown sequence an input set of continuous frames.

You can select files in `input` directory under `OPTINIST_DIR`.
To put files there, see next [](DirectorySetting) section.

<br>
<p align="center">
<img width="400px" src="../_static/workflow/imageList.png" alt="Image file selection" />
</p>
<br/>

```{eval-rst}
.. note::
   Currently, only image files with {.tif, .TIF, .tiff, .TIFF} extensions are supported.
```

(DirectorySetting)=
#### Directory Setting
OptiNiSt uses `OPTINIST_DIR` for retrieving data and saving results. OptiNiSt searches for input data in the 'input' directory within `OPTINIST_DIR`. The default `OPTINIST_DIR` is `/tmp/studio` on your computer.

Choosing a folder makes all the TIFF files in the shown sequence an input set of continuous frames.

You may not want to modify your original data folder, or you may want to make your data folder visible and accessible to OptiNiSt because imaging data can be large and take time to copy. You can take either strategy in assigning your data path:

1. **Upload from GUI**

    Click on the LOAD button on the node. The LOAD button copies the selected file to your `OPTINIST_DIR/input`.

    **By this method, you cannot upload multiple files or folder at once**.
    - If you want to upload multiple files or folder at once, use the method below.

2. **Copy files to `OPTINIST_DIR`**

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

3. **Change the setting of `OPTINIST_DIR`**

    This requires modifying source codes. See [](each-platforms-for-developer) installation guide section.
    `OPTINIST_DIR` is defined in `optinist/studio/app/dir_path.py`. Change line for `OPTINIST_DIR`, `INPUT_DIR`, and `OUTPUT_DIR` according to your demand. Changing `dir_path.py` may also be necessary when running the pipeline on your cluster computers. Also, you can quickly change `OPTINIST_DIR` by changing the environment variable before launching. The change is effective after relaunching.

#### Other Data Formats As The Input

<br>
<p align="center">
<img width="300px" src="../_static/workflow/csv_connect.png" alt="CSV Connect" />
</p>

eta, cca, correlation, cross_correlation, granger, glm, lda, and svm assume the input neural data shape is frames x cells matrix.
Because the output of CaImAn and Suite2P on the pipeline is cell x frames, the default setting for neural data for these analyses is set to transpose.
Pca and tsne can be done in either direction depending on your purpose.
The function assumes their input to be samples x features.
Rows and columns can be specified by `settings` appearing after selecting the csv data.
Note that the number of data points has to be the same as the number of frames of image data.
Fluo data node is for cell's fluorescence timecourse data given as .csv.

Another data format prepared is hdf5. This format is compatible with the nwb data format.
CSV and hdf5 nodes have black output connectors.
The edge connected to the black output connector can be connected to any input connector. Be careful; this means that it does not check the format correspondence between input and output.


### Adding Nodes
Select Data or Algorithm node from the treeview on the left.
Clicking the + button adds the analysis nodes to the Workflow field.

<br>
<p align="center">
  <img width="600px" src="../_static/workflow/add_algorithm.png" alt="Add Algorithm" />
</p>
<br>

The left side of the window displays all available analysis methods. ROI detection tools (currently Suite2P, CaImAn and LCCD) are in the "Algorithm" category, and all other pre-installed analyses are in the "optinist" category.

### Algorithm Nodes

Each algorithm node has PARAM button and OUTPUT button.

<br>
<p align="center">
<img width="300px" src="../_static/workflow/algorithm_node_buttons.png" alt="Algo node buttons" />
</p>
<br/>

- **PARAM**

  You can see or edit parameters for the algorithm.
  Clicking param button, the parameters are displayed in the right side of the window.

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


## Running pipelines
Click the RUN button at the top right.
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

## SNAKEMANE settings
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


## NWB setting
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
   So, they are not used for analysis.

   But, only ``imaging_plane.imaging_rate`` is used as parameter for frame rate.
   See details in :ref:`Switch time course plot units <SwitchTimeUnit>`.
```
