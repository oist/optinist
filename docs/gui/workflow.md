Workflow
=================

<br>
<p align="center">
<img width="400px" src="../_static/workflow/whole.png" alt="workflow" />
<br/>
  
  - OptiNiSt makes it easy to create analysis pipelines using the GUI.
- In the workflow field, you can:
    - Select the data and the algorithms or analysis methods (nodes).
    - Connect these nodes to define the order of processing (pipelines).
- The analysis pipeline can be parallel or bifurcating.

## Creating workflow
### Setting Input data
  
  
<br>
<p align="center">
<img width="200px" src="../_static/tutorials/fig3_imagenode.png" alt="imageNode" />
<img width="200px" src="../_static/workflow/components/imageList.png" alt="imageNode" />
</p>
<br/>

By default, an Image node is displayed. This node defines the path to the data to use.

Once the data is accessible, you can view it by following these steps:

1. Click on the SELECT IMAGE button on the Image node.
2. Select a file or a folder. Choosing a folder makes all the TIFF files in the shown sequence an input set of continuous frames.

  
**Note:** Currently, image files with {.tif, .TIF, .tiff, .TIFF} extensions are accepted. Other extensions will be added on request.

#### Directory Setting
OptiNiSt uses `OPTINIST_DIR` for retrieving data and saving results. OptiNiSt searches for input data in the 'input' directory within `OPTINIST_DIR`. The default `OPTINIST_DIR` is `/tmp/optinist` on your computer.

Choosing a folder makes all the TIFF files in the shown sequence an input set of continuous frames.

You may not want to modify your original data folder, or you may want to make your data folder visible and accessible to OptiNiSt because imaging data can be large and take time to copy. You can take either strategy in assigning your data path:

1. **Copy your original data file to `OPTINIST_DIR` ** To copy the data to `OPTINIST_DIR`, click on the LOAD button on the node. The LOAD button copies the selected file to your `OPTINIST_DIR/input`. This can be done from the GUI.

2. **Change the setting of `OPTINIST_DIR` ** `OPTINIST_DIR` is defined in `optinist/optinist/api/dir_path.py`. Change line for `OPTINIST_DIR`, INPUT_DIR, and OUTPUT_DIR according to your demand. Changing `dir_path.py` may also be necessary when running the pipeline on your cluster computers. Also, you can quickly change OPTINIST_DIR by changing the environment variable before launching. The change is effective after relaunching.

#### Other Data Formats As The Input

<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/csv_connect.png" alt="CSV Connect" />
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
<br>
<p align="center">
  <img width="300px" src="../_static/workflow/components/add_algorithm.png" alt="Add Algorithm" />
</p>
<br>

Select algorithms or analysis methods from the treeview on the left by clicking "+" button. 

The left side of the window displays all available analysis methods. Clicking on the + mark adds the analysis nodes to the Workflow field. ROI detection tools (currently Suite2P, CaImAn and LCCD) are in the "Algorithm" category, and all other pre-installed analyses are in the "optinist" category.

Let's start with sample TIFF data (`mouse2p_2_long.tiff`) and try Suite2P ROI detection. 
First, you need to determine the image you will use. Select your image as explained [above](#assigning-input-data-path). 
Once it is selected, the name of the files is shown in the Image node.

### Parameter button and output button on the node

<br>
<p align="center">
<img width="300px" src="../_static/tutorials/fig10.1_buttons.png" alt="Set Parameter" />
</p>
<br/>

Each node has PARAM button and OUTPUT button. 

- **Editing Parameters:** Click on the PARAM button to view the parameters. You can edit them as needed. The names, types, and default values of the parameters are the same as the original algorithms. Refer to the original documentation to confirm the meaning of the parameters. The link list is available at [Implemented Analysis](https://optinist.readthedocs.io/en/latest/utils/implemented_analysis.html).

- **Checking Results:** The OUTPUT button is for a quick check of the results. The button becomes active after the successful execution of the pipeline. For details about the charts, see [Inspecting the Images and the Plots on Visualize](#inspecting-the-images-and-the-plots-on-visualize).

### Connecting Nodes 

<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/connect_edge.png" alt="Connect Algorithm" />
</p>
<br/>

Connect colored connectors of the nodes by dragging your cursor from the output connector to the next input connector to create connecting edges. The color of the connector indicates the data type of the input and the output.
You can only connect the input and output connectors of the same color. 

**DataType List**
- <span style="color: red; ">ImageData</span>
- <span style="color: green; ">Suite2pData</span>
- <span style="color: orange; ">Fluorescence</span>
- <span style="color: yellow; ">Behavior</span>
- <span style="color: blue; ">Iscell</span>


### Removing Nodes or Connects
Clicking on the x mark on a node or on an edge removes it from the workflow field. 


## Running pipelines
<br>
<p align="left">
<img width="200px" src="../_static/tutorials/fig11_runall.png" alt="Whole" />
</p>
<br/>

Click the RUN button at the top right to see two dropdown choices: RUNALL and RUN. 

- **RUNALL:**
    - Runs the entire process.
    - Assigns a new folder for saving the results.　This folder name is only for the user’s convenience. The actual folder name is long digit random letter+number. <!--Further information about the structure of the saved results is [here](https://optinist.readthedocs.io/en/latest/gui/record.html).)
-->

- **RUN:**
    - Skips already executed processes.
    - Checks for differences from the previous pipeline, including the existence of results and parameter changes.
    - If differences are detected, the downstream process is re-executed.
    - Results are overwritten in the existing folder.

- **Cancel:**
    - Abort the running pipeline immediately.

## SNAKEMANE and NWB SETTING

<br>
<p align="left">
<img width="400px" src="../_static/tutorials/fig13_nwbsnakemake.png" alt="Whole" />
</p>

SNAKEMAKE and NWB SETTING buttons are for parameters for snakemake and output NWB file.

- **SNAKEMAKE SETTING:**
    - OptiNiSt's pipeline construction is based on [Snakemake](https://snakemake.readthedocs.io/en/stable/), a pipeline controlling tool for Python scripts.
    - The Snakemake parameter setting is following.
      - use_conda: If this is on, snakemake uses conda environment.
      - cores: Specifies the number of cores to use. If not specified, snakemake uses number of available cores in the machine. 
      - forceall: Flag to indicate the execution of the target regardless of already created output.
      - forcetargets: Users may not want to change this. 
      - lock: Users may not want to change this.
    - For details about snakemake parameters please refer to [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html)

- **NWB SETTING:**
    - Defines the metadata associated with the input data.
    - By configuring this, the output NWB file includes the information set here.
    - The parameter you set here is only for your record and is not used for calculations within OptiNiSt.
    - You can leave this setting as the default.
    - The details of NWB setting in OptiNiSt are folling.
      - session_description: a description of the session where this data was generated
      - identifier: a unique text identifier for the file
      - experiment_description: general description of the experiment
      - device: device used to aquire the data (information such as manufacturer, firmware version, model etc.)
      - optical_channel: information about the optical channel used to acquire the data
      - imaging_plane: information about imaging such as sampling rate, excitation wave length, calcium indicator.
      - image_serises: information about imaing time
      - ophys: general information about imaging
    - For general information about NWB, refer to the [official documentation](https://www.nwb.org/getting-started/).

