Tutorials
=================

* [Opening the browser](#opening-the-browser)
* [Making pipelines on WORKFLOW](#making-pipelines-on-workflow)
  * [Assigning input data path](#assigning-input-data-path)
  * [Selecting analysis methods](#selecting-analysis-methods)
  * [Creating pipelines](#creating-pipelines)
  * [Parameter button and output button on the node](#parameter-button-and-output-button-on-the-node)
  * [running pipelines](#running-pipelines)
  * [SNAKEMANE and NWB SETTING](#snakemane-and-nwb-setting)
  * [Time series analyses after ROI extraction](#time-series-analyses-after-roi-extraction)
  * [Additional information on WORKFLOW](#additional-information-on-workflow)
        * [setting OPTINIST_DIR](#setting-optinist_dir)
      * [about the assumed data shape](#about-the-assumed-data-shape)
        * [snakemake settings](#snakemake-settings)
        * [NWB settings](#nwb-settings)
* [Inspecting the images and the plots on VISUALIZE](#inspecting-the-images-and-the-plots-on-visualize)
  * [Checking movies](#checking-movies)
  * [Showing ROI and time courses](#showing-roi-and-time-courses)
  * [Savind plots](#savind-plots)
* [Managing pipelines on RECORD](#managing-pipelines-on-record)

## Opening the browser
To start OptiNiSt, you need to open a console and activate optinist environment `conda activate optinist` and type `run_optinist` or, change to optinist directory `cd ~/optinist/` and run main script `python main.py`. 
The console shows the log once the startup is completed.

<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig1_console.png" alt="Whole" />
</p>

Once you see this, open your web browser (Google Chrome is recommended) at localhost:8000.
You are ready to start if the OptiNiSt page appears.

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig2_open.png" alt="Whole" />
</p>

OptiNiSt has three different pages, WORKFLOW, VISUALIZE, and RECORD. You can toggle these by clicking on the tag. 
<br>
<p align="left">
<img width="300px" src="./_static/tutorials/fig2.2_tags.png" alt="Whole" />
</p>

## Making pipelines on WORKFLOW 
After launching, the first page you see is the workflow page. The workflow page is a place to define the analysis pipeline. You determine the data you will analyze, select the type of the algorithm or analysis method you use, and set the parameters and the order of analysis.  

### Assigning input data path
As a default, it shows an image node. This node defines the path to the data to use.  

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig3_imagenode.png" alt="Whole" />
</p>

OptiNiSt uses OPTINIST_DIR for retrieving data and saving results. OptiNiSt searches input data from the 'input' directory in OPTINIST_DIR. A default OPTINIST_DIR is `/tmp/optinist` in your computer.

You may not want to change anything in your original data folder, or you may wish to make your data folder visible and accessible to OptiNist because the imaging data is sometimes huge and takes time to copy. You can take either strategy in assigning your data path.

1. Copy your original data file to OPTINIST_DIR and assign the data path to the copied data.   
Clicking on the UPLOAD button on the node opens the file explorer or finder so that you can select the data file. UPLOAD button copies the selected file to your OPTINIST_DIR/input. This can be done from the GUI.  
  
2. Change the setting of OPTINIST_DIR by editing dir_path.py file. See [setting optinist directory](#setting-optinist_dir). Change is effective after re-launching.

Once the data is made accessible, the SELECT IMAGE button on the image node can assign the file as the input to the pipeline. You can select a file or a folder. Choosing a folder makes all the tiff files in the shown sequence an input set of continuous frames. See the case input file is [other than tiff] (../docs/gui/workflow.md)/

### Selecting analysis methods
The left side of the window shows all available analysis methods. Clicking on the + mark adds the analysis nodes to the workflow field. ROI detection tools (currently suite2P and CaImAn ) are in ‘Algorithm’ category, and all other pre-installed analyses are in ‘optinist’ category.

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig4_algorithms.png" alt="Whole" />
</p>

Let’s start with sample tiff data (mouse2p_2_long.tiff) and try suite2p ROI detection.
First, you need to determine the image you use. Select your image as explained [above](#assigning-input-data-path).
Once it is selected, it shows the name of files in the image node. 

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig5_imagenode2.png" alt="Whole" />
</p>

### Creating pipelines

Then, connect the analysis nodes in the order you like to process. Drugging from an output connector to an input connector creates an edge. The color of the connector indicates the format. For example, red is the image type format. You can only connect the same color. (Exception: black is an undefined data format. You can connect the black connector with any other connector, but be careful it does not check the consistency of input and output).

<br>
<p align="left">
<img width="800px" src="./_static/tutorials/fig6_suite2p.png" alt="Whole" />
</p>

As for Suite2P, you might not use suite2P_registration (motion correction). In that case, you can connect the suite2p_file_convert to suite2p_roi directly. 

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig7_suite2p2.png" alt="Whole" />
</p>

Also, you can perform motion correction of CaImAn (caiman_mc) and then perform suite2P_roi.

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig8_suite2pcaiman.png" alt="Whole" />
</p>

The nodes should be connected as long as the input and the output are of the same format type (same color).
Also, you can branch the flow. In the example, the two caiman_mc with different parameter settings are created, and the downstream from caiman_mc is also different. Each node's results are saved in a separate folder (See [RECORD part](#managing-pipelines-on-record)). 

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig9_workflowbranch.png" alt="Whole" />
</p>

### Parameter button and output button on the node
Each node has PARAM button and OUTPUT button. 

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig10.1_buttons.png" alt="Whole" />
</p>

Clicking on PARAM shows the parameters. Edit this as you like. The names, types and the default values of the parameters are the same as the original algorithms. Refer to the original documentation to confirm the meaning of the parameters. The link list is on [README](../../README.md).

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig10_parameters.png" alt="Whole" />
</p>

OUTPUT button is for the quick check of the result. The button becomes effective after the successful execution of the pipeline. [Here](#inspecting-the-images-and-the-plots-on-visualize) explains the details of the charts. 

### running pipelines

Now you are ready to run the pipeline.
RUN button at the right top shows two pulldown choices. RUNALL runs all the process. RUNALL assigns a new folder for saving the results. On the other hand, RUN skips the already ran processes. It checks the difference from the previous pipeline, including the existence of the results and the parameter change. If they are detected, the downstream process after that process is re-executed. The results are overwritten into the existing folder.

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig11_runall.png" alt="Whole" />
</p>

When you click on the RUNALL, it shows the window to determine the folder name. This folder name is only for the user’s convenience. The actual folder name is long digit random letter+number. Further information about the structure of the saved results is [here](LINKTO /docs/gui/record.md #setting-input-images).


Next to the RUN button, there is the CANCEL button. You can abort the running pipeline with this button. It immediately cancels the current execution.

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12_cancel.png" alt="Whole" />
</p>



### SNAKEMANE and NWB SETTING

SNAKEMAKE and NWB SETTING buttons are for parameters for snakemake and output NWB file.
The pipeline construction of Optinist is based on snakemake (ref), which is the pipeline controlling tool for python scripts. The SNAKEMAKE parameter setting is [here](#snakemake-settings).

<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig13_nwbsnakemake.png" alt="Whole" />
</p>

NWB SETTING defines the metadata for the NWB file as an output. The parameter you set here is only for your record and not used for the calculation inside OptiNiSt. You can leave this as default. The details of NWB setting in OptiNiSt is [here](#nwb-settings). Also, general info about NWB is [here](https://www.nwb.org/getting-started/).




### Time series analyses after ROI extraction
OptiNiSt offers some basic time-series analysis functions. For example, event-triggered averaging can be applied to the ROI time-series data created by OptiNiSt. Assuming that you have the result of ROI extraction, here explains how to create the pipeline. Because the ROI time-series is in NWB format, the hdf5 data node is appropriate as the input node. 

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12.00_hdfnode.png" alt="Whole" />
</p>

Add the hdf5 node to the field. Upload the data to the OPTINIST_DIR. In addition to UPLOAD and SELECT to assign the file, you need to indicate the position of the fluorescence data in the HDF5 structure (STRUCTURE button appeared after you SELECT HDF5).

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12.01_hdfnode2.png" alt="Whole" />
</p>

NWB structure of Suite2P and CaImAn is different because OptiNiSt inherits each algorithm's original NWB output format. You will find the colums and rows are opposite between Suite2P outputs and CaImAn outputs. You can re-assign the rows and columns in the parameter setting of the analysis node (in this case eta node). 

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12.02_hdfs2p.png" alt="Whole" />
</p>

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12.03_hdfcaiman.png" alt="Whole" />
</p>

In this example, the behavioral data format is .csv. The csv data node or behavior data node is used for behavior input. 

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12.04_behaviornode.png" alt="Whole" />
</p>

Once you SELECT CSV, SETTINGS button appears in the behavior node. This button confirms the inside of csv data and makes it possible to transpose the matrix if needed. If your csv includes the headers, you can also assign it to ignore it in creating the matrix. Set Index adds index columns to the matrix.

<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig12.05_behaviornode2.png" alt="Whole" />
</p>


Add event tirggered averaging (eta) node and connect fluorescence and behavior nodes to eta node. And Run the workflow. 



<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig12.06_etaworkflow.png" alt="Whole" />
</p>

After finishing the process, you can quickly confirm your event-triggered average plot by clicking the OUTPUT button on the eta node. This figure is also available at VISUALIZE page.

<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig12.07_etamean.png" alt="Whole" />
</p>

The plots are for quick confirmation of the results. If you want to look into the results more in detail,   available variables are all saved in the OptiNiSt output in NWB format. They are saved in processing/optinist. The NWB file is easily retrieved at RECORD page with just one click. To inspect the data, [HDFView](https://www.hdfgroup.org/downloads/hdfview/) is convenient. 

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig12.08_resulthdf.png" alt="Whole" />
</p>






### Additional information on WORKFLOW

##### setting OPTINIST_DIR
The file assigning the OPTINIST_DIR is optinist/optinist/api/dir_path.py. Change line for OPTINIST_DIR, INPUT_DIR, and OUTPUT_DIR according to your demand. Changing dir_path.py may also be necessary when running the pipeline on your cluster computers. Also, you can quickly change OPTINIST_DIR by changing the environment variable by typing 'export OPTINIST_DIR="your_saving_dir"' before launching.

#### about the assumed data shape 
eta, cca, correlation, cross_correlation, granger, glm, lda, and svm assume the input neural data shape is frames x cells matrix. Because the output of CaImAn and Suite2P on the pipeline is cell x frames, the default setting for neural data for these analyses is set to transpose. 

Pca and tsne can be done in either direction depending on your purpose. The function assumes their input to be samples x features.  


##### snakemake settings
For details about snakemake parameters please refer to [here](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
use_conda: If this is on, snakemake uses conda environment.<br>
cores: Specifies the number of cores to use. If not specified, snakemake uses number of available cores in the machine. <br>
forceall: Flag to indicate the execution of the target regardless of already created output.<br>
forcetargets:  <br>
lock:  <br>

##### NWB settings
For detais about NWB please refer to [here](https://pynwb.readthedocs.io/en/latest/pynwb.file.html)
session_description: a description of the session where this data was generated <br>
identifier: a unique text identifier for the file  <br>
experiment_description: general description of the experiment <br>
device: device used to aquire the data (information such as firmware version, model etc.) <br>
optical_channel: channel used to acqure the data <br>
imaging_plane:  <br>
image_serises: <br>
ophys: <br>





## Inspecting the images and the plots on VISUALIZE
After executing the pipeline, you may want to check and compare the results.
VISUALIZE page is the place to work on this. You can replay the tiff time-series, see the cell ROI images, the plot of cell fluorescence or spike time-series, and other plots showing the results of analyses. See [here](../gui/visualize.md) for basic usage.


### Checking movies
You may want to check some frames of the multi-page tiff files. Visualize page offers the way to check. After creating a plot box by clicking on + mark, Select the image using the SELECT IMAGE button on the left top.
You can select the range of the frame by assigning 1st and last frame numbers. LOAD button starts loading the data.

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig21_loadmovie.png" alt="Whole" />
</p>

Click on the PLAY button within the plotting box to play the loaded movie.
The number indicated on the right of PAUSE button is the frame interval in milliseconds. 

<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig22_movie.png" alt="Whole" />
</p>


### Showing ROI and time courses
After running the ROI detection algorithms, the most often created plots are extracted cells' shape and fluorescence time series. To show the plot, prepare two plotting boxes.

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig23_twobox.png" alt="Whole" />
</p>

In one plotting box (ex, the one with ID:0), select a background image such as meanimg from the Select Item pulldowns.

<br>
<p align="left">
<img width="100px" src="./_static/tutorials/fig24_selectitem.png" alt="Whole" />
</p>

In the same plotting box, select cell_roi from the Select Roi pull-downs. Both Suite2P and CaImAn include the process to drop the extracted ROIs that do not meet the criteria. In OptiNiSt, the cell ID is given to all the ROIs. Cell_roi is the ROIs that passed the criteria. 

<br>
<p align="left">
<img width="100px" src="./_static/tutorials/fig25_selectroi.png" alt="Whole" />
</p>

The plotting box (ID:0) shows the background image and detected cells.
<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig26_roi.png" alt="Whole" />
</p>

In another plotting box (ex, the one with ID:1), select fluorescence from the Select Item pulldown.
And select 0(same ID with the plotting box of your ROI image) from the ref image pull down. By doing this,  the two plotting boxes are linked. 

<br>
<p align="left">
<img width="400px" src="./_static/tutorials/fig27_fluo.png" alt="Whole" />
</p>

Now you can explore the ROI and time course. The color of ROI and corresponding time course is matched. You will know the cell ID by letting your mouse over the cell in the image. Clicking on the cell automatically adds the fluorescence time course of the clicked cell. 
<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig28_roifluo.png" alt="Whole" />
</p>

If it is tiring to select the cell by clicking one by one, turn on the drag select button on the right in the plotting box of ROI. It enables selecting all the cells within the rectangular area you define.

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig29_dragselect.png" alt="Whole" />
</p>

### Savind plots
You can save created plots in svg, png, jpeg, or webp format. Please select the format, decide the saving name in the lower area on the left panel, and click the camera mark in the plotting box. Svg format saves the plot as a vector-based graphical format which may be convenient when you need high-resolution figures.

<br>
<p align="left">
<img width="200px" src="./_static/tutorials/fig30_saving.png" alt="Whole" />
</p>

## Managing pipelines on RECORD
RECORD section keeps your analysis pipeline easy to organize and easy to retrieve. For the basic usage of the RECORD page, see also [here](../gui/record.md) The RECORD page shows the summary of current status of OPTINIST_DIR/output. 

<br>
<p align="left">
<img width="600px" src="./_static/tutorials/fig40_recordall.png" alt="Whole" />
</p>

Clicking the Reproduce arrow retrieves the pipeline onto the workflow. This function is convenient when you re-start analysis after closing the browser. The reproduced pipeline needs to be RUN again (not ALL RUN) to make plots available.

<br>
<p align="left">
<img width="100px" src="./_static/tutorials/fig41_reproduce.png" alt="Whole" />
</p>

The Download buttons for the workflow column and the NWB column copy the snakemake config or NWB file to your download folder. The snakemake config file contains the workflow information and parameters for each node. The NWB file contains the data and its analysis results. This function is convenient when users share the same analysis pipeline or inspect the output results.

<br>
<p align="left">
<img width="150px" src="./_static/tutorials/fig42_workflownwb.png" alt="Whole" />
</p>
