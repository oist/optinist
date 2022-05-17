Getting Started
=================

# Opening the browser

To start OptiNiSt, you need to open a console and activate optinist environment `conda activate optinist` and type `run_optinist` or,

`cd ~/optinist/`  
`python main.py`  


The console shows the log (Fig1) once the startup is completed.

<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig1.png" alt="Whole" />
</p>

Once you see this, open your web browser (Google Chrome is recommended) at localhost:8000.
You are ready to start if the OptiNiSt page (Fig2) appears.


<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>


# Making Pipelines
After launching, the first page you see is the workflow page. The workflow page is a place to define the analysis pipeline. You determine the data you will analyze, select the type of the algorithm or analysis method you use, and set the parameters and the order of analysis.  


<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig3.png" alt="Whole" />
</p>

 
## Assigning input data path
As a default, it shows an image node (Fig3 arrow, Fig4). 
This node defines the path to the data to use.  


<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>


OptiNiSt uses OPTINIST_DIR for retrieving data and saving results. OptiNiSt searches input data from the OPTINIST_DIR/input directory and save the results to the OPTINIST_DIR/output directory. 

You may not want to change anything in your original data folder, or you may wish to make your data folder visible to OptiNist because the imaging data is sometimes huge and takes time to copy. You can take either strategy in assigning your data path.

1. Copy your original data file to OPTINIST_DIR and assign the data path to the copied data. See [uploading data to OPTNIST_DIR](#uploading-data-to-optinist_dir) This can be done from the GUI.  
2. Set OPTINIST_DIR as the data path by changing the dir_path.py file. See [setting optinist directory](#setting-optinist_dir)

Once the data is made accessible, the SELECT IMAGE button on the image node becomes possible to assign the file as the input to the pipeline. You can select a file or a folder. Choosing a folder makes all the tiff files in the shown sequence an input set of continuous frames.


<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>

### Selecting analysis methods
The left side of the window shows pull-downs for all available analysis methods (Fig5). Clicking on the + mark adds the analysis nodes to the workflow field. ROI detection tools (currently suite2P and CaImAn ) are in ‘Algorithm’ category, and all other pre-installed analyses are in ‘optinist’ category.

Let’s start with sample tiff data (mouse2p_2_long.tiff) and try suite2p ROI detection.
First, you need to determine the image you use. Select your image as explained [above](#assigning-input-data-path).
Once it is selected, it shows the name of files in the image node. 

<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>


Then, connect the analysis nodes in the order you like to process. 
The color of the connector indicates the input and the output of the node. Color corresponds to the format of data. For example, red is the image type format. You can only connect the same color. (Exception: black is an undefined data format. You can connect the black connector with any other connector, but be careful it does not check the consistency of input and output).


<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>

As for Suite2P, you might not use suite2P_registration (motion correction). In that case, you can connect the suite2p_file_convert to suite2p_roi directly. 

Also, you can perform motion correction of CaImAn (caiman_mc) and then perform suite2P_roi.

The nodes should be connected as long as the input and the output are of the same format type (same color).
Also, you can branch the flow. In the example, the two caiman_mc with different parameter settings are created, and the downstream from caiman_mc is also different. Each node's results are saved in a separate folder (See [RECORD](#record) part). 


### parameter button and output button on the node
Each node has PARAM button and OUTPUT button. 
PARAM shows the parameters. Edit this as you like. 
OUTPUT button is for the quick check of the result. The button becomes effective after the successful execution of the pipeline.
It shows the results as a graph. [VISUALIZE](#visualzing-images-and-plots) page offers a detailed comparison of the charts. 


Now you are ready to run the pipeline.
RUN button at the right top shows the pulldown choices. RUNALL runs all the process. RUNALL assigns a new folder for saving the results. On the other hand, RUN skips the already ran processes. It checks the difference from the previous pipeline, including the existence of the results and the parameter change. If they are detected, the downstream process after that process is re-executed. The results are overwritten into the existing folder.

<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>

When you click on the RUNALL, it shows the window to determine the folder name. This folder name is only for the user’s convenience. The actual folder name is long digit random letter+number. Further information about the structure of the saved results is here (LINK TO RECORD).

<br>
<p align="left">
<img width="200px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>

Next to the RUN button, there is the CANCEL button. You can abort the running pipeline with this button. It immediately cancels the current execution.


SNAKEMAKE and NWB SETTING buttons are for parameters for snakemake and nwb.
The pipeline construction of Optinist is based on snakemake (ref), which is the pipeline controlling tool for python scripts. The SNAKEMAKE parameter setting is [here](#snakemake-settings).


NWB SETTING defines the metadata for the NWB file as an output. The parameter you set here is only for your record and not used for the calculation inside Optinist. You can leave this as default.





#### uploading data to OPTNIST_DIR
Clicking on the UPLOAD button on the node opens the file explorer or finder so that you can select the data file. UPLOAD button copies the selected file to your OPTINIST_DIR/input. 


#### setting OPTINIST_DIR
The file assiging the OPTINIST_DIR is optinist/optinist/api/dir_path.py. Change line for OPTINIST_DIR, INPUT_DIR and OUTPUT_DIR according to your demand. Changing dir_path.py may also be necessary when running the pipeline on your cluster computers.

#### snakemake settings
use_conda:  ADD COMMENTS!
cores:  ADD COMMENTS!
forceall:  ADD COMMENTS!
forcetargets:  ADD COMMENTS!
lock:  ADD COMMENTS!







<!-- 
As a default, `OPTINIST_DIR` is assigned to /tmp/ in your computer.  
To specifically assign your OPTINIST_DIR in which OptiNiSt copies the data and saves the results,
type `export OPTINIST_DIR=~/your_dir` or change the line in optinist/api.dir_path.py directly

from
`_DEFAULT_DIR = '/tmp/optinist'`
to 
_DEFAULT_DIR = '/your_dir'`
-->



How it works
=================
