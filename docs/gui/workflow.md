Workflow
=================

OptiNiSt can help creating your analysis pipelines easily on the GUI. In this workflow field, you select the data and the algorithms or analysis methods (nods). Connecting these nods defines the order of processing (pipelines). The analysis pipeline can be parallel or bifurcating.
<br>
<p align="center">
<img width="400px" src="../_static/workflow/whole.png" alt="workflow" />

<br/>

## Creating workflow
### Setting Input Images

1. Click `UPLOAD` button to upload your image files to `OPTINIST_DIR`.  
<!--   
(Large data files take long time to upload, so it copies to `OPTINIST_DIR` directly and load file more quickly. )-->
  If your files are in a remote place and your analysis computer is local, uploading (copying) input files to the local directory saves time for accessing them. 
  
2. Click `SELECT IMAGE` button to set the path to the data as the input. You can select one file or all the files in a folder. All the image tiff files in a folder are concatenated if you choose a folder.

** Currentlly, image files with {.tif, .TIF, .tiff, .TIFF} extensions are accepted. Other extensions will be added on request.


<br>
<p align="center">
<img width="200px" src="../_static/workflow/components/imageNode.png" alt="imageNode" />
<img width="200px" src="../_static/workflow/components/imageList.png" alt="imageNode" />
</p>
<br/>

### Adding Nodes
Select algorithms or analysis methods from the treeview on the left by clicking "+" button. 
<br>
<p align="center">
  <img width="300px" src="../_static/workflow/components/add_algorithm.png" alt="Add Algorithm" />
</p>

<br>

### Setting Parameters
Click `PARAM` button of the node to change parameters.
<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/setparam.png" alt="Set Parameter" />
</p>

<br/>

### Connecting Nods 
Connect colored connectors of the nodes by dragging your cursor from the output connector to the next input connector to create connecting edges. The color of the connector indicates the data type of the input and the output.
You can only connect the input and output connectors of the same color. 

**DataType List**
- <span style="color: red; ">ImageData</span>
- <span style="color: green; ">Suite2pData</span>
- <span style="color: orange; ">Fluorescence</span>
- <span style="color: yellow; ">Behavior</span>
- <span style="color: blue; ">Iscell</span>

<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/connect_edge.png" alt="Connect Algorithm" />
</p>

<br/>

### Removing Nods 
Clicking on the x mark on a node or on an edge removes it from the workflow field. 


### Other Data Formats As The Input
Some of the analyses, such as event-triggered averaging or GLM need timecourse of behavioral data as an input. A behavior node is best for it.
The format should be .csv. Rows and columns can be specified by `settings` appearing after selecting the csv data. Note that the number of data points has to be the same as the number of frames of image data. 
Fluo data node is for cell's fluorescence timecourse data given as .csv. 

Another data format prepared is hdf5. This format is compatible with the nwb data format.
CSV and hdf5 nodes have black output connectors. The edge connected to the black output connector can be connected to any input connector. Be careful; this means that it does not check the format correspondence between input and output.

<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/csv_connect.png" alt="CSV Connect" />
</p>



## Running pipelines
### RUN ALL
After making the workflow pipeline, click "RUN ALL" button to execute the computation.
RUN ALL command creates a new output folder in the OPTINIST_DIR and executes all the pipelines. Naming it helps you identify individual workflow in the RECORD field. Note the name you specify here is just a tag. The actual folder name is automatically assigned as 32 digit letter and number. See the RECORD field for more detail.

<br>
<p align="center">
<img width="350px" src="../_static/workflow/components/run_start.png" alt="Run Start" />
<img width="150px" src="../_static/workflow/components/run_name.png" alt="Run Name" />
</p>

Runninng...
During running, the indicators of the nodes show the progress of the computation. Green checkmarks indicate the termination of the process without errors.
Red checkmarks indicate the abortion of the process because of the error. You can check the error messages by moving the cursor over the checkmark.

<p align="center">
<img width="300px" src="../_static/workflow/components/running.png" alt="Running" />
</p>

### RUN
Once you `RUN ALL`, you may want to change the last pipeline's parameters or add a new node to the current pipeline.
Use the `RUN` button in such a case. It automatically detects the change and runs only the nodes which were changed. The additional results are overwritten into the same folder. If you leave the previous results for comparison, use `RUN ALL` instead.

<p align="center">
<img width="300px" src="../_static/workflow/components/run_checkpoint.png" alt="Run Checkpoint" />
</p>

<br />

### OUTPUT
Once the green checkmark appears, you can quickly visualize the output by clicking the `OUTPUT` button on the node.
`OUTPUT` button shows a new window for the resulting graph. See VISUALIZE for the details of the chart.

<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/run_finish.png" alt="Finish Run" />
<img width="100px" src="../_static/workflow/components/run_output.png" alt="Finish Run" />
</p>



<br/>

