Workflow
=================

# Workflow Page
OptiNiSt can make analysis pipelines through connecting nodes and run on GUI. There are many analysis flow combinations. It selects nodes and connects edges to create a pipelne flow.
<br>
<p align="center">
<img width="400px" src="../images/workflow/whole.png" alt="workflow" />

<br/>

## Create workflow
### Set Image Input
Set tiff image in your working directory and set as input.
1. Click `UPLOAD` button and upload file to `OPTINIST_DIR`.  
(An large data file takes long time tWWo upload, so it copies to `OPTINIST_DIR` directly and load file more quickly. )
2. Click `SELECT IMAGE` button and set as input images. 

** Currentlly, input images need to have {.tif, .TIF, .tiff, .TIFF} extension. (developing other extension)
<br>
<p align="center">
<img width="200px" src="../images/workflow/components/imageNode.png" alt="imageNode" />
<img width="200px" src="../images/workflow/components/imageList.png" alt="imageNode" />
</p>
<br/>

### Add Algorithm
Next, add an algorithm from left treeview.
1. Select an algorithm from the treeview
2. click "+" button to flowchart.
<br>
<p align="center">
  <img width="300px" src="../images/workflow/components/add_algorithm.png" alt="Add Algorithm" />
</p>

<br>

### Set Parameter
Click `PARAM` button and change parameters.
<br>
<p align="center">
<img width="300px" src="../images/workflow/components/setparam.png" alt="Set Parameter" />
</p>

<br/>

### Connect Algorithm Edge
To make workflow, connect edges between node to node by drag&drop.  
You can only connect the same color edges (correspond to arguments and returns type, imageType, timeSeries type, ...).  

**DataType List**
- <span style="color: red; ">ImageData</span>
- <span style="color: green; ">Suite2pData</span>
- <span style="color: orange; ">Fluorescence</span>
- <span style="color: yellow; ">Behavior</span>
- <span style="color: blue; ">Iscell</span>

<br>
<p align="center">
<img width="300px" src="../images/workflow/components/connect_edge.png" alt="Connect Algorithm" />
</p>

<br/>

### CSV Input
Csv input could be any types, so it is shown in black color and could be connected with any arguments.
<br>
<p align="center">
<img width="300px" src="../images/workflow/components/csv_connect.png" alt="CSV Connect" />
</p>



## Run workflow
### Run all workflow
After you make workflow graph, click "Run" button and run workflow.
<br>
<p align="center">
<img width="350px" src="../images/workflow/components/run_start.png" alt="Run Start" />
<img width="150px" src="../images/workflow/components/run_name.png" alt="Run Name" />
</p>

Runninng...
<p align="center">
<img width="300px" src="../images/workflow/components/running.png" alt="Running" />
</p>

### Run workflow from checkpoint
Add new algorithms and run from a checkpoint.
<p align="center">
<img width="300px" src="../images/workflow/components/run_checkpoint.png" alt="Run Checkpoint" />
</p>

<br />

### Finish workflow
If the node earned a select tab in the algorithm node, it means that the node finished running.  
You can check result image, timseries or heatmap data in the visualize page.
<br>
<p align="center">
<img width="300px" src="../images/workflow/components/run_finish.png" alt="Finish Run" />
<img width="100px" src="../images/workflow/components/run_output.png" alt="Finish Run" />
</p>



<br/>

