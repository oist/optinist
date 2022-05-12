Workflow
=================

OptiNiSt can create your analysis pipelines easily on the GUI by connecting nodes and running them. Analysis flow combinations can be parallel or bifurcating. 
<br>
<p align="center">
<img width="400px" src="../_static/workflow/whole.png" alt="workflow" />

<br/>

## Creating workflow
### Setting Input Images

1. Click `UPLOAD` button to upload your image files to `OPTINIST_DIR`.  
  この部分不明
(Large data files take long time to upload, so it copies to `OPTINIST_DIR` directly and load file more quickly. )
  If your files are in the remote place and your analysis computer is local, uploading (copying) input files to local directory saves time for accessing to them. 
2. Click `SELECT IMAGE` button to set the path as the input _static. 
 

** Currentlly, files with {.tif, .TIF, .tiff, .TIFF} extensions are accepted. Other extensions will be added on the list if requested.
<br>
<p align="center">
<img width="200px" src="../_static/workflow/components/imageNode.png" alt="imageNode" />
<img width="200px" src="../_static/workflow/components/imageList.png" alt="imageNode" />
</p>
<br/>

### Adding Algorithms
Select an algorithm from the treeview on the left by clicking "+" button. Each of the BOX defines one step (node) of the analysis.
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

### Connecting BOX Edges
To create workflow pipeline, connect edges between node to node by drag&drop.  
You can only connect the same color edges (correspond to arguments and returns type, imageType, timeSeries type, ...).  

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

### CSV Input
Csv input could be any types, so it is shown in black color and could be connected with any arguments.
<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/csv_connect.png" alt="CSV Connect" />
</p>



## Run workflow
### Run all workflow
After you make workflow graph, click "Run" button and run workflow.
<br>
<p align="center">
<img width="350px" src="../_static/workflow/components/run_start.png" alt="Run Start" />
<img width="150px" src="../_static/workflow/components/run_name.png" alt="Run Name" />
</p>

Runninng...
<p align="center">
<img width="300px" src="../_static/workflow/components/running.png" alt="Running" />
</p>

### Run workflow from checkpoint
Add new algorithms and run from a checkpoint.
<p align="center">
<img width="300px" src="../_static/workflow/components/run_checkpoint.png" alt="Run Checkpoint" />
</p>

<br />

### Finish workflow
If the node earned a select tab in the algorithm node, it means that the node finished running.  
You can check result image, timseries or heatmap data in the visualize page.
<br>
<p align="center">
<img width="300px" src="../_static/workflow/components/run_finish.png" alt="Finish Run" />
<img width="100px" src="../_static/workflow/components/run_output.png" alt="Finish Run" />
</p>



<br/>

