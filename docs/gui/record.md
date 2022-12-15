Record
=================
In the RECORD field, you can check the workflow status in your studio_DIR and manage your analysis pipeline. The table lists all the pipelines in your studio_DIR. By clicking on the mark on the 2nd column, you can show the details of each pipeline.

studio records and reproduces past workflow pipelines. It can download results in nwb format.
<br>
<p align="center">
<img width="400px" src="../_static/record/whole.png" alt="Whole"/>
</p>


studio can:
- Record past executed workflow
- Reproduce past workflow
- Download results as NWB files
- Download snakemake config file

### Record Table
- **Timestamp**: latest execution date and time.
- **ID**: unique ID. This is the directory name for the whole results of the pipeline.
- **Name**: user-defined workflow name.
- **Success**: success or failure (abortion with error) of execution.
- **Reproduce**: button to reproduce the workflow graph to the WORKFLOW field. 
- **SnakeFile**: button to copy snakemake config file to your download folder on your computer.
- **NWB**: button to copy the analysis results as NWB file to your download folder on your computer.
- **Delete**: button to delete the workflow from the studio_DIR.

<p align="center">
<img width="400px" src="../_static/record/components/table.png" alt="Table"/>
</p>


### Details
- **Function**: names of the nods (function).
- **nodeID**: unique ID. This is the directory name for the results of each node.
- **Success**: success or failure (abortion with error) of execution of the node.
- **NWB**: button to copy the analysis results as NWB file to your download folder on your computer.


<p align="center">
<img width="400px" src="../_static/record/components/details.png" alt="Details"/>
</p>
