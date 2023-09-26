Record
=================
<br>
<p align="center">
<img width="400px" src="../_static/record/whole.png" alt="Whole"/>
</p>


In the RECORD field, you can check the workflow status in your OPTINIST_DIR and manage your analysis pipeline. The table lists all the pipelines in your OPTINIST_DIR. By clicking on the mark on the 2nd column, you can show the details of each pipeline.

OptiNiSt records and reproduces past workflow pipelines. It can download results in nwb format.


OptiNiSt can:
- Record past executed workflow
- Reproduce past workflow
- Download workflow config file
- Download snakemake config file
- Download results as NWB files

### Record Table

<p align="center">
<img width="400px" src="../_static/record/components/table.png" alt="Table"/>
</p>

- **Timestamp**: latest execution date and time.
- **ID**: unique ID. This is the directory name for the whole results of the pipeline.
- **Name**: user-defined workflow name.
- **Success**: success or failure (abortion with error) of execution.
- **Reproduce**: button to reproduce the workflow graph to the WORKFLOW field.
- **Workflow**: button to copy the workflow composition to your download folder on your computer.
- **SnakeFile**: button to copy snakemake config file to your download folder on your computer.
- **NWB**: button to copy the analysis results as NWB file to your download folder on your computer.
- **Delete**: button to delete the workflow from the OPTINIST_DIR.



### Details

<p align="center">
<img width="400px" src="../_static/record/components/details.png" alt="Details"/>
</p>

- **Function**: names of the nods (function).
- **nodeID**: unique ID. This is the directory name for the results of each node.
- **Success**: success or failure (abortion with error) of execution of the node.
- **NWB**: button to copy the analysis results as NWB file to your download folder on your computer.


### Reproduce Button

<br>
<p align="left">
<img width="100px" src="../_static/tutorials/fig41_reproduce.png" alt="Whole" />
</p>

Clicking the Reproduce arrow retrieves the pipeline onto the workflow. This function is convenient when you restart the analysis after closing the browser. The reproduced pipeline needs to be `RUN` again (not `RUN ALL`) to make plots available.

### Download Buttons

<br>
<p align="left">
<img width="300px" src="../_static/tutorials/fig42_workflownwb.png" alt="Whole" />
</p>

You can download 3 types of config files here.

- Workflow: the workflow and parameters information to reproduce in OptiNiSt GUI.
- Snakemake: the workflow and parameters information used by Snakemake.
- NWB: nwbfile which contains the data and its analysis results.

This function is convenient when users want to share the same analysis pipeline or inspect the output results.
