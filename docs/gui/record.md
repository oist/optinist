Record
=================
<br>
<p align="center">
<img width="600px" src="../_static/record/whole.png" alt="Whole"/>
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
<img width="600px" src="../_static/record/components/table.png" alt="Table"/>
</p>

| Header | Description |
| --- | --- |
| Timestamp | Latest execution timestamp. <br>It shows start time, end time and elapsed time of the workflow. |
| ID | Workflow's unique id. <br>This is same as the directory name for the whole results of the workflow. |
| Name | User-defined workflow name. You can edit the name by clicking the name. |
| Success | Workflow's status. Success, error or runnnig. |
| Reproduce | Button to reproduce the workflow to the WORKFLOW field. <br>You can visualize the results for the workflow by clicking the button. |
| Workflow | Button to download the workflow config yaml file. <br>This file can be used on import workflow button on WORKFLOW tab. <br>See details in [](ImportWorkflowYaml). |
| SnakeFile | Button to download the snakemake config file. |
| NWB | Button to download the analysis results as NWB file. |
| Delete | Button to delete the workflow record from the OPTINIST_DIR. |

### Details
You can check the results of each node by clicking arrow on the left of the table.

<p align="center">
<img width="600px" src="../_static/record/components/details.png" alt="Details"/>
</p>

| Header | Description |
| --- | --- |
| Function | Name of the node. |
| nodeID | Unique id of the node. <br>This is same as the directory name for the results of the node. |
| Success | Node's status. Success, error or runnnig. |
| NWB | Button to download the analysis results for the algorithm as NWB file. |

If status is error, you can see the error message by clicking status icon.

<p align="center">
<img width="600px" src="../_static/record/components/error_message.png" alt="ErrorMsg"/>
</p>

### Reproduce Button
Clicking the Reproduce arrow retrieves the pipeline onto the workflow.

This function is convenient when you want to check the results of the past workflow or reuse the same analysis pipeline.

<br>
<p align="left">
<img width="100px" src="../_static/tutorials/fig41_reproduce.png" alt="Whole" />
</p>

```{eval-rst}
.. note::
   You needed to RUN the workflow to visualize the results until before version 1.1.0.
   From version 1.1.0, you can visualize the results without running the workflow again.
```

### Download Buttons

<br>
<p align="left">
<img width="300px" src="../_static/tutorials/fig42_workflownwb.png" alt="Whole" />
</p>

You can download 3 types of config files here.

| FileType | Description |
| --- | --- |
| Workflow | It includes the workflow and parameters information to reproduce in OptiNiSt GUI. <br>This file can be used on import workflow button on WORKFLOW tab. <br>See details in [](ImportWorkflowYaml). |
| SnakeFile | It includes the workflow and parameters information used by Snakemake. |
| NWB | NWBfile which contains the data and its analysis results. |
