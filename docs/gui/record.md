Record
=================

OptiNiSt records and reproduces past workflow pipelines. It can download results in nwb format.
<br>
<p align="center">
<img width="400px" src="../_static/record/whole.png" alt="Whole"/>
</p>

OptiNiSt can:
- Record past executed workflow
- Reproduce past workflow
- Download results as NWB files
- Download snakemake config file

### Record Table
- **Timestamp**: last executing time.
- **ID**: Unique ID. Saving directory name.
- **Name**: workflow name, registered at executing.
- **Success**: execution was successful or not.
- **Reproduce**: reproduce result to workflow.
- **SnakeFile**: Downalod snakemake config file.
- **NWB**: Downalod result as NWB file.
- **Delete**: Delete past workflow.

<p align="center">
<img width="400px" src="../_static/record/components/table.png" alt="Table"/>
</p>


### Details
- **Function**: Executing function name.
- **nodeID**: Function unique ID.
- **Success**: execution was successful or not.
- **NWB**: Downalod result as NWB files.


<p align="center">
<img width="400px" src="../_static/record/components/details.png" alt="Details"/>
</p>
