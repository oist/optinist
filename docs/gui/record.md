Record
=================

OptiNiSt records and reproduces past workflow pipelines. It can download results as nwb format.
<br>
<p align="center">
<img width="400px" src="../_static/record/whole.png" alt="Whole"/>
</p>

Can
- Record past executing workflow
- Reproduce past workflow
- Download result as NWB
- Download snakemake config file

### Record Table
- **Timestamp**: last executing time.
- **ID**: Unique ID. Saving directory name.
- **Name**: workflow name, registered at executing.
- **Success**: executing is success or not.
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
- **Success**: executing is success or not.
- **NWB**: Downalod result as NWB file.


<p align="center">
<img width="400px" src="../_static/record/components/details.png" alt="Details"/>
</p>
