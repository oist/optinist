Getting Started
=================

To start OptiNiSt, you need to open a console and type `run_optinist` or,

`cd ~/optinist/`  
`python main.py`  

As a default, `OPTINIST_DIR` is assigned to /tmp/ in your computer.  
To specifically assign your OPTINIST_DIR in which OptiNiSt copies the data and saves the results,
type `export OPTINIST_DIR=~/your_dir` or change the line in optinist/api.dir_path.py directly

from
`_DEFAULT_DIR = '/tmp/optinist'`
to 
_DEFAULT_DIR = '/your_dir'`


The console shows the log (Fig1) once the startup is completed.

<br>
<p align="center">
<img width="400px" src="../_static/tutorials/Fig1.png" alt="Whole" />
</p>

Once you see this, open your web browser (Google Chrome is recommended) at localhost:8000.
You are ready to start if the Optinist window (Fig2) appears.


<br>
<p align="center">
<img width="400px" src="../_static/tutorials/Fig2.png" alt="Whole" />
</p>

Accepted data types
=================
Currently,  image data (multiple page-tiff), behavior data (.csv), and .hdf5 based file such as .nwb formats are supported. 

As the input to Suite2P or CaImAn,  .tiff files or .hdf5 (or .nwb) files can be used.
As the input to post-roi-detection analysis, .csv, and  .hdf5 format is supported.

Using WORKFLOW 
=================



How it works
=================
