(data_nodes)=
Data Nodes
=================

OptiNiSt accepts a variety of data types. Data are added in "nodes", depending on their type. Here we summarise the accepted data of each node.

#### Image (TIFF)
- Format: 3D or 4D array shape (time, y, x) or (time, z, y, x).
- File Extension: ".tif", ".tiff", ".TIF", ".TIFF"
#### CSV
- Format: 2D data structure
  - Can have an optional header.
  - Can be transposed.
 <!-- Check ',' ';', Excel 'carriage return/line feed', Mac 'CR', UTF-8, UTF-16, and UTF-32-->
- File Extension: ".csv"
#### Fluo / Behavior
- Format: 2D data structure (similar to CSV)
  - Rows typically represent time points
  - Columns typically represent individual cells or ROIs
  - Can have an optional header
  - Can be transposed if specified in parameters
- File Extension: ".csv"
#### HDF5 / [NWB](https://nwb-schema.readthedocs.io/en/latest/format_description.html)
- Format: Contains Groups, Datasets, Attributes (metadata).
  - A Group is similar to a folder and may contain an arbitrary number of other groups and datasets,
  - A Dataset describes an n-dimensional array and provides the primary means for storing data,
  - An Attribute* is a small dataset that is attached to a specific group or dataset and is typically used to store metadata specific to the object they are associated with.
 - ".hdf5", ".nwb", ".HDF5", ".NWB"
#### Matlab
 - Format:
 - File Extension: ".mat"
#### Microscope
 - File Extensions: ".nd2", ".oir", ".isxd", ".thor.zip"

##### Inscopix(.isxd)
- We currently use `py_isx <https://github.com/inscopix/py_isx>`_ package with version 1.0.* for Inscopix data. This library only supports CellSet and Movie file types.

##### NIKON(.nd2)

##### Olympus(.oir)
- We currently implement reading NIKON and Olympus files using a C library. This is only available on Linux and Windows (not on MacOS yet, sorry).
- If you use Linux, use OptiNiSt by :ref:`developer mode <native_platforms_developer>`.
     - open ``studio/app/optinist/wrappers/optinist/conda/microscope.yaml`` and uncomment the ``- gcc=12`` line.

<!-- Data classes:
    ImageData:

Defined in image.py

Can accept data as:
A string filepath
A list of string filepaths
A numpy array

For numpy array input:
Saves data as a .tif file
Expects 3D or 4D array shape (time, y, x) or (time, z, y, x)


    TimeSeriesData:
Defined in timeseries.py
Expects 1D or 2D numpy array
If 1D, converts to 2D with shape (1, length)
Can also accept CSV filepath as input

    CsvData:
Defined in csv.py
Can accept:
String filepath to CSV
1D or 2D numpy array
If 1D array, converts to 2D with shape (1, length)

    HeatMapData:
Defined in heatmap.py
Expects 2D numpy array

    HistogramData:
Defined in histogram.py
Expects 1D numpy array
Reshapes to 2D array with shape (1, length)

    ScatterData:
Defined in scatter.py
Expects 1D or 2D numpy array
Transposes input data

    BarData:
Defined in bar.py
Expects 1D or 2D numpy array
If 1D, converts to 2D with shape (1, length)

    LineData:
Defined in line.py
Expects 2D numpy array

    PieData:
Defined in pie.py
Expects 1D numpy array
Reshapes to 2D array with shape (1, length)

    PolarData:
Defined in polar.py
Expects 2D numpy array

    HTMLData:
Defined in html.py
Expects string HTML content -->
