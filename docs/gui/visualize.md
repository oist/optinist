Visualize
=================

<p align="center">
<img width="400px" src="../_static/visualize/whole.png" alt="Whole" />
</p>

OptiNiSt visualizes the analysis results by plotly. 



## Adding a display box

<p align="center">
<img width="400px" src="../_static/visualize/components/set_display_box.png" alt="SetDisplayBox" />
</p>

Click the "+" button to create a display box.   
To add another display box below, click the "+" button shown below. 

<p align="center">
<img width="400px" src="../_static/visualize/components/add_column.png" alt="AddColumn" />
</p>

To add another display box to the right, click **︙** at the upper right side of the box and show the pulldown. Then, click `Insert into next columns`.

### Selecting an item to show

<p align="center">
<img width="400px" src="../_static/visualize/components/select_output_item.png" alt="SelectOutputItem" />
</p>

Pull down of the `Select Item` shows the available item to show. Select one of these items.

### Customizing visualization parameters

<p align="center">
<img width="400px" src="../_static/visualize/components/customize_param.png" alt="CustomizeParam" />
</p>

Select one of the display boxes by clicking inside of the box. The blue highlight of the box indicates the current selection of the display box. The parameters shown on the left are attached to the currently selected box. Details of the parameter are explained here.

## ROI and timecourse 
There is a way to link ROI plots and fluorescence time series. 
Create One box showing ROI and another box showing fluorescence. You can link two boxes by setting `ref image` in the fluorescence box to be the ID of the ROI box. (ID of the box is on the left upper side). By clicking on the ROI of a cell, you can visualize the corresponding fluorescence time course in the fluorescence box. By turning on the `drag select` in ROI box, you can select multiple cells in the image at once. 


## Ref Image

<p align="center">
<img width="500px" src="../_static/visualize/components/ref_image.png" alt="RefImage" />
</p>

TimeSeries plot can refer to a image plot to synthronize cell number index. Click or drag image plot so that cell number indexes are synthoronized in corresponding timeseries plot.

## Editing ROI

<p align="left">
<img width="400px" src="../_static/tutorials/edit-roi/box.png" alt="Whole" />
</p>

To the edit roi, prepare a plotting box. 

<p align="left">
<img width="200px" src="../_static/tutorials/fig24_selectitem.png" alt="Whole" />
</p>

In one plotting box (ex, the one with ID:0), select a background image such as meanImg from the Select Item pulldowns.
In the same plotting box, select cell_roi from the Select Roi pull-downs. 

<br>
<p align="left">
<img width="200px" src="../_static/tutorials/fig25_selectroi.png" alt="Whole" />
</p>
The plotting box (ID:0) shows the background image and detected cells.
<br>
<p align="left">
<img width="400px" src="../_static/tutorials/edit-roi/cell_roi_selected.png" alt="Whole" />
</p>

<p align="left">
<img width="400px" src="../_static/tutorials/edit-roi/add_roi_clicked.png" alt="Whole" />
</p>


You can click the <strong>Add ROI</strong> button then drag drop, resize the white cirle to change the new ROI position and size.
Press <strong>OK</strong> or <strong>Cancel</strong> button to Add or No

<p align="left">
<img width="400px" src="../_static/tutorials/edit-roi/roi_selected_merge_or_delete.png" alt="Whole" />
</p>

Or click on each cell ROI to delete ROI or merge ROIs (when you select 2 or more ROI cells)
Press <strong>Merge ROI</strong> or <strong>Delete ROI</strong> or <strong>Cancel</strong> button to Merge or Delete or No.


## saving plots
For your record, you can save the plots. Select the format and give the filename at the bottom of the parameters. Clicking on the camera mark in the box saves the plot figures to your download folder on your computer.

