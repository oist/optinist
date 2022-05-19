Visualize
=================


OptiNiSt visualizes the analysis results by plotly. 
<br>
<p align="center">
<img width="400px" src="../_static/visualize/whole.png" alt="Whole" />
</p>



### Adding a display box
Click the "+" button to create a display box.   
To add another display box below, click the "+" button shown below. 

<p align="center">
<img width="400px" src="../_static/visualize/components/set_display_box.png" alt="SetDisplayBox" />
</p>

To add another display box to the right, click **︙** at the upper right side of the box and show the pulldown. Then, click `Insert into next columns`.

<p align="center">
<img width="400px" src="../_static/visualize/components/add_column.png" alt="AddColumn" />
</p>


### Selecting an item to show
Pull down of the `Select Item` shows the available item to show. Select one of these items.

<p align="center">
<img width="400px" src="../_static/visualize/components/select_output_item.png" alt="SelectOutputItem" />
</p>


### Customizing visualization parameters
Select one of the display boxes by clicking inside of the box. The blue highlight of the box indicates the current selection of the display box. The parameters shown on the left are attached to the currently selected box. Details of the parameter are explained here.


### Add display box into columns right
Click **︙** and `Insert into next columns` to add another display box.
<p align="center">
<img width="400px" src="../_static/visualize/components/customize_param.png" alt="CustomizeParam" />
</p>


## ROI and timecourse 
There is a way to link ROI plots and fluorescence time series. 
Create One box showing ROI and another box showing fluorescence. You can link two boxes by setting `ref image` in the fluorescence box to be the ID of the ROI box. (ID of the box is on the left upper side). By clicking on the ROI of a cell, you can visualize the corresponding fluorescence time course in the fluorescence box. By turning on the `drag select` in ROI box, you can select multiple cells in the image at once. 


## Ref Image
TimeSeries plot can refer to a image plot to synthronize cell number index. Click or drag image plot so that cell number indexes are synthoronized in corresponding timeseries plot.
<p align="center">
<img width="500px" src="../_static/visualize/components/ref_image.png" alt="RefImage" />
</p>

## saving plots
For your record, you can save the plots. Select the format and give the filename at the bottom of the parameters. Clicking on the camera mark in the box saves the plot figures to your download folder on your computer.

