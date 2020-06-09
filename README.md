# Tracing_quantification_and_analysis

pipelines for detection and analysis of RV+ neurons and AAV+ pixels. 

## Prerequisites
* AAV tracings imaged with SP5 or SP8 laser scanning confocal microscope
* RV tracings imaged with epifluorescent microscope 
* images need to be single sections 

## Requirements
* FIJI (Fiji is just ImageJ, NIH)
* Python (written in Python 3.6)

## Workflow for AAV tracings

1.	```autonomous_pixel_detection.ijm```

2.	```counting_AAV.ijm``` and ```ROIset```

* ROIset needs to be re-adjusted for every image separately. 
  If section cutting is uneven, a mix between ROIsets needs to be used.
  If ROIset does not fit over the section even after adjustments, new ROIs need to be drawn
  
3.	```summary_AAV.py```
script needs to be in the folder of counting data 

4.	```analysis_AAV.py```


## Workflow for RV tracings

1.	```autonomous_neuron_detection.ijm```
Before, train ‘Trainable Weka Segmentation’ on 3 – 4 sections and save settings
Select path in line 34

2.	```counting_RV.ijm``` and ```ROIset```
ROIset needs to be re-adjusted for every image separately. 
If section cutting is uneven, a mix between ROIsets needs to be used.
If ROIset is off, new ROIs need to be drawn

3.	```summary_RV.py```
script needs to be in the folder of counting data 

4.	```analysis_RV.py```

## Outputs
* autonomous_xx_detection: stack of DAPI + GFP + segmented image
* counting_RV: two excel tables
  1) excel sheet containing table of ROI labels, counts, total area, average size and %area. 
  2) excel sheet contains a table with area and X and Y values of for single cells.
* counting_AAV: excel sheet containing table of ROI labels, counts, total area and pixel density [Pixel/um²]
* summary_xx: csv. file of all counting data combined
analysis_xx: csv. files with data for pivot tables for plots along rostro-caudal axis

#### 

_Author of macros and scripts Daniel Gehrlach of the 
[Gogolla Lab](https://www.neuro.mpg.de/gogolla) at 
[Max-Planck-Institute of Neurobiology](https://www.neuro.mpg.de/en)._

_Repository created by Caroline Weiand of the 
[Gogolla Lab](https://www.neuro.mpg.de/gogolla) at 
[Max-Planck-Institute of Neurobiology](https://www.neuro.mpg.de/en)._
