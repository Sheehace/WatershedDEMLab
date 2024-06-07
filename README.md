# WatershedDEMLab
WatershedDEMLab is a simple, Landlab-based routine for constructing landscape evolution models of real-world watersheds. Users supply a digital elevation model of a single watershed in which all pixels are assumed to flow to a common pour point. The routine then allows users to specify a range of parameters including fluvial erosion / sediment transport, hillslope / soil processes, and spatially / temporally variable uplift, climate, rock type, and land use conditions.


The WatershedDEMLab-master folder contains three subfolders:


-**Demdata:** A folder containing example DEMs. It is also a convenient location for users to store their own DEMs.

**Input:** A folder containing the Input.csv file. This file feeds parameter values into the WatershedDEMlab.py script.

**WatershedDEMLab:** A folder containing the WatershedDEMLab.py script. Running this script creates a landscape evolution model using the parameters prescribed in the Input.csv file.


The WatershedDEMlab routine is contained entirely within the WatershedDEMLab-master folder.  Users should follow this procedure:


1)	Download the WatershedDEMLab-master folder onto your computer. 

2)	Open the Input.csv file, specify your parameter values, and save the file. One of the parameters is the path to a DEM.

3)	Run the Watershed.py script. The script will run a landscape evolution model and save the resulting data in a directory specified in the Input.csv file. 


WatershedDEMLab is still in development and has many limitations. Future (before 2025) iterations will introduce additional and more flexible model parameters, detailed descriptions of parameter values, and example / tutorial Jupyter notebooks. It will also be repackaged into a single python function (or set of functions) that remove the mandatory file nesting scheme outlined above. 
