# Work-File
## AOS proj ipynb
Containing ipynb for AOS proj (contains code for different monthly plots of AOS data)
## TestBook ipynb
Also contains TestBook ipynb node book containing code which sorts through sonding data in a given directory, filters out files with neutral boundary layer, strips time of launch, height of boundary layer, and 90min post/pre launch before writing this info to new netcdf file.

Goes on to use the neutral boundary layer information to subset datastreams by the neutral ABL times/heights before saving those subsets to new netcdf files.
## NeuABL.py
This file does what the TestBook does, but it is put into a python script as opposed to the notebook.
## argTest.py
Acts as a driver script for calling functions that will either create daily pblhgt sonding files or subset a given product by hgt/times of neutral boundary layer. Used modules are incorporated into the script. Output files will have site and facility layer depth. Multiple daily files will now be generated if multiple daily files already exist in datastream (sans pblhtsonde1cfarl).
## PblIdxTest.py
Called by argTest.py if user wishes to create daily pblHgt files with launch times and pbl hgt.
## generic_proc_idx.py
called by argTest.py to subset a given aos product by the pbl hgt/time over a given period
