# Work-File
## AOS proj ipynb
Containing ipynb for AOS proj (contains code for different monthly plots of AOS data)
## TestBook ipynb
Also contains TestBook ipynb node book containing code which sorts through sonding data in a given directory, filters out files with neutral boundary layer, strips time of launch, height of boundary layer, and 90min post/pre launch before writing this info to new netcdf file.

Goes on to use the neutral boundary layer information to subset datasets of AOSSMPS and MICROBASEKAPLUS by the neutral ABL times/heights before saving those subsets to new netcdf files.
## NeuABL.py
This file does what the TestBook does, but it is put into a python script as opposed to the notebook.
