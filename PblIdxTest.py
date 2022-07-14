import numpy as np
import os
import xarray as xr
from pathlib import Path
import act
import matplotlib.pyplot as plt
import matplotlib
import glob
import datetime
from datetime import datetime
from datetime import timedelta

def pbl_idx(
	out_dir,
	site,
	facility,
	lvl,
	s_yr,
	s_mon,
	s_day,
	e_yr,
	e_mon,
	e_day
):
    """
    Function to query the ARM research server archive 
    and consolodate pblhtsondemcfarl1 from a given day 
    into a single NETCDF file while dropping all variables except
    launch time, regime type, PBL height, and the associated 
    QC variables. These daily files can then be used to index 
    other data products by the the time/height of a neutral regime.
    If no pblhtsondemcfarl1 files are available for a given day for 
    the specified site/facility/data level, an error message will be 
    printed giving this information and the program will move on to the
    next day in the provided time range.

    Parameters
    ----------
    out_dir: string or Path object
    	path to the directory at which a subdirectory (neut_pbl_hgt) will
    	be created and generated NETCDF files will be stored.
    site: string
    	site of the facility at which the pblhtsonde data was recorded. 
    	e.g., hou, sgp, ena, nsa, etc.
    facility: string
    	facility at which the pblhtsonde data was recorded.
    	e.g., E13, C!, M1, etc.
    lvl: string
    	data level of the pblhtsondemcfarl1 files to use
    	e.g., a0,c1,b1,etc.
    s_yr: int
    	starting year for the range of daily files to be created
    	e.g., 2010, 2019, 2022
    s_mon: int
    	starting month for the range of daily files to be created
    	e.g., 01, 7, 07, 11
    s_day: day
    	starting day for the range of daily files to be created
    	e.g., 01, 7, 22, 31
    e_yr: int
    	ending year for the range of daily files to be created
    	e.g., 2010, 2019, 2022
    e_mon: int
    	ending month for the range of daily files to be created
    	e.g., 01, 7, 07, 11
    e_day: int
    	ending day of the range of daily files to be created (NOT INCLUSIVE)
    	e.g., 02, 8, 23, 31

    Returns
    -------
    None: Nothing is returned by the function. Files are generated and
    placed into an output directory.

    Example
    -------
		..code-block:: python

			from PblIdxTest import pbl_idx
    		pbl_idx(out_dir='/home/mishler/Output_files',site='sgp',facility='C1',lvl='c1',s_yr=2019,
    				s_mon=5, s_day=01, e_yr=2019, e_mon=6, e_day=01)

    """

    import xarray as xt
    import numpy as np
    import datetime
    from pathlib import Path
    #instantiating variables with parameter info
    curr_yr = s_yr
    curr_mo = s_mon
    cur_day = s_day
    end_yr = e_yr
    end_mo = e_mon
    end_day = e_day
    site = site
    facility = facility
    lvl = lvl
    #instantiating path string the server archive
    archive_path = r'/data/archive' #-------------------------------------------------------------------------------------------USE PATH OBJ!!!!
    out_dir = out_dir
    #Instantiating path to output directory for generated files
    pbl_out = Path(out_dir, 'neut_pbl_hgt') 
    newdir = pbl_out
    #Creating the directory if it does not yet exist
    newdir.mkdir(parents=True, exist_ok=True)
    #instantiating string to located data files on server
    site_loc = f'{site}pblhtsonde1mcfarl{facility}.{lvl}'
    #Setting path string to directory of pblhtsondemcfarl1 files
    pbl_in = f'{archive_path}/{site}/{site_loc}'#-------------------------------------------------------------------------------USE PATH OBJ!!!!
    #Setting a tag name for outputfiles
    site_tag = f'{site}{facility}_'
    #Setting varaible to hold date of final day in range of generated files (NONINCLUSIVE)
    end_day = datetime.datetime(year=end_yr, month=end_mo, day=end_day)
    #List of variables to drop of pblhtsondemcfarl1 files
    pbl_vars = ['height_ss', 'layer', 'base_time', 'time_offset', 'atm_pres', 'air_temp', 'wspd', 'rh',
                'pbl_heigh_heffter', 'qc_pbl_height_heffter',
                'pbl_height_bulk_richardson_pt25', 'qc_pbl_height_bulk_richardson_pt25', 'pressure_gridded',
                'lapserate_theta_ss',
                'lapserate_theta_smoothed', 'atm_pres_ss', 'theta_ss', 'wspd_ss', 'richardson_number',
                'virtual_theta_ss',
                'bottom_inversion', 'top_inversion', 'lapserate_max', 'delta_theta_max', 'level_1_liu_liang',
                'level_2_liu_liang','pbl_height_bulk_richardson_pt5', 'qc_pbl_height_bulk_richardson_pt5']
    #Setting variable to hold date of the first day in range of files
    init_day = datetime.datetime(year=curr_yr, month=curr_mo, day=cur_day)
    #Setting varaible to hold updated day as files are generated
    curr_date = init_day
    #Setting information to tell when a new month has begun processing
    comp_str = init_day.strftime("%Y%m%d")
    used_mon = curr_mo
    #looping over each day in the provided date range
    count = 0
    while curr_date < end_day:
    	#Setting variable with string of the current date
        day_str = (init_day + datetime.timedelta(days=count)).strftime("%Y%m%d")

 		#Counting the number pblhtsondemcfarl files available for the given day
 		#If none available, no new file will be generated for this day and program
 		# will move on to the next day in date range.
        files = Path(pbl_in).glob('*.' + day_str + '.*')
        num_files = 0
        for file in files:
            num_files += 1
        if num_files == 0:
            print("------No sonding files for", day_str)
            #Checking if a new month is being processed
            if comp_str[5] != day_str[5]:
                print("**Started processing month", used_mon + 1)
                used_mon += 1
            comp_str = day_str
            curr_date = curr_date + datetime.timedelta(days=1)
            count += 1
            continue
        #creating a dataset with the daily pblhtsondemcfarl1 files
        ds_object = xr.open_mfdataset(f'{pbl_in}/*.{day_str}.*', drop_variables=pbl_vars) #---------------------------CREATE PATH OBJ FOR THIS!!!

        #Finding the indices with UNIQUE pbl heights -- corresponds with a given launch time
        _, index = np.unique(ds_object['pbl_height_liu_liang'], return_index=True)
        
        #subsetting the dataset by the unique indices found above
        ds_object = ds_object.isel(time=index)
        
        #If, somehow, two launch time have the EXACT same pbl hgt, this will flag
        if len(ds_object['time']) != num_files:
            print("**ERROR!! Lost data on day:", day_str)

        #Ensuring values are of the correct type
        ds_object['pbl_regime_type_liu_liang'] = ds_object['pbl_regime_type_liu_liang'].astype(np.int8)
        flag_values = ds_object['pbl_regime_type_liu_liang'].attrs['flag_values'].split(',')
        flag_values = np.array(flag_values, dtype=np.int8)
        ds_object['pbl_regime_type_liu_liang'].attrs['flag_values'] = flag_values
        ds_object['pbl_height_liu_liang'] = ds_object['pbl_height_liu_liang'].astype(np.float32)

        #Prevents a file with dimension corruption from breaking code. 
        exit_flag = False
        for var_name in ['lat', 'lon', 'alt']:
            try:
                ds_object[var_name] = ds_object[var_name].isel(time=0).drop('time')
            except:
                exit_flag = True
                print("------An error occured while making lat/lon/alt scalr on", day_str)
                if comp_str[5] != day_str[5]:
                    print("**Started processing month", used_mon + 1)
                    used_mon += 1
                    comp_str = day_str
                break
        if exit_flag:
            curr_date = curr_date + datetime.timedelta(days=1)
            count += 1
            continue
        for var_name in ds_object.data_vars:
            units = ds_object[var_name].attrs['units']
            if units == 'unitless':
                ds_object[var_name].attrs['units'] = '1'

        #Making sure dataset times are in order
        ds_object = ds_object.sortby('time') 
        
        #Creating NETCDF file with daily dataset and storing file in 'neut_pbl_hgt' subdirectory
        # within the provided directory
        ds_object.to_netcdf(f'{pbl_out}/{site_tag}{day_str}_pbl_info.nc', format='NETCDF4')#--------------------CREATE PATH OBJ FOR THIS!!!!!
        #Checking if a new month has begun processing
        if comp_str[5] != day_str[5]:
            print("**Started processing month", used_mon + 1)
            used_mon += 1
        comp_str = day_str
        curr_date = curr_date + datetime.timedelta(days=1)
        count += 1
    #Alerting user that the program has finished
    print("Sonde files read and Neutral ABL Datasets Created...Program Ends!!")

 