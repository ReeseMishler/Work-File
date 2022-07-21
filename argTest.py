#! /usr/bin/env python3
import argparse
import datetime
from pathlib import Path
import xarray as xr
import act
import numpy as np
import os
import generic_proc_idx 
from sys import exit

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
    other datastreams by the the time/height of a neutral regime.
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
    #archive_path = r'/data/archive' #-------------------------------------------------------------------------------------------USE PATH OBJ!!!!
    archive_path = Path('/data','archive')
    out_dir = out_dir
    #Instantiating path to output directory for generated files
    pbl_out = Path(out_dir, 'neut_pbl_hgt') 
    newdir = pbl_out
    #Creating the directory if it does not yet exist
    newdir.mkdir(parents=True, exist_ok=True)
    #instantiating string to located data files on server
    site_loc = f'{site}pblhtsonde1mcfarl{facility}.{lvl}'
    #Setting path string to directory of pblhtsondemcfarl1 files
    #pbl_in = f'{archive_path}/{site}/{site_loc}'#-------------------------------------------------------------------------------USE PATH OBJ!!!!
    pbl_in = Path(archive_path,site,site_loc)
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
        files = pbl_in.glob('*.' + day_str + '.*')
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
        read_in = Path(pbl_in, f'*.{day_str}.*') 
        
        ds_object = act.io.armfiles.read_netcdf(str(read_in), drop_variables=pbl_vars)
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
        path_out = Path(pbl_out,f'{site_tag}{day_str}_pbl_info.nc')
        #ds_object.to_netcdf(f'{pbl_out}/{site_tag}{day_str}_pbl_info.nc', format='NETCDF4')#--------------------CREATE PATH OBJ FOR THIS!!!!!
        ds_object.to_netcdf(path_out, format='NETCDF4')
        print(f"Created {site_tag}{day_str}_pbl_info.nc at path: {pbl_out}")
        #Checking if a new month has begun processing
        if comp_str[5] != day_str[5]:
            print("**Started processing month", used_mon + 1)
            used_mon += 1
        comp_str = day_str
        curr_date = curr_date + datetime.timedelta(days=1)
        count += 1
    #Alerting user that the program has finished
    print("Sonde files read and Neutral ABL Datasets Created...Program Ends!!")

def datastream_idx(
    out_dir,
    datastream,
    site,
    facility,
    lvl,
    sond_fac,
    s_yr,
    s_mon,
    s_day,
    e_yr,
    e_mon,
    e_day
): 
    """
    Function to query the ARM research server archive for a given datastream
    and subset said datastream by neutral pbl regime times/heights from a given day 
    into a NETCDF. These daily files can then be used to index 
    other data datastreams by the the time/height of a neutral regime.
    If no daily datastream files or pblhtsondemcfarl1 files (for which to subset the
    datastream) are available for a given day for the specified site/facility/data level,
    an error message will be printed giving this information and the program will move
    on to the next day in the provided time range.

    Parameters
    ----------
    out_dir: string or Path object
        path to the directory at which a subdirectory (neut_{datastream}) will
        be created and generated NETCDF files will be stored.
    datastream: string
        data datastream being subset with the neutral PBL time/height info
    site: string
        site of the facility at which the datastream data was recorded. 
        e.g., hou, sgp, ena, nsa, etc.
    facility: string
        facility at which the datastream data was recorded.
        e.g., E13, C!, M1, etc.
    lvl: string
        data level of the datastream files to use
        e.g., a0,c1,b1,etc.
    sond_fac: string
        the facility where the pblhtsonde files were recorded
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

            from generic_proc_idx import datastream_idx
            datastream_idx(out_dir='/home/mishler/Output_files',datastream='aossmps',site='sgp',facility='E13',lvl='b1',
                    sond_fac='C1',s_yr=2019,s_mon=5, s_day=01, e_yr=2019, e_mon=6, e_day=01)

    """

    #instantiating variables with parameter info
    keep_vars = ['total_N_conc','qc_total_N_conc','base_time', 'liquid_water_content', 'liquid_water_content_uncertainty_random',
    			'qc_liquid_water_content', 'ice_water_content', 'ice_water_content_uncertainty_random', 'qc_ice_water_content', 'liq_effective_radius',
				'liq_effective_radius_uncertainty_random', 'qc_liq_effective_radius', 'ice_effective_radius', 'ice_effective_radius_uncertainty_random', 
				'qc_ice_effective_radius', 'aqc_precip', 'lat', 'lon', 'alt']

    datastream = datastream
    site = site
    facility = facility
    sonde_facility = sond_fac
    lvl = lvl
    curr_yr = s_yr
    curr_mo = s_mon
    curr_day = s_day
    end_year = e_yr
    end_mo = e_mon
    end_day = e_day
    out_dir = out_dir
    #instantiating string to locate data files on server archive
    site_loc = f'/{site}{datastream}{facility}.{lvl}'#-------------------------------------------------------USE PATH OBJECT WHERRE THIS GOES TO GET RIF OF '/'
    site_loc2 = f'{site}{datastream}{facility}.{lvl}'
    sonde_site_loc = f'{site}{sonde_facility}'  
    #instantiating path string the server archive
    micro_in = Path('/data','archive',site,site_loc2)
    
    #Setting a path to the directory containing the neutral hgt/time files
    pbl_directory = Path(out_dir, 'neut_pbl_hgt') #----------------------------------------------------------------NEED TO CHECK AND SEE IF THIS EXISTS!!!!
    if Path(pbl_directory).exists()==False:
    	print("No directory exists constaining the daily pbl files for neutral regime indexing. Please create these first.")
    	return
    	#pbl_directory = Path('/home/mishler/Output_files', 'neut_pbl_hgt')
    
    #Tag for the output netcdf files
    micro_out_dir = f'neut_{datastream}'
    #Setting path to the directory netcdf files will be placed and created directory if it doesnt exists
    out_dir = Path(out_dir, micro_out_dir)
    newdir = out_dir
    newdir.mkdir(parents=True, exist_ok=True)
    #Setting varaible to hold date of final day in range of generated files (NONINCLUSIVE)
    end_date = datetime.datetime(year=end_year, month=end_mo, day=end_day)

    #Setting variable to hold date of the current day in range of files (set to start day initially)
    date = datetime.datetime(year=curr_yr, month=curr_mo, day=curr_day)
    #looping over each day in the provided date range
    while date < end_date:
        day_str = date.strftime("%Y%m%d")
        #tag for finding daily product file in archive
        out_tag = f'{site}{datastream}{facility}_{day_str}' #------------------------------------------WHERE DOES THIS GO?!? USE PATH OBJECT TO GET RID OF '/'
        micro_in_day = f'{site_loc}.{day_str}*' #f'{site_loc2}.{day_str}*'
        #Using day to name output NETCDF file
        micro_out_fn = f'{out_tag}_neuPbl.nc' 

        ##################################### Variables to drop from AOSSMPS & MICROBASE #############################################
        if datastream == 'microbasekaplus':
            drp_vars = ['time_offset', 'time_bounds', 'height_bounds', 'aqc_liquid_water_content',
                        'aqc_ice_water_content','aqc_liq_effective_radius', 'aqc_ice_effective_radius',
                        'aqc_retrieval', 'aqc_clear_cloud','mwr_scale_factor', 'aqc_stat2_lwp']
        elif datastream == 'aossmps':
            drp_vars = ['diameter_mobility', 'base_time', 'time_offset', 'time_bounds',
                        'diameter_mobility_bounds','lower_size', 'dN_dlogDp', 'qc_dN_dlogDp', 'total_SA_conc',
                        'qc_total_SA_conc','dD_to_dSA','total_V_conc', 'qc_total_V_conc', 'dD_to_dV', 'aerosol_flow',
                        'bypass_flow','sheath_flow','delay_time', 'geometric_mean', 'geometric_std', 'mean', 'median',
                        'mode','sample_temperature','sample_relative_humidity', 'sample_pressure', 'mean_free_path',
                        'gas_viscosity','reference_gas_temperature','reference_gas_pressure', 'reference_mean_free_path',
                        'reference_gas_viscosity','sutherland_constant','diffusion_correction', 'multiple_charge_correction',
                        'nanoparticle_agglomerate_mobility_analysis','status_flag', 'd50', 'low_voltage', 'high_voltage',
                        'hv_polarity', 'tube_diameter','tube_length','DMA_inner_radius', 'DMA_outer_radius', 'DMA_characteristic_length']
        else:
            drp_vars = []

        #Counting the number product files available for the given day
        #If none available, no new file will be generated for this day and program
        # will move on to the next day in date range.
        
        #micro_in_path = f'{micro_in}{micro_in_day}' #-------------------------------------------------------------USE PATH OBJECT wherever this goes!!!!!!!!1
        micro_in_path = Path(micro_in,micro_in_day)
        
        files = micro_in.glob('*.' + day_str + '.*')
        num_files = 0
        for file in files:
            num_files += 1
        if num_files == 0:
            print(f"------No {site}{datastream}{facility}.{lvl} files for {day_str} for path {micro_in_path}") 
            date = date + datetime.timedelta(days=1)
            continue

        #creating a dataset with the daily datastream files
        dsA = act.io.armfiles.read_netcdf(str(micro_in_path), drop_variables=drp_vars)
        
        #Setting a path string to the pblht files corresponding to this day
        #pbl_path = f'{pbl_directory}/{sonde_site_loc}_{day_str}*' #--------------------------------------------------------USE PATH OBJECT!!!!!
        
        pbl_path = Path(pbl_directory, f'{sonde_site_loc}_{day_str}*')

        #Attempting to create a dataset with pblht files. If none exists, the datastream
        # cannot be subset by neutral pbl information and will thus be skipped
        try:
        	pbl_obj = act.io.armfiles.read_netcdf(str(pbl_path))
        except:
        	print(f"------No {site}{sonde_facility} pbl info files for {day_str}")
        	date = date + datetime.timedelta(days=1)
        	continue
        
        #Ensuring the pblht file is in order WRT time, applying DQR for regime type and pbl height, and applying QC
        
        pbl_obj = pbl_obj.sortby('time')
        pbl_obj.clean.cleanup()
        var2 = 'pbl_regime_type_liu_liang'
        pbl_obj = act.qc.arm.add_dqr_to_qc(pbl_obj, variable=var2)
        pbl_obj.clean.normalize_assessment()
        pbl_obj.qcfilter.datafilter(variables='qc_pbl_regime_type_liu_liang', rm_assessments=['Bad', 'Incorrect'])
        var3 = 'pbl_height_liu_liang'
        pbl_obj = act.qc.arm.add_dqr_to_qc(pbl_obj, variable=var3)
        pbl_obj.clean.normalize_assessment()
        pbl_obj.qcfilter.datafilter(variables='qc_pbl_height_liu_liang', rm_assessments=['Bad', 'Incorrect'])
        # SELECT NEUTRAL REGIME ONLY from oblht dataset
        
        pbl_obj = pbl_obj.sel(time=pbl_obj['pbl_regime_type_liu_liang'] == 0)
        
        

        #######SETTING UP INTITAL PIECE OF NEUTRALLY SUBSET DATA datastream AND THEN ADDING PIECES TO IT#########
        # subsetting datastream by neutral boundary layer info
        # We are using 90 minutes before and after a given launch time w/ neutral regime
        try:
        	aaa = pbl_obj['time'][0].values - np.timedelta64(90, 'm')  # 90 min before launch
        	bbb = pbl_obj['time'][0].values + np.timedelta64(90, 'm')  # 90 min after launch
        	hgt = pbl_obj['pbl_height_liu_liang'][0].values 
        except:
        	print(f"No neutral regime times for {day_str}")
        	date = date + datetime.timedelta(days=1)
        	continue
        #microbasekaplus can be subset be the pblhgt as well as time
        if datastream == 'microbasekaplus':
            datastream_obj = dsA.sel(time=slice(aaa, bbb), height=slice(0, hgt))
        else:
            datastream_obj = dsA.sel(time=slice(aaa, bbb))
        

        ## Setting up an end to the loop
        last = len(pbl_obj['time'])

        # looping through the remaining neutral boundary layer times for day & concatenating into one dataset 
        for idx in range(1, last):
            # Check to see if the new dataset start time occurs BEFORE the ending time of the previous dataset
            if (pbl_obj['time'][idx].values - np.timedelta64(90, 'm')) <= (
                    pbl_obj['time'][idx - 1].values + np.timedelta64(90, 'm')):
                # if the start time of an iteration overlaps with the endtime of the previous,
                # we set the new start time to 1 minute after the previous end time
                strt = (pbl_obj['time'][idx - 1].values + np.timedelta64(90, 'm')) + np.timedelta64(1, 'm')
                end = pbl_obj['time'][idx].values + np.timedelta64(90, 'm')
                hgt = pbl_obj['pbl_height_liu_liang'][idx].values
                if datastream == 'microbasekaplus':
                	temp = dsA.sel(time=slice(strt, end), height=slice(0, hgt))
                else:
                	temp = dsA.sel(time=slice(strt, end))                
                datastream_obj = xr.concat([datastream_obj, temp], dim='time')
            else:
                strt = pbl_obj['time'][idx].values - np.timedelta64(90, 'm')
                end = pbl_obj['time'][idx].values + np.timedelta64(90, 'm')
                hgt = pbl_obj['pbl_height_liu_liang'][idx].values
                if datastream == 'microbasekaplus':
                	temp = dsA.sel(time=slice(strt, end), height=slice(0, hgt))
                else:
                	temp = dsA.sel(time=slice(strt, end))
                datastream_obj = xr.concat([datastream_obj, temp], dim='time')
        ## Getting into correct format for DQR application
        datastream_obj.clean.cleanup()

        ## Applying DQR's to each data variable available
        dqr_vars = None
        datastream_obj = act.qc.arm.add_dqr_to_qc(datastream_obj, variable=dqr_vars)

        # Normalizing quality control assessments
        datastream_obj.clean.normalize_assessment()

        # Output the new neutral abl datastream containing the applied dqr to neut_{datastream} subdirectory
        # of given output directory
        datastream_path = Path(out_dir,micro_out_fn)
        #product_path = f'{out_dir}/{micro_out_fn}' #-----------------------------------------------------------USE PATH OBJECT
        datastream_obj.to_netcdf(datastream_path, format='NETCDF4')
        print(f"Created {micro_out_fn} at path: {out_dir}")
        date = date + datetime.timedelta(days=1)

    print("End Program!!")  



# Setting the code that should be run when script is run as main function
############# Main Program ################
if __name__ == "__main__":

	p = Path.cwd()

	parser = argparse.ArgumentParser()
	parser.add_argument("data_stream", help="The data_stream to subset with neutral PBL info -- pblhtsonde1mcfarl,aossmps,microbasekaplus,etc. (String)")
	parser.add_argument("site", help="The site to use -- SGP, ENA, HOU, NSA, ect. (String)")
	parser.add_argument("facility", help="The facility for the given site -- C1,M1,E13,etc. (String)")
	parser.add_argument("data_lvl", help="The data level to use -- a0, b1, c1, etc. (String)")
	parser.add_argument("strt_yr", help="The year of the start date -- ex. 2019 (Int)", type=int)
	parser.add_argument("strt_month", help="The month of the start date (Int 1-12)", type=int)
	parser.add_argument("strt_day", help="The day of start date (Int 1-31).", type=int)
	parser.add_argument("end_yr", help="The year of the end date -- ex. 2019 (Int)", type=int)
	parser.add_argument("end_month", help="The month of the end date (Int 1-12)", type=int)
	parser.add_argument("end_day", help="The day of end date (Int 1-31)", type=int)
	parser.add_argument("sonde_facility", help="The facility the pblht file came from (same as facility if product is pblhtsonde1mcfarl) -- C1,M1,E13,etc. (String)")
	parser.add_argument("home_dir", nargs='?', help="Chosen path for output files (/home/mishler). Will add 'Output_files' subdirectory to provided path", default=p)
	args = parser.parse_args()
	home_path = Path(args.home_dir)
	home_path.mkdir(parents=True, exist_ok=True)

	output_drectory = Path(home_path, 'Output_files')


	site_loc = f'{args.site.lower()}{args.data_stream}{args.facility.upper()}.{args.data_lvl.lower()}'
	end_date = (datetime.datetime(year=args.end_yr, month=args.end_month, day=args.end_day)).strftime("%Y%m%d")
	strt_date = (datetime.datetime(year=args.strt_yr, month=args.strt_month, day=args.strt_day)).strftime("%Y%m%d")

	print(f'We will be generating neutral regime {args.data_stream} sonding files from {site_loc} from {strt_date} to {end_date}.\n'
	          f'Using output directory {output_drectory}')
	newdir = output_drectory
	newdir.mkdir(parents=True, exist_ok=True)
	
	if args.data_stream == 'pblhtsonde1mcfarl':
	    pbl_idx(out_dir=newdir,site=args.site.lower(),facility=args.facility.upper(),lvl=args.data_lvl.lower(),s_yr=args.strt_yr,
	    	s_mon=args.strt_month,s_day=args.strt_day,e_yr=args.end_yr, e_mon=args.end_month, e_day=args.end_day)
	else:
	    datastream_idx(out_dir=output_drectory,datastream=args.data_stream.lower(), site=args.site.lower(), facility=args.facility.upper(),lvl=args.data_lvl.lower(),
	    	sond_fac=args.sonde_facility.upper(),s_yr=args.strt_yr,s_mon=args.strt_month, s_day=args.strt_day, e_yr=args.end_yr, e_mon=args.end_month,
	     	e_day=args.end_day)	

