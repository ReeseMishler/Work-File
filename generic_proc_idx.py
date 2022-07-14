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


def product_idx(
    out_dir,
    product,
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
    Function to query the ARM research server archive -----------------------------------------FIX THIS DESCRIPTION
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
        path to the directory at which a subdirectory (neut_{product}) will
        be created and generated NETCDF files will be stored.
    product: string
        data product being subset with the neutral PBL time/height info
    site: string
        site of the facility at which the product data was recorded. 
        e.g., hou, sgp, ena, nsa, etc.
    facility: string
        facility at which the product data was recorded.
        e.g., E13, C!, M1, etc.
    lvl: string
        data level of the product files to use
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

            from generic_proc_idx import product_idx
            product_idx(out_dir='/home/mishler/Output_files',product='aossmps',site='sgp',facility='E13',lvl='b1',
                    sond_fac='C1',s_yr=2019,s_mon=5, s_day=01, e_yr=2019, e_mon=6, e_day=01)

    """


    import xarray as xr
    import act
    import numpy as np
    import datetime
    from pathlib import Path
    #instantiating variables with parameter info
    product = product
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
    site_loc = f'/{site}{product}{facility}.{lvl}'#-------------------------------------------------------USE PATH OBJECT WHERRE THIS GOES TO GET RIF OF '/'
    sonde_site_loc = f'{site}{sonde_facility}'  
    #instantiating path string the server archive
    micro_in = r"/data/archive/" + site + site_loc #--------------------------------------------------------------------USE PATH OBJECT
    #Setting a path to the directory containing the neutral hgt/time files
    pbl_directory = Path(out_dir, 'neut_pbl_hgt') #----------------------------------------------------------------NEED TO CHECK AND SEE IF THIS EXISTS!!!!
    #Tag for the output netcdf files
    micro_out_dir = f'neut_{product}'
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
        out_tag = f'{site}{product}{facility}_{day_str}' #------------------------------------------WHERE DOES THIS GO?!? USE PATH OBJECT TO GET RID OF '/'
        micro_in_day = f'{site_loc}.{day_str}*'
        #Using day to name output NETCDF file
        micro_out_fn = f'{out_tag}_neuPbl.nc' 

        ##################################### Variables to drop from AOSSMPS & MICROBASE #############################################
        if product == 'microbasekaplus':
            drp_vars = ['time_offset', 'time_bounds', 'height_bounds', 'aqc_liquid_water_content',
                        'aqc_ice_water_content','aqc_liq_effective_radius', 'aqc_ice_effective_radius',
                        'aqc_retrieval', 'aqc_clear_cloud','mwr_scale_factor', 'aqc_stat2_lwp']
        elif product == 'aossmps':
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
        micro_in_path = f'{micro_in}{micro_in_day}' #-------------------------------------------------------------USE PATH OBJECT wherever this goes!!!!!!!!1
        files = Path(micro_in).glob('*.' + day_str + '.*')
        num_files = 0
        for file in files:
            num_files += 1
        if num_files == 0:
            print(f"------No {site}{product}{facility}.{lvl} files for {day_str}") 
            date = date + datetime.timedelta(days=1)
            continue 
        #Microbasekaplus files take a while to create so this is to output progress
        if (product == 'microbasekaplus'):
        	print(f"Creating {product} dataset for {day_str}")

        #creating a dataset with the daily product files
        dsA = act.io.armfiles.read_netcdf(micro_in_path, drop_variables=drp_vars)
        
        #Setting a path string to the pblht files corresponding to this day
        pbl_path = f'{pbl_directory}/{sonde_site_loc}_{day_str}*' #--------------------------------------------------------USE PATH OBJECT!!!!!
        
        #Attempting to create a dataset with pblht files. If none exists, the product
        # cannot be subset by neutral pbl information and will thus be skipped
        try:
        	pbl_obj = act.io.armfiles.read_netcdf(pbl_path)
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
        

        #######SETTING UP INTITAL PIECE OF NEUTRALLY SUBSET DATA PRODUCT AND THEN ADDING PIECES TO IT#########
        # subsetting Product by neutral boundary layer info
        # We are using 90 minutes before and after a given launch time w/ neutral regime
        aaa = pbl_obj['time'][0].values - np.timedelta64(90, 'm')  # 90 min before launch
        bbb = pbl_obj['time'][0].values + np.timedelta64(90, 'm')  # 90 min after launch
        hgt = pbl_obj['pbl_height_liu_liang'][0].values
        #microbasekaplus can be subset be the pblhgt as well as time
        if product == 'microbasekaplus':
            product_obj = dsA.sel(time=slice(aaa, bbb), height=slice(0, hgt))
        else:
            product_obj = dsA.sel(time=slice(aaa, bbb))
        

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
                if product == 'microbasekaplus':
                	temp = dsA.sel(time=slice(strt, end), height=slice(0, hgt))
                else:
                	temp = dsA.sel(time=slice(strt, end))                
                product_obj = xr.concat([product_obj, temp], dim='time')
            else:
                strt = pbl_obj['time'][idx].values - np.timedelta64(90, 'm')
                end = pbl_obj['time'][idx].values + np.timedelta64(90, 'm')
                hgt = pbl_obj['pbl_height_liu_liang'][idx].values
                if product == 'microbasekaplus':
                	temp = dsA.sel(time=slice(strt, end), height=slice(0, hgt))
                else:
                	temp = dsA.sel(time=slice(strt, end))
                product_obj = xr.concat([product_obj, temp], dim='time')
        ## Getting into correct format for DQR application
        product_obj.clean.cleanup()

        ## Applying DQR's to each data variable available
        dqr_vars = None
        product_obj = act.qc.arm.add_dqr_to_qc(product_obj, variable=dqr_vars)

        # Normalizing quality control assessments
        product_obj.clean.normalize_assessment()

        # Output the new neutral abl product containing the applied dqr to neut_{product} subdirectory
        # of given output directory
        product_path = f'{out_dir}/{micro_out_fn}' #-----------------------------------------------------------USE PATH OBJECT
        product_obj.to_netcdf(product_path, format='NETCDF4')
        date = date + datetime.timedelta(days=1)

    print("End Program!!")


