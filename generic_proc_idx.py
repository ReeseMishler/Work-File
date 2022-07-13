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


def product_idx(out_dir,product, site, facility, lvl, sond_fac, s_yr, s_mon, s_day, e_yr, e_mon, e_day):
    import xarray as xr
    import act
    import numpy as np
    import datetime
    from pathlib import Path
    product = product
    site = site
    facility = facility
    sonde_facility = sond_fac
    lvl = lvl
    site_loc = f'/{site}{product}{facility}.{lvl}' 
    sonde_site_loc = f'{site}{sonde_facility}'  # ----THIS MIGHT BE WRONG because its the sonde file right? so facility might differ?
    micro_in = r"/data/archive/" + site + site_loc
    out_dir = out_dir
    pbl_directory = Path(out_dir, 'neut_pbl_hgt') #------WE WILL LIKELY NEED TO CHECK TO SEE IF THIS DIRECTORY EXISTS
    micro_out_dir = f'neut_{product}'
    out_dir = Path(out_dir, micro_out_dir)
    newdir = out_dir
    newdir.mkdir(parents=True, exist_ok=True)
    curr_yr = s_yr
    curr_mo = s_mon
    curr_day = s_day
    end_year = e_yr
    end_mo = e_mon
    end_day = e_day
    end_date = datetime.datetime(year=end_year, month=end_mo, day=end_day)
    date = datetime.datetime(year=curr_yr, month=curr_mo, day=curr_day)

    while date < end_date:
        day_str = date.strftime("%Y%m%d")
        out_tag = f'/{site}{product}{facility}_{day_str}'
        micro_in_day = f'{site_loc}.{day_str}*'
        micro_out_fn = f'{out_tag}_neuPbl.nc'

        ##################################### Variables to drop from AOSSMPS #############################################
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

        micro_in_path = f'{micro_in}{micro_in_day}' 
        files = Path(micro_in).glob('*.' + day_str + '.*')
        num_files = 0
        for file in files:
            num_files += 1
        if num_files == 0:
            print(f"------No {site}{product}{facility}.{lvl} files for {day_str}") 
            date = date + datetime.timedelta(days=1)
            continue
        # path from which to read in AOSSMPS files (1 month)
        
        if (product == 'microbasekaplus'):
        	print(f"Creating {product} dataset for {day_str}")

        #ADD TRY CASE? ERROR THROWN IF NOTHING IS GENERATED?
        dsA = act.io.armfiles.read_netcdf(micro_in_path, drop_variables=drp_vars)

        pbl_path = f'{pbl_directory}/{sonde_site_loc}_{day_str}*'
        #print(f'Path to Neut PBL Files: {pbl_path}')
        try:
        	pbl_obj = act.io.armfiles.read_netcdf(pbl_path)
        except:
        	print(f"------No {site}{sonde_facility} pbl info files for {day_str}")
        	date = date + datetime.timedelta(days=1)
        	continue
        
        ########### WE NEED TO APPLY THE DQR AND QC HERE!!!! ############################
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
        ## -------------SELECT NEUTRAL REGIME ONLY------------------------------
        pbl_obj = pbl_obj.sel(time=pbl_obj['pbl_regime_type_liu_liang'] == 0)

        ## subsetting Product by neutral boundary layer info
        aaa = pbl_obj['time'][0].values - np.timedelta64(90, 'm')  # 90 min before launch
        bbb = pbl_obj['time'][0].values + np.timedelta64(90, 'm')  # 90 min after launch
        hgt = pbl_obj['pbl_height_liu_liang'][0].values
        if product == 'microbasekaplus':
            tester = dsA.sel(time=slice(aaa, bbb), height=slice(0, hgt))
        else:
            tester = dsA.sel(time=slice(aaa, bbb))
        

        ## Setting up an end to the loop
        last = len(pbl_obj['time'])
        #print("Creating Neutral ABL MICROBASE Dataset...")
        # looping through the remaining neutral boundary layer times & concatenating into one dataset
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
                tester = xr.concat([tester, temp], dim='time')
            else:
                strt = pbl_obj['time'][idx].values - np.timedelta64(90, 'm')
                end = pbl_obj['time'][idx].values + np.timedelta64(90, 'm')
                hgt = pbl_obj['pbl_height_liu_liang'][idx].values
                if product == 'microbasekaplus':
                	temp = dsA.sel(time=slice(strt, end), height=slice(0, hgt))
                else:
                	temp = dsA.sel(time=slice(strt, end))
                tester = xr.concat([tester, temp], dim='time')
        ## Getting into correct format for DQR application
        tester.clean.cleanup()
        product_obj = tester
        ## Applying DQR's for total number concentration --------------------WHAT DQR TO APPLY AND VARIABLES TO DROP?!?!
        #print("Applying DQR to product dataset...")
        dqr_vars = None
        product_obj = act.qc.arm.add_dqr_to_qc(product_obj, variable=dqr_vars)

        # Normalizing quality control assessments
        product_obj.clean.normalize_assessment()

        # Output the new neutral abl microbase containing the applied dqr and qc variable
        product_path = f'{out_dir}/{micro_out_fn}'
        #print("Creating NETCDF of neutral PBL MICROBASE...")
        product_obj.to_netcdf(product_path, format='NETCDF4')
        date = date + datetime.timedelta(days=1)

    print("End Program!!")


