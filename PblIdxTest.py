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

def pbl_idx(out_dir, site, facility, lvl, s_yr, s_mon, s_day, e_yr, e_mon, e_day):
    import xarray as xt
    import numpy as np
    import datetime
    from pathlib import Path
    curr_yr = s_yr
    curr_mo = s_mon
    cur_day = s_day
    end_yr = e_yr
    end_mo = e_mon
    end_day = e_day
    site = site
    facility = facility
    lvl = lvl

    archive_path = r'/data/archive'
    out_dir = out_dir
    pbl_out = Path(out_dir, 'neut_pbl_hgt')
    print(f'pbl_out: {pbl_out}')
    newdir = pbl_out
    #newdir.mkdir(parents=True, exist_ok=True)
    site_loc = f'{site}pblhtsonde1mcfarl{facility}.{lvl}'
    pbl_in = f'{archive_path}/{site}/{site_loc}'
    site_tag = f'{site}{facility}_'
    end_day = datetime.datetime(year=end_yr, month=end_mo, day=end_day)
    pbl_vars = ['height_ss', 'layer', 'base_time', 'time_offset', 'atm_pres', 'air_temp', 'wspd', 'rh',
                'pbl_heigh_heffter', 'qc_pbl_height_heffter',
                'pbl_height_bulk_richardson_pt25', 'qc_pbl_height_bulk_richardson_pt25', 'pressure_gridded',
                'lapserate_theta_ss',
                'lapserate_theta_smoothed', 'atm_pres_ss', 'theta_ss', 'wspd_ss', 'richardson_number',
                'virtual_theta_ss',
                'bottom_inversion', 'top_inversion', 'lapserate_max', 'delta_theta_max', 'level_1_liu_liang',
                'level_2_liu_liang','pbl_height_bulk_richardson_pt5', 'qc_pbl_height_bulk_richardson_pt5']

    init_day = datetime.datetime(year=curr_yr, month=curr_mo, day=cur_day)
    curr_date = init_day
    comp_str = init_day.strftime("%Y%m%d")
    used_mon = curr_mo
    count = 0
    while curr_date < end_day:
        day_str = (init_day + datetime.timedelta(days=count)).strftime("%Y%m%d")
 
        #files = Path(pbl_in).glob('*.' + day_str + '.*')
        #num_files = 0
        #for file in files:
        #    num_files += 1
        #if num_files == 0:
        #    print("------No sonding files for", day_str)
        #    if comp_str[5] != day_str[5]:
        #        print("**Started processing month", used_mon + 1)
        #        used_mon += 1
        #    comp_str = day_str
        #    curr_date = curr_date + datetime.timedelta(days=1)
        #    count += 1
        #    continue
        ##print(f'{pbl_in}/*.{day_str}.*')
        #ds_object = xr.open_mfdataset(f'{pbl_in}/*.{day_str}.*', drop_variables=pbl_vars)

        #_, index = np.unique(ds_object['pbl_height_liu_liang'], return_index=True)

        #ds_object = ds_object.isel(time=index)

        #if len(ds_object['time']) != num_files:
        #    print("**ERROR!! Lost data on day:", day_str)

        #ds_object['pbl_regime_type_liu_liang'] = ds_object['pbl_regime_type_liu_liang'].astype(np.int8)
        #flag_values = ds_object['pbl_regime_type_liu_liang'].attrs['flag_values'].split(',')
        #flag_values = np.array(flag_values, dtype=np.int8)
        #ds_object['pbl_regime_type_liu_liang'].attrs['flag_values'] = flag_values

        #ds_object['pbl_height_liu_liang'] = ds_object['pbl_height_liu_liang'].astype(np.float32)
        #exit_flag = False
        #for var_name in ['lat', 'lon', 'alt']:
        #    try:
        #        ds_object[var_name] = ds_object[var_name].isel(time=0).drop('time')
        #    except:
        #        exit_flag = True
        #        print("------An error occured while making lat/lon/alt scalr on", day_str)
        #        if comp_str[5] != day_str[5]:
        #            print("**Started processing month", used_mon + 1)
        #            used_mon += 1
        #            comp_str = day_str
        #        break
        #if exit_flag:
        #    curr_date = curr_date + datetime.timedelta(days=1)
        #    count += 1
        #    continue
        #for var_name in ds_object.data_vars:
        #    units = ds_object[var_name].attrs['units']
        #    if units == 'unitless':
        #        ds_object[var_name].attrs['units'] = '1'

        #ds_object = ds_object.sortby('time') 
        ##ds_object.to_netcdf(f'{out_dir}/{pbl_out}/{site_tag}{day_str}_pbl_info.nc', format='NETCDF4')
        #ds_object.to_netcdf(f'{pbl_out}/{site_tag}{day_str}_pbl_info.nc', format='NETCDF4')
        print(f'{pbl_out}/{site_tag}{day_str}_pbl_info.nc')
        ##print(f'{out_dir}/{pbl_out}/{site_tag}{day_str}_pbl_info.nc')
        if comp_str[5] != day_str[5]:
            print("**Started processing month", used_mon + 1)
            used_mon += 1
        comp_str = day_str
        curr_date = curr_date + datetime.timedelta(days=1)
        count += 1

    print("Sonde files read and Neutral ABL Datasets Created...Program Ends!!")

 