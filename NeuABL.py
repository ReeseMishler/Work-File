import numpy as np
import os
import xarray as xr
from pathlib import Path
import act
import matplotlib.pyplot as plt
import matplotlib
import glob
from datetime import datetime
from datetime import timedelta
#%matplotlib inline

time_arr = []
b4_time = []
aftr_time = []
pbl_hgt = []
# Iterate over files in the directory
directory = r"C:\Users\Reese\Desktop\AOS Proj Data\pblht"
for filename in os.listdir(directory):
    #Setting path for the given sounding file
    f = os.path.join(directory, filename)
    # checking if it is a file
    if os.path.isfile(f):
        ds1 = xr.open_dataset(f)
        #Check if regime is neutral
        if ds1['pbl_regime_type_liu_liang'].values == 0.0:
            givent_time = ds1['time'][0].values # sonde launch time
            time_arr.append(givent_time)
            b = ds1['time'][0].values - np.timedelta64(90,'m') #subtracting 90 min frm launch time
            b4_time.append(b)
            a = ds1['time'][0].values + np.timedelta64(90,'m') #adding 90 min to launch time
            aftr_time.append(a)
            p = ds1['pbl_height_liu_liang'].values
            pbl_hgt.append(p)
#Creating dataset holding Neutral ABL info
ds8 = xr.Dataset({
    'launch_time': xr.DataArray(
            data   = time_arr,
            dims   = ['time'],
            coords = {'time': time_arr}
            ),
    'b4_time': xr.DataArray(
                data   = b4_time,
                dims   = ['time'],
                coords = {'time': time_arr}
                ),
    'aftr_time': xr.DataArray(
                data   = aftr_time,
                dims   = ['time'],
                coords = {'time': time_arr}
                ),
    'pbl_hgt': xr.DataArray(
                data   = pbl_hgt,
                dims   = ['time'],
                coords = {'time': time_arr}
                )
            },
        attrs = {'example_attr': 'this is a global attribute'}
    )
#Setting path for output netcdf file
pbl_out_path = r'C:\Users\Reese\Desktop\AOS Proj Data\Ouput_Files\May2019_Neut_Pbl.nc'
ds8.to_netcdf(pbl_out_path,format='NETCDF4')

#Varaibles/dimentions to drop from AOSsmps
drop_vars = ['diameter_mobility','base_time','time_offset','time_bounds','diameter_mobility_bounds', 'lower_size', 'dN_dlogDp', 'qc_dN_dlogDp', 'total_SA_conc', 'qc_total_SA_conc', 'dD_to_dSA',
            'total_V_conc', 'qc_total_V_conc', 'dD_to_dV', 'aerosol_flow', 'bypass_flow', 'sheath_flow', 'delay_time', 'geometric_mean',
            'geometric_std', 'mean', 'median', 'mode', 'sample_temperature', 'sample_relative_humidity', 'sample_pressure', 'mean_free_path',
            'gas_viscosity', 'reference_gas_temperature', 'reference_gas_pressure', 'reference_mean_free_path', 'reference_gas_viscosity',
            'sutherland_constant', 'diffusion_correction', 'multiple_charge_correction', 'nanoparticle_agglomerate_mobility_analysis',
            'status_flag', 'd50', 'low_voltage', 'high_voltage', 'hv_polarity', 'tube_diameter', 'tube_length', 'DMA_inner_radius',
            'DMA_outer_radius', 'DMA_characteristic_length']
#path from which to read in AOSSMPS files (1 month)
smps_in_path = r"C:\Users\Reese\Desktop\AOS Proj Data\smps\*.nc"
dsA = act.io.armfiles.read_netcdf(smps_in_path, drop_variables=drop_vars)

#subsetting AOSsmps by neutral boundary layer info
aaa = ds8['b4_time'][0].values #90 min before launch
bbb = ds8['aftr_time'][0].values #90 min after launch
tester = dsA.sel(time=slice(aaa, bbb))

#Setting up an end to the loop
last = len(ds8['b4_time'])
#looping through the remaining neutral boundary layer times & concatenating into one dataset
for idx in range(1,last):
    #Check to see if the new dataset start time occurs BEFORE the ending time of the previous dataset
    if ds8['b4_time'][idx].values <= ds8['aftr_time'][idx-1].values:
        # if the start time of an iteration overlaps with the endtime of the previous,
        # we set the new start time to 1 minute after the previous end time
        strt = ds8['aftr_time'][idx-1].values + np.timedelta64(1, 'm')
        end = ds8['aftr_time'][idx].values
        temp = dsA.sel(time=slice(strt, end))
        tester = xr.concat([tester, temp], dim='time')
    else:
        strt = ds8['b4_time'][idx].values
        end = ds8['aftr_time'][idx].values
        temp = dsA.sel(time=slice(strt, end))
        tester = xr.concat([tester, temp], dim='time')
#Getting into correct format for DQR application
tester.clean.cleanup()
#Applying DQR's for total number concentration
var2 = 'total_N_conc'
obj2 = act.qc.arm.add_dqr_to_qc(tester, variable = var2)
#Normalizing quality control assessments
obj2.clean.normalize_assessment()
smps_qc_inc_path = r'C:\Users\Reese\Desktop\AOS Proj Data\Ouput_Files\SMPSqcInc.nc'
obj2.to_netcdf(smps_qc_inc_path,format='NETCDF4')
#We also want to make a file with the qc value applied and tossed.
obj2.qcfilter.datafilter(variables = 'total_N_conc', rm_assessments=['Bad', 'Incorrect', 'Indeterminate', 'Suspect'])
smps_qc_applied = r'C:\Users\Reese\Desktop\AOS Proj Data\Ouput_Files\SMPSqcUsed.nc'
obj2.to_netcdf(smps_qc_applied,format='NETCDF4')
#-----------------------------------------------------------------------------------------------------------------------
#We next begin working with the Microbase product
micro_in_path=r'C:\Users\Reese\Desktop\AOS Proj Data\microbase\*.nc'
micro_drp_vars = ['base_time'] #Add or subtract as needed
micro = act.io.armfiles.read_netcdf(micro_in_path, drop_variables=micro_drp_vars)
#Setting start and end times for an initial neutral pbl dataset that can be added to
strt = ds8['b4_time'][0].values
end = ds8['aftr_time'][0].values
hgt = ds8['pbl_hgt'][0].values
#slicing the full microbase dataset by the start time, endtime, and pbl height of the first index of neutral ABL dataset
micro_pbl_neutrl = micro.sel(time=slice(strt, end), height=slice(0,hgt))

#Create a full dataset of neutral abl MICROBASE info
end = len(ds8['b4_time']) #THIS MUST BE SWITCHED IN BELOW LOOP AS END OF RANGE WHEN USING A MONTH
# Looping over the indexes of the neutral pbl file to concat with initial dataset
for idx in range(1,20): #Only looping through twenty because we are using three days of data!! switch out for 'end'
    #Check to see if the new dataset start time occurs BEFORE the ending time of the previous dataset
    if ds8['b4_time'][idx].values <= ds8['aftr_time'][idx-1].values:
        # Setting the new start time to be 1min after the previous end time
        strt = ds8['aftr_time'][idx-1].values + np.timedelta64(1, 'm')
        end = ds8['aftr_time'][idx].values
        hgt = ds8['pbl_hgt'][idx].values
        temp = micro.sel(time=slice(strt, end), height=slice(0,hgt))
        micro_pbl_neutrl = xr.concat([micro_pbl_neutrl, temp], dim='time')
    else:
        strt = ds8['b4_time'][idx].values
        end = ds8['aftr_time'][idx].values
        hgt = ds8['pbl_hgt'][idx].values
        temp = micro.sel(time=slice(strt, end), height=slice(0,hgt))
        micro_pbl_neutrl = xr.concat([micro_pbl_neutrl, temp], dim='time')
#Saving the fully subset (by neutral abl times/heights) to netCDF file
micro_out_path = r'C:\Users\Reese\Desktop\AOS Proj Data\Ouput_Files\Microbase_Neu_ABL.nc'
micro_pbl_neutrl.to_netcdf(micro_out_path,format='NETCDF4')