
import argparse
import datetime
from pathlib import Path
import os
import generic_proc_idx
import PblIdxTest
from sys import exit

#p = '/home/mishler/Output_files'
p = Path.cwd()
#days = range(1,32)
#months = range(1,13)
#years = range(1990,2025)
parser = argparse.ArgumentParser()
parser.add_argument("product", help="The product to subset with neutral PBL info -- pblhtsonde1mcfarl,aossmps,microbasekaplus,etc. (String)")
parser.add_argument("site", help="The site to use -- SGP, ENA, HOU, NSA, ect. (String)")
parser.add_argument("facility", help="The facility for the given site -- C1,M1,E13,etc. (String)")
parser.add_argument("data_lvl", help="The data level to use -- a0, b1, c1, etc. (String)")
parser.add_argument("strt_yr", help="The year of the start date -- ex. 2019 (Int)", type=int)
parser.add_argument("strt_month", help="The month of the start date (Int 1-12)", type=int)#,choices=months)
parser.add_argument("strt_day", help="The day of start date (Int 1-31).", type=int)#,choices=days)
parser.add_argument("end_yr", help="The year of the end date -- ex. 2019 (Int)", type=int)
parser.add_argument("end_month", help="The month of the end date (Int 1-12)", type=int)#,choices=months)
parser.add_argument("end_day", help="The day of end date (Int 1-31)", type=int)#, choices=days)
parser.add_argument("sonde_facility", help="The facility the pblht file came from (same as facility if product is pblhtsonde1mcfarl) -- C1,M1,E13,etc. (String)")
parser.add_argument("home_dir", nargs='?', help="Chosen path for output files (/home/mishler). Will add 'Output_files' subdirectory to provided path", default=p)
args = parser.parse_args()
if Path(args.home_dir).exists()==False:
	print("The directory you provided does not exist.")
	exit()

output_drectory = Path(args.home_dir, 'Output_files')
#out_path = Path(args.out_dir)
site_loc = f'{args.site.lower()}{args.product}{args.facility.upper()}.{args.data_lvl.lower()}'
end_date = (datetime.datetime(year=args.end_yr, month=args.end_month, day=args.end_day)).strftime("%Y%m%d")
strt_date = (datetime.datetime(year=args.strt_yr, month=args.strt_month, day=args.strt_day)).strftime("%Y%m%d")

print(f'We will be generating neutral regime {args.product} sonding files from {site_loc} from {strt_date} to {end_date}.\n'
          f'Using output directory {output_drectory}')
newdir = output_drectory#args.home_dir.joinpath("Test_files")
newdir.mkdir(parents=True, exist_ok=True)

if args.product == 'pblhtsonde1mcfarl':
    PblIdxTest.pbl_idx(out_dir=newdir,site=args.site.lower(),facility=args.facility.upper(),lvl=args.data_lvl.lower(),s_yr=args.strt_yr,
    	s_mon=args.strt_month,s_day=args.strt_day,e_yr=args.end_yr, e_mon=args.end_month, e_day=args.end_day)
else:
    generic_proc_idx.product_idx(out_dir=newdir,product=args.product.lower(), site=args.site.lower(), facility=args.facility.upper(),lvl=args.data_lvl.lower(),
    	sond_fac=args.sonde_facility.upper(),s_yr=args.strt_yr,s_mon=args.strt_month, s_day=args.strt_day, e_yr=args.end_yr, e_mon=args.end_month,
     	e_day=args.end_day)	

