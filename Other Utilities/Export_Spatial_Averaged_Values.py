# Copyright 2018-2018 Ming Han - ming.han(at)uwaterloo.ca
#
# License
# This file is part of Ming Han's personal code library.
#
# Ming Han's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ming Han's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Ming Han's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#



import datetime as dt  # Python standard library datetime  module
import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
import os
import arcpy
from arcpy import env
from arcpy.sa import *
from simpledbf import Dbf5
from dbfpy import dbf
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

############### part1 create empty nc sample file with all possible variable
WorkingFolder = sys.argv[1]
ncfilename = sys.argv[2]
############### Change the input in this section
#WorkingFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/RDRS'
#ncfilename = 'RDRS_CaPA24hr_forcings_final.nc'
#catply = 'finalcat.shp'
###############################################
varnm   = ['TRAF_60268832','ALAT_0','O1_0','PR_0','AHFL_0','WT_59868832','R2CH'] ##variable name that will be included in the NetCDF file
Outdata =   
ncfile = os.path.join(WorkingFolder,ncfilename)
dsin2 =  Dataset(ncfile,'r')


ncin =  Dataset(ncfile,'r')
datevar = []
t_unit  = ncin.variables['time'].units
nctime = ncin.variables['time'][:]
t_cal =ncin.variables['time'].calendar
datevar = [num2date( x ,units = t_unit,calendar = t_cal) for x in nctime]

outdatahr_var = pd.DataFrame(np.full((len(datevar),len(varnm)),np.nan),index = datevar,columns=varnm)


for ivar in range(0,len(varnm)):
    varname         = varnm[ivar]
    print(varname)
    variable_values = ncin.variables[varname][:,:,:]
    Sumed           = np.sum(variable_values,axis = 0)
    print(Sumed[0])

ncin.close()


