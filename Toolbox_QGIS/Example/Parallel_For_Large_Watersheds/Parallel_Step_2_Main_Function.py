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



###########################################################################################
#
# This is a python script to use BasinMaker to call another script which will define 
# routing structure for each subregion in Parallel way 
#
###########################################################################################


def Generaterouting_subregion(i):
	os.system("python Parallel_Step_2_Script_For_Each_Sub_Region.py  " + str(i))


import time
from joblib import Parallel, delayed
import os
import pandas as pd
import timeit


############ Variable needs to be modified to run this example ######     

### Define a output folder where to store subregion information 
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"

### The BasinMaker folder 
BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

### Define number of processors 
ncores = 3
########### Variable needs to be modified to run this example ######


### Define derived folder 
DataBase_Folder = os.path.join(BasinMaker_Folder,'Toolbox_QGIS','tests','testdata','Required_data_to_start_from_dem')
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder,'SubRegion_info')


### Start timer 
start = timeit.default_timer()

### Read sub-region information 
SubReg_info_main = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv'))
arg_instances = range(0,len(SubReg_info_main))
### run watershed delineation for each subregion 
Parallel(n_jobs=ncores, verbose=1, backend="threading")(delayed(Generaterouting_subregion)(i) for i in arg_instances)
End = timeit.default_timer()

print("use    ",start - End)
##  3372.4267685 sec