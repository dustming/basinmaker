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
# This is a python script to use BasinMaker to divide a large watershed domain into
# several sub-regions
#
###########################################################################################


import os
import shutil
import sys
import tempfile
import timeit

import pandas as pd
from simpledbf import Dbf5

from basinmaker import basinmaker

############ Variable needs to be modified to run this example ######

### Define a output folder where to store subregion information
num  = str(np.random.randint(1, 10000 + 1))
path_working_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_" +num,"work")
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"
### The BasinMaker folder
datafolder = "../../tests/testdata/Required_data_to_start_from_dem"
########### Variable needs to be modified to run this example ######
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder, "SubRegion_info")

### Start timer
start = timeit.default_timer()

### Define input paths
na_hydem = os.path.join(datafolder, "DEM_big_merit.tif")  #'HydroSHED15S.tif')#
in_lake = os.path.join(datafolder, "hylake.shp")

### Initialize the BasinMake
basinmaker = basinmaker(path_working_folder=path_working_folder)
 
### define extent 
basinmaker.define_project_extent_method(
    mode="using_dem", path_dem_in=na_hydem
)

# divide domain in to subregions 
basinmaker.divide_domain_into_sub_regions_method(
    path_lakefile_in = in_lake,
    lake_attributes = ["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    Min_Num_Domain=9,
    Max_Num_Domain=30,
    Initaial_Acc=250000,
    Delta_Acc=20000,
    CheckLakeArea=10,
    fdr_path = '#',
    Acc_Thresthold_stream=2000,
    max_memory=2048*3,
    Out_Sub_Reg_Folder=Out_Sub_Reg_Dem_Folder,
)

End = timeit.default_timer()

print("use    ", start - End)
## -510.85583060000005 for woods lake watershed
