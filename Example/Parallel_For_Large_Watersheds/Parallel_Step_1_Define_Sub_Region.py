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

from ToolboxClass import LRRT

############ Variable needs to be modified to run this example ######

### Define a output folder where to store subregion information
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"

### The BasinMaker folder
BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

########### Variable needs to be modified to run this example ######


### Define derived folder
DataBase_Folder = os.path.join(
    BasinMaker_Folder,
    "Toolbox_QGIS",
    "tests",
    "testdata",
    "Required_data_to_start_from_dem",
)
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder, "SubRegion_info")

### Start timer
start = timeit.default_timer()

### Define input paths
na_hydem = os.path.join(DataBase_Folder, "DEM_big_merit.tif")  #'HydroSHED15S.tif')#
in_wd = os.path.join(DataBase_Folder, "Bkfullwidth_depth.shp")
in_lake = os.path.join(DataBase_Folder, "HyLake.shp")
in_obs = os.path.join(DataBase_Folder, "obs.shp")
landuse = os.path.join(DataBase_Folder, "landuse.tif")
landuseinfo = os.path.join(DataBase_Folder, "Landuse_info.csv")


### Initialize the BasinMaker
RTtool = LRRT(
    dem_in=na_hydem,
    Lakefile=in_lake,
    OutputFolder=Outputfolder,
    Path_Sub_Reg_Out_Folder=Out_Sub_Reg_Dem_Folder,
    Is_Sub_Region=-1,
)
### Define Region of Interest
RTtool.Generatmaskregion()
### Define sub-region
RTtool.Generatesubdomain(
    Min_Num_Domain=1,
    Max_Num_Domain=1000000000,
    Initaial_Acc=500000,
    Delta_Acc=500000,
    max_memory=2048 * 8,
    Acc_Thresthold_stream=2000,
    CheckLakeArea=10,
)
### Generate subregion output and define routing structure between subregions
RTtool.Generatesubdomainmaskandinfo(Out_Sub_Reg_Dem_Folder=Out_Sub_Reg_Dem_Folder)

End = timeit.default_timer()

print("use    ", start - End)
## -510.85583060000005 for woods lake watershed
