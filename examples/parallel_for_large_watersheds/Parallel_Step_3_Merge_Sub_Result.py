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
# This is a python script to use BasinMaker to combine lake river routing structure of
# each subregion
#
###########################################################################################


import os
import shutil
import sys
import tempfile
import timeit

import pandas as pd

from basinmaker import basinmaker

############ Variable needs to be modified to run this example ######

### Define a output folder where to store subregion information
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"

########### Variable needs to be modified to run this example ######
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder, "SubRegion_info")
Sub_Region_OutputFolder = os.path.join(Out_Sub_Reg_Dem_Folder, "subregion_result")
OutputFolder_Combined = os.path.join(Out_Sub_Reg_Dem_Folder, "Combined")
path_subregion_inlet = os.path.join(Out_Sub_Reg_Dem_Folder, "sub_reg_inlet.shp")
####
SubReg_info = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_reg_info.csv"))
path_working_folder = Out_Sub_Reg_Dem_Folder

start = timeit.default_timer()
### initialize the toolbox
basinmaker = basinmaker(path_working_folder=Out_Sub_Reg_Dem_Folder)

basinmaker.combine_sub_region_results_method(
    path_sub_region_info=SubReg_info,
    sub_region_outputfolder=Sub_Region_OutputFolder,
    outputfolder=OutputFolder_Combined,
    is_final_result=False,
    path_subregion_inlet=path_subregion_inlet,
    gis_platform="qgis",
)

End = timeit.default_timer()

print("use    ", start - End)
### -73.6744739
