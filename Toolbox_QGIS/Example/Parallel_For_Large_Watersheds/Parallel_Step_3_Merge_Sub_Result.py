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


from ToolboxClass import LRRT
import os
import pandas as pd
import tempfile
import shutil
import sys
import tempfile
import timeit


############ Variable needs to be modified to run this example ######     

### Define a output folder where to store subregion information 
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"

### The BasinMaker folder 
BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

############ Variable needs to be modified to run this example ######  


### Define derived folder 
DataBase_Folder = os.path.join(BasinMaker_Folder,'Toolbox_QGIS','tests','testdata','Required_data_to_start_from_dem')
Sub_Region_OutputFolder = Outputfolder
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder,'SubRegion_info')
OutputFolder_Combined = os.path.join(Outputfolder,'Combined')

####
SubReg_info = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv'))

start = timeit.default_timer()

RTtool=LRRT()
RTtool.Combine_Sub_Region_Results(Sub_Region_info = SubReg_info,Sub_Region_OutputFolder = Sub_Region_OutputFolder,
                                  Path_Down_Stream_Points = os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_Reg_Outlets_Down.shp'),
                                  Is_Final_Result = True,
                                  OutputFolder = OutputFolder_Combined)
End = timeit.default_timer()

print("use    ",start - End)
### -73.6744739