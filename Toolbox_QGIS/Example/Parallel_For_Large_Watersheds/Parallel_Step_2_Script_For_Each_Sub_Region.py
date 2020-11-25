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
# This is a python script to use BasinMaker to define lake river routing structure for  
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



############ Variable needs to be modified to run this example ######     

### Define a output folder where to store subregion information 
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"

### The BasinMaker folder 
BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

########### Variable needs to be modified to run this example ######


### Define derived folder 
DataBase_Folder = os.path.join(BasinMaker_Folder,'Toolbox_QGIS','tests','testdata','Required_data_to_start_from_dem')
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder,'SubRegion_info')

### the id of each subregion 
ith_subregion = int(sys.argv[1])


### Inputs
na_hydem = os.path.join(DataBase_Folder,'DEM_big_merit.tif')#'HydroSHED15S.tif')#
in_wd = os.path.join(DataBase_Folder,'Bkfullwidth_depth.shp')
in_lake = os.path.join(DataBase_Folder,'HyLake.shp')
in_obs = os.path.join(DataBase_Folder,'obs.shp')
landuse = os.path.join(DataBase_Folder,'landuse.tif')
landuseinfo = os.path.join(DataBase_Folder,'Landuse_info.csv')
#CA_HYDAT =os.path.join(DataBase_Folder,'Shapefiles','Observation','Hydat.sqlite3')

### open log file
if ith_subregion == 0:
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder,'Log'), 'w')
    file_object.write('SubRegion Id,           Status' + '\n')
else:
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder,'Log'), 'a')


SubReg_info = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv'))


### defineinputs for ith subregion
basinid = int(SubReg_info['Sub_Reg_ID'].values[ith_subregion])
ProjectName = SubReg_info['ProjectNM'].values[ith_subregion]
Path_Sub_Polygon = os.path.join(Out_Sub_Reg_Dem_Folder,SubReg_info['Ply_Name'].values[ith_subregion])

### run for each sub region

### initialize the toolbox 
RTtool=LRRT(Lakefile = in_lake,Landuse = landuse,Landuseinfo = landuseinfo,
            obspoint = in_obs,OutputFolder = os.path.join(Outputfolder,ProjectName),
            Path_Sub_Reg_Out_Folder = Out_Sub_Reg_Dem_Folder,Is_Sub_Region = 1,
            WidDep = in_wd,debug = True,TempOutFolder = os.path.join(tempfile.gettempdir(),ProjectName))

### define processing extent    
RTtool.Generatmaskregion(Path_Sub_Polygon = Path_Sub_Polygon)

### prepare input data within the processing extent     
RTtool.Generateinputdata()

### delineate watershed without lakes 
RTtool.WatershedDiscretizationToolset(max_memroy = 1024*5)

### add lakes into routing network
RTtool.AutomatedWatershedsandLakesFilterToolset(Thre_Lake_Area_Connect = 0,Thre_Lake_Area_nonConnect = 0,max_memroy = 1024*2)

### calcuate parameters 
RTtool.RoutingNetworkTopologyUpdateToolset_riv(projection = 'EPSG:3573',Outlet_Obs_ID = basinid)

### Megrge catchment covered by the same lake 
RTtool.Define_Final_Catchment(OutputFolder = os.path.join(Outputfolder,ProjectName),Path_final_rivply = os.path.join(Outputfolder,ProjectName,'finalriv_info_ply.shp'),Path_final_riv = os.path.join(Outputfolder,ProjectName,'finalriv_info.shp'))

if os.path.exists(os.path.join(os.path.join(Outputfolder,ProjectName),'finalcat_info.shp')):
    file_object.write(str(basinid) + '          Successful ' + '\n')
    file_object.close()
else:
    shutil.rmtree(os.path.join(Outputfolder,ProjectName),ignore_errors=True)
    file_object.write(str(basinid) + '          Failed ' + '\n')
    file_object.close()
