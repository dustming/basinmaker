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
#This is a python script to use An automated ArcGIS toolbox for watershed delineation with lakes v1.tbx
#This script show how to use this tool box to generate lake-river routing sturcture with HydroSHEDS dataset.
#
###########################################################################################

from ToolboxClass import LRRT
import os
import pandas as pd 
import tempfile 
import shutil
import sys


########### Define datafolder and output folder
#Datafolder where you save the HydroSHEDS data 
#Outputfolder wihere you save the generated lake-river routing structure
###########

Datafolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Data_Lakeofwoods/"
Outputfolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Outputs1/"
####

ith_subregion = int(sys.argv[1])

###########Define Inputs the same for each region
na_hyshedply = "#"
na_hydem = Datafolder + "merit_dem.tif" #"Woodlake_example.tif" #"merit_dem.tif"
na_dir = "#"
na_acc = "#"
in_wd = Datafolder + "narivs_new.shp"
in_lake = Datafolder +  "HyLake2.shp"
in_obs = Datafolder + "obsfinal.shp"
landuse = Datafolder + "landuse2"
landuseinfo = Datafolder + "landuse.csv"
CA_HYDAT =os.path.join(Datafolder,'Hydat.sqlite3')
Out_Sub_Reg_Dem_Folder = os.path.join(Datafolder,'SubRegion') 


### open log file 
if ith_subregion == 0:
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder,'Log'), 'w')
    file_object.write('SubRegion Id,           Status' + '\n')
else:
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder,'Log'), 'a')


SubReg_info = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder,'Sub_reg_info.csv'))

	##defineinputs for ith subregion
basinid = int(SubReg_info['Sub_Reg_ID'].values[ith_subregion])            
ProjectName = SubReg_info['ProjectNM'].values[ith_subregion]
Path_Sub_Polygon = os.path.join(Out_Sub_Reg_Dem_Folder,SubReg_info['Ply_Name'].values[ith_subregion])
    
    ## run for each sub region
try:
    RTtool=LRRT(Lakefile = in_lake,Landuse = landuse,Landuseinfo = landuseinfo,obspoint = in_obs,OutputFolder = Outputfolder, ProjectNM = ProjectName,Path_Sub_Reg_Out_Folder = Out_Sub_Reg_Dem_Folder,Is_Sub_Region = 1)
    RTtool.Generatmaskregion(Path_Sub_Polygon = Path_Sub_Polygon)
    RTtool.Generateinputdata()
    RTtool.WatershedDiscretizationToolset(max_memroy = 1024*2)
    RTtool.AutomatedWatershedsandLakesFilterToolset(Thre_Lake_Area_Connect = 5,Thre_Lake_Area_nonConnect = 5,max_memroy = 1024*2)
    RTtool.RoutingNetworkTopologyUpdateToolset_riv('EPSG:3573',Outlet_Obs_ID = basinid)
    RTtool.Output_Clean(clean = 'False')

    Datafolder_final = os.path.join(Outputfolder,ProjectName)
    RTtool.Define_Final_Catchment(Datafolder = Datafolder_final)
    file_object.write(str(basinid) + '          Successful ' + '\n')
    file_object.close()
except:
    shutil.rmtree(os.path.join(Outputfolder,ProjectName),ignore_errors=True)
    print(basinid)
    file_object.write(str(basinid) + '          Failed ' + '\n')
    file_object.close()
    pass

