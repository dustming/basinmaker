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

##### Define subregion 
ProjectName = 'mid'  #'Example2' #'Example'
OutletPointxy = [-91.961,48.975] # [-92.387,49.09] #mid [-92.387,49.09] #[-91.961,48.975]  ##small[-92.387,49.09]	
RTtool=LRRT(dem_in = na_hydem, Lakefile = in_lake,obspoint = in_obs,OutputFolder = Outputfolder,ProjectNM = ProjectName,Path_Sub_Reg_Out_Folder = Out_Sub_Reg_Dem_Folder)
RTtool.Generatmaskregion(OutletPoint = OutletPointxy)
RTtool.Generatesubdomain(Min_Num_Domain = 9,Max_Num_Domain = 13,Initaial_Acc = 10000,Delta_Acc = 2000,max_memory=2048,Acc_Thresthold_stream=50,CheckLakeArea = 1)
RTtool.Generatesubdomainmaskandinfo(Out_Sub_Reg_Dem_Folder = Out_Sub_Reg_Dem_Folder)
