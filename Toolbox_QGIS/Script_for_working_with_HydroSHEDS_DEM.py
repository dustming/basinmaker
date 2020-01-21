# Copyright 2018-2018 Ming Han - ming.han(at)uwaterloo.ca
#
# License
# This file is part of Ming Han's personal code library.
#
# Ming Han's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation+ "     " +  either version 3 of the License+ "     " +  or
# (at your option) any later version.
#
# Ming Han's personal code library is distributed in the hope that it will be useful+ "     " + 
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Ming Han's personal code library.  If not+ "     " +  see <http://www.gnu.org/licenses/>.
#
###########################################################################################
#
#This is a python script to use An automated ArcGIS toolbox for watershed delineation with lakes v1.tbx
#This script show how to use this tool box to generate lake-river routing sturcture with HydroSHEDS dataset.
#
###########################################################################################
Scriptfolder = "C:/Users/dustm/Documents/GitHub/RoutingTool/QGIS_version/"

import os
import sys 
sys.path.insert(1, "C:/Users/dustm/Documents/GitHub/RoutingTool/QGIS_version/")
import Generateinputs
########### Define datafolder and output folder
#Datafolder where you save the HydroSHEDS data 
#Outputfolder wihere you save the generated lake-river routing structure
###########

Datafolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Database/"
Outputfolder = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Outputs1/"


QGISpython = "C:/Python27/ArcGISx6410.6/python.exe"
Scriptfolder = "C:/Users/dustm/Documents/GitHub/RoutingTool/QGIS_version/"
###Import the toolbox

####################################### Delineate catchment without lake

###########Prepare Inputs 
na_hyshedply = Datafolder + "Shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev12_v1c.shp"
na_hydem = Datafolder + "Rasters/na_dem_15s/na_dem_15s"
na_dir = Datafolder + "Rasters/na_dir_15s/na_dir_15s"
na_acc = Datafolder + "Rasters/na_acc_15s/na_acc_15s"
in_wd = Datafolder + "Shapefiles/Width_Depth/narivs_new.shp"
na_shddbf = Datafolder +  "Shapefiles/hybas_na_lev01-12_v1c/hybas_na_lev12_v1c.dbf"
in_lake = Datafolder +  "Shapefiles/HyLake.shp"
in_lakearea = "0"
in_obs = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Data/POI.shp"
OutHyID2 = "-1"
landuse = Datafolder + "Rasters/landuse"
landuseinfo = Datafolder + "Rasters/landuse.csv"

accthres = "500"
OutHyID1 = 7120968990
#OutHyID1 = 7120298860
#OutHyID1 = 7120528350
basename = 'RoutingProductExample123' 
Outputfoldercase = Outputfolder+basename #+'_ACC_'+accthres
if not os.path.exists(Outputfoldercase):
	os.makedirs(Outputfoldercase)

step1command = QGISpython + "   " + Scriptfolder + "Generateinputs.py" +  "    " + na_hydem+ "     " + na_dir+ "     " + na_acc+ "     " + na_hyshedply+ "     " + in_wd+ "     " + in_lake+ "     " + landuse+ "     " + landuseinfo+ "     " + in_obs+ "     " + str(OutHyID1)+ "     " + OutHyID2+ "     " + Outputfoldercase
