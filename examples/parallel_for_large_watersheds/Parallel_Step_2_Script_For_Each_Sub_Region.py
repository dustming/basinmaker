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


import os
import shutil
import sys
import tempfile

import pandas as pd

from basinmaker import basinmaker

############ Variable needs to be modified to run this example ######

### Define a output folder where to store subregion information
Outputfolder = "C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Petawawa/lake_of_woods/"

### The BasinMaker folder
datafolder = "../../tests/testdata/Required_data_to_start_from_dem"

########### Variable needs to be modified to run this example ######
Out_Sub_Reg_Dem_Folder = os.path.join(Outputfolder, "SubRegion_info")

### the id of each subregion
ith_subregion = int(sys.argv[1])


### Inputs
na_hydem = os.path.join(Out_Sub_Reg_Dem_Folder, "sub_reg_dem.pack") 

### open log file
if ith_subregion == 0:
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder, "Log"), "w")
    file_object.write("SubRegion Id,           Status" + "\n")
else:
    file_object = open(os.path.join(Out_Sub_Reg_Dem_Folder, "Log"), "a")


SubReg_info = pd.read_csv(os.path.join(Out_Sub_Reg_Dem_Folder, "Sub_reg_info.csv"))


### defineinputs for ith subregion
basinid = int(SubReg_info["Sub_Reg_ID"].values[ith_subregion])
ProjectName = SubReg_info["ProjectNM"].values[ith_subregion]
Path_Sub_Polygon = os.path.join(
    Out_Sub_Reg_Dem_Folder, SubReg_info["Ply_Name"].values[ith_subregion]
)

path_working_folder = os.path.join("C:/Users/dustm/Documents","testsub",ProjectName)

path_output_folder = os.path.join(Out_Sub_Reg_Dem_Folder,'subregion_result',ProjectName)

### run for each sub region

### initialize the toolbox
basinmaker = basinmaker(
   path_working_folder = path_working_folder
)

basinmaker.define_project_extent_method(
    mode="using_provided_ply", path_dem_in=na_hydem,path_extent_ply=Path_Sub_Polygon
)

basinmaker.watershed_delineation_without_lake_method(
    mode="usingsubreg",
    max_memroy=1024 * 4,
    subreg_acc_path = os.path.join(
        Out_Sub_Reg_Dem_Folder, "sub_reg_acc.pack"
    ),
    subreg_fdr_path = os.path.join(
        Out_Sub_Reg_Dem_Folder, "sub_reg_nfdr_arcgis.pack"
    ),
    subreg_str_r_path = os.path.join(
        Out_Sub_Reg_Dem_Folder, "sub_reg_str_r.pack"
    ),
    subreg_str_v_path = os.path.join(
        Out_Sub_Reg_Dem_Folder, "sub_reg_str_v.pack"
    ),
    gis_platform="qgis",
)

basinmaker.watershed_delineation_add_lake_control_points(
    path_lakefile_in=os.path.join(datafolder, "hylake.shp"),
    lake_attributes=["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    threshold_con_lake = 0,
    threshold_non_con_lake = 0,
    path_obsfile_in=os.path.join(datafolder, "obs.shp"),
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    path_sub_reg_outlets_v = os.path.join(
        Out_Sub_Reg_Dem_Folder, "Sub_Reg_Outlet_v.pack"
    ),
    max_memroy=1024 * 4,
    gis_platform="qgis",
)

basinmaker.add_attributes_to_catchments_method(
    path_bkfwidthdepth=os.path.join(datafolder, "bkf_wd.shp"),
    bkfwd_attributes=["WIDTH", "DEPTH", "Q_Mean", "UP_AREA"],
    path_landuse=os.path.join(datafolder, "landuse_modis_250.tif"),
    path_landuse_info=os.path.join(datafolder, "Landuse_info3.csv"),
    gis_platform="qgis",
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    lake_attributes =["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"] ,
    outlet_obs_id=basinid,
    output_folder=path_output_folder,
    path_sub_reg_outlets_v = os.path.join(
        Out_Sub_Reg_Dem_Folder, "outlet_pt_info.shp"
    ),
    k_in = SubReg_info["k"].values[ith_subregion],
    c_in = SubReg_info["c"].values[ith_subregion],
)
#
basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)

if os.path.exists(
    os.path.join(path_output_folder, "finalcat_info.shp")
):
    file_object.write(str(basinid) + "          Successful " + "\n")
    file_object.close()
else:
    shutil.rmtree(path_output_folder, ignore_errors=True)
    file_object.write(str(basinid) + "          Failed " + "\n")
    file_object.close()
