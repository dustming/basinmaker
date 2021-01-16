import os
import tempfile 
import numpy as np 
from basinmaker import basinmaker

#
# BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

###Floder where store the inputs for tests function
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")
path_output_folder = os.path.join("../../tests/testdata", "test3")
path_working_folder = os.path.join("../../tests/testdata", "test4")

basinmaker = basinmaker(
    path_output_folder=path_output_folder, path_working_folder=path_working_folder
)

basinmaker.define_project_extent_method(
    mode="using_outlet_pt",
    path_dem_in=os.path.join(datafolder, "DEM_big_merit.tif"),
    outlet_pt=[-91.977, 48.25],
    buffer_distance=0.03,
)

basinmaker.watershed_delineation_without_lake_method(
    acc_thresold=1000,
    mode="usingdem",
    max_memroy=1024 * 4,
    gis_platform="qgis",
)


basinmaker.watershed_delineation_add_lake_control_points(
    path_lakefile_in=os.path.join(datafolder, "HyLake.shp"),
    lake_attributes=["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    threshold_con_lake = 0,
    threshold_non_con_lake = 0,
    path_obsfile_in=os.path.join(datafolder, "obs.shp"),
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    max_memroy=1024 * 4,
    gis_platform="qgis",
)

basinmaker.add_attributes_to_catchments_method(
    path_bkfwidthdepth=os.path.join(datafolder, "Bkfullwidth_depth.shp"),
    bkfwd_attributes=["WIDTH", "DEPTH", "Q_Mean", "UP_AREA"],
    path_landuse=os.path.join(datafolder, "landuse.tif"),
    path_landuse_info=os.path.join(datafolder, "Landuse_info3.csv"),
    gis_platform="qgis",
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    lake_attributes =["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"] ,
    outlet_obs_id=1,
    path_sub_reg_outlets_v="#",
    output_folder=path_output_folder,
)


basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)
