import os
import tempfile 
import numpy as np 
from basinmaker import basinmaker

#############################################
# define working folder, output folder amd data folder  
#############################################
num  = str(np.random.randint(1, 10000 + 1))
path_output_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_oih" +num,"output")
path_working_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_oih" +num,"work")
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")

#############################################
# initialize basinmaker with working folder    
#############################################
basinmaker = basinmaker(
    path_working_folder=path_working_folder
)


############################################
# define extent of the processing domain  
############################################
basinmaker.define_project_extent_method(
    mode="using_dem",
    path_dem_in=os.path.join(datafolder, "oih_30_dem.tif"),
    gis_platform="qgis",
)


#############################################
# generate a watershed delineation without considering lakes 
#############################################
basinmaker.watershed_delineation_without_lake_method(
    acc_thresold=5000,
    mode="usingdem",
    max_memroy=1024 * 4,
    gis_platform="qgis",
)


############################################
# add lake and obs control points in the existing watershed delineation  
############################################
basinmaker.watershed_delineation_add_lake_control_points(
    path_lakefile_in=os.path.join(datafolder, "hylake.shp"),
    lake_attributes=["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    threshold_con_lake = 0,
    threshold_non_con_lake = 0,
    path_obsfile_in=os.path.join(datafolder, "obs.shp"),
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    max_memroy=1024 * 4,
    gis_platform="qgis",
)


#############################################
# add hydrological attributes to existing watershed delineation  
#############################################
basinmaker.add_attributes_to_catchments_method(
    path_bkfwidthdepth=os.path.join(datafolder, "bkf_wd.shp"),
    bkfwd_attributes=["WIDTH", "DEPTH", "Q_Mean", "UP_AREA"],
    path_landuse=os.path.join(datafolder, "landuse_modis_250.tif"),
    path_landuse_info=os.path.join(datafolder, "Landuse_info3.csv"),
    gis_platform="qgis",
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    lake_attributes =["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"] ,
    outlet_obs_id=1,
    path_sub_reg_outlets_v="#",
    output_folder=path_output_folder,
)


#############################################
# combine catchments covered by the same lakes 
#############################################
basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)
