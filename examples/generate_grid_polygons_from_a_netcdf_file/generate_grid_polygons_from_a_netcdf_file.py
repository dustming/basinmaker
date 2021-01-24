import os
import numpy as np
import pytest
import tempfile
from basinmaker import basinmaker

#############################################
# define working folder, output folder amd data folder
#############################################
num = str(np.random.randint(1, 10000 + 1))
num = "101"
path_output_folder = os.path.join(
    tempfile.gettempdir(), "basinmaker_exp_" + num, "output"
)
path_working_folder = os.path.join(
    tempfile.gettempdir(), "basinmaker_exp_" + num, "work"
)
datafolder = os.path.join(
    "../../tests/testdata", "existing_lake_river_routing_structure"
)
HRU_Folder = os.path.join("../../tests/testdata/", "HRU")
demfolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")


#############################################
# initialize basinmaker with working folder
#############################################
basinmaker = basinmaker(path_working_folder=path_working_folder)


#############################################
# generate grid weights from netcdf file 
#############################################

basinmaker.obtain_grids_polygon_from_netcdf_file_method(
    netcdf_path=os.path.join(HRU_Folder,'test.nc'),
    output_folder=path_output_folder,
    coor_x_nm='lon',
    coor_y_nm='lat',
    is_rotated_grid=1,
    r_coor_x_nm='lon',
    r_coor_y_nm='lat',
    spatial_ref='EPSG:4326',
    x_add=-360,
    y_add=0,
    gis_platform='qgis',
)


basinmaker.generate_hrus_methods(
    OutputFolder=path_output_folder,
    Path_Subbasin_Ply=os.path.join(datafolder, "finalcat_info.shp"),
    Path_Connect_Lake_ply=os.path.join(datafolder, "sl_connected_lake.shp"),
    Path_Non_Connect_Lake_ply=os.path.join(datafolder, "sl_non_connected_lake.shp"),
    Lake_Id="Hylak_id",
    Path_Landuse_Ply="#",
    Landuse_ID="Landuse_ID",
    Path_Soil_Ply="#",
    Soil_ID="Soil_ID",
    Path_Veg_Ply="#",
    Veg_ID="Veg_ID",
    Path_Other_Ply_1="#",
    Other_Ply_ID_1="O_ID_1",
    Path_Other_Ply_2="#",
    Other_Ply_ID_2="O_ID_2",
    Landuse_info=os.path.join(HRU_Folder, "landuse_info.csv"),
    Soil_info=os.path.join(HRU_Folder, "soil_info.csv"),
    Veg_info=os.path.join(HRU_Folder, "veg_info.csv"),
    DEM=os.path.join(demfolder, "oih_30_dem.tif"),
    gis_platform="qgis",
)



basinmaker.generate_area_weight_of_two_polygons_method(
    target_polygon_path=os.path.join(path_output_folder, "finalcat_hru_info.shp"),
    mapping_polygon_path=os.path.join(path_output_folder, "Gridncply.shp"),
    col_nm="HRU_ID",
    output_folder=path_output_folder,
    gis_platform="qgis",
)

