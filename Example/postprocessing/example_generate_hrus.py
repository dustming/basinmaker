import os

import pytest

from basinmaker import BasinMakerQGIS

#
# BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

###Floder where store the inputs for tests function
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")
path_output_folder = os.path.join("../../tests/testdata", "test3", "select_product")
path_working_folder = os.path.join("../../tests/testdata", "test3")
HRU_Folder = os.path.join("../../tests/testdata/", "HRU")

basinmaker = BasinMakerQGIS(
    path_output_folder=path_output_folder, path_working_folder=path_working_folder
)

basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_working_folder,
    Path_final_rivply=os.path.join(
        path_working_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_working_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)

basinmaker.generate_hrus_methods(
    OutputFolder=path_working_folder,
    Path_Subbasin_Ply=os.path.join(
        path_working_folder, "finalcat_info.shp"
    ),
    Path_Connect_Lake_ply=os.path.join(
        path_working_folder, "sl_connected_lake.shp"
    ),
    Path_Non_Connect_Lake_ply=os.path.join(
        path_working_folder, "sl_non_connected_lake.shp"
    ),
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
    DEM="#",
    gis_platform="qgis",
)

