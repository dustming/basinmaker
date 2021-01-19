import os
import numpy as np 
import pytest
import tempfile
from basinmaker import basinmaker

#############################################
# define working folder, output folder amd data folder  
#############################################
num  = str(np.random.randint(1, 10000 + 1))
path_output_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_" +num,"output")
path_working_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_" +num,"work")
datafolder = os.path.join("../../tests/testdata", "existing_lake_river_routing_structure")
HRU_Folder = os.path.join("../../tests/testdata/", "HRU")


#############################################
# initialize basinmaker with working folder    
#############################################
basinmaker = basinmaker(path_working_folder=path_working_folder)


#############################################
# obtain simplified catchments and river network before merging lakes  
#############################################
basinmaker.simplify_routing_structure_by_drainage_area_method(
    OutputFolder=path_output_folder,
    Path_final_riv_ply=os.path.join(
        datafolder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(datafolder, "river_without_merging_lakes.shp"),
    Path_Con_Lake_ply=os.path.join(datafolder, "sl_connected_lake.shp"),
    Path_NonCon_Lake_ply=os.path.join(datafolder, "sl_non_connected_lake.shp"),
    Area_Min=500,
    gis_platform="qgis",
)


#############################################
# combine catchments covered by the same lake 
#############################################
basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)

#############################################
# generate HRUS, lake and land hru  
#############################################
# basinmaker.generate_hrus_methods(
#     OutputFolder=path_output_folder,
#     Path_Subbasin_Ply=os.path.join(
#         path_output_folder, "finalcat_info.shp"
#     ),
#     Path_Connect_Lake_ply=os.path.join(
#         path_output_folder, "sl_connected_lake.shp"
#     ),
#     Path_Non_Connect_Lake_ply=os.path.join(
#         path_output_folder, "sl_non_connected_lake.shp"
#     ),
#     Lake_Id="Hylak_id",    
#     Path_Landuse_Ply="#",
#     Landuse_ID="Landuse_ID",
#     Path_Soil_Ply="#",
#     Soil_ID="Soil_ID",
#     Path_Veg_Ply="#",
#     Veg_ID="Veg_ID",
#     Path_Other_Ply_1="#",
#     Other_Ply_ID_1="O_ID_1",
#     Path_Other_Ply_2="#",
#     Other_Ply_ID_2="O_ID_2",
#     Landuse_info=os.path.join(HRU_Folder, "landuse_info.csv"),
#     Soil_info=os.path.join(HRU_Folder, "soil_info.csv"),
#     Veg_info=os.path.join(HRU_Folder, "veg_info.csv"),
#     DEM="#",
#     gis_platform="qgis",
# )