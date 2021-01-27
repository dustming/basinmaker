import os
import numpy as np
import pytest
import tempfile
from basinmaker import basinmaker


#############################################
# define working folder, output folder amd data folder
#############################################
num = str(np.random.randint(1, 10000 + 1))
num = "100"
path_output_folder = os.path.join(
    tempfile.gettempdir(), "basinmaker_exp_" + num, "output"
)
path_working_folder = os.path.join(
    tempfile.gettempdir(), "basinmaker_exp_" + num, "work"
)
datafolder = os.path.join(
    os.path.abspath("../../../tests/testdata/"), "existing_lake_river_routing_structure"
)
HRU_Folder = os.path.join(os.path.abspath("../../../tests/testdata/"), "HRU")

#############################################
# initialize basinmaker with working folder
#############################################
basinmaker = basinmaker(path_working_folder=path_working_folder)


#############################################
# obtain part of the routing sturcture by gauge name or
# subbasin ids
#############################################
basinmaker.select_part_of_routing_product_methods(
    OutputFolder=path_output_folder,
    Path_Catchment_Polygon=os.path.join(
        datafolder, "catchment_without_merging_lakes.shp"
    ),
    Path_River_Polyline=os.path.join(datafolder, "river_without_merging_lakes.shp"),
    Path_Con_Lake_ply=os.path.join(datafolder, "sl_connected_lake.shp"),
    Path_NonCon_Lake_ply=os.path.join(datafolder, "sl_non_connected_lake.shp"),
    Gauge_NMS=["#"],  # ['02KB001'],
    mostdownid=1388,
    mostupid=1198,
    Path_Points="#",
    gis_platform="arcgis",
)
 
# basinmaker.combine_catchments_covered_by_the_same_lake_method(
#     OutputFolder=path_output_folder,
#     Path_final_rivply=os.path.join(
#         path_output_folder, "catchment_without_merging_lakes.shp"
#     ),
#     Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
#     gis_platform="qgis",
# )
