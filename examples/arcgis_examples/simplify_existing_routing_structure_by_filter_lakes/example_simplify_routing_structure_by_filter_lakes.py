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
# simplify existing routing structure by filtering lakes
#############################################
basinmaker.simplify_routing_structure_by_filter_lakes_method(
    OutputFolder=path_output_folder,
    Path_final_riv_ply=os.path.join(datafolder, "catchment_without_merging_lakes.shp"),
    Path_final_riv=os.path.join(datafolder, "river_without_merging_lakes.shp"),
    Path_Con_Lake_ply=os.path.join(datafolder, "sl_connected_lake.shp"),
    Path_NonCon_Lake_ply=os.path.join(datafolder, "sl_non_connected_lake.shp"),
    Thres_Area_Conn_Lakes=5,
    Thres_Area_Non_Conn_Lakes=1,
    Selection_Method="ByArea",
    Selected_Lake_List_in=[],
    gis_platform="arcgis",
)

basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="arcgis",
)
