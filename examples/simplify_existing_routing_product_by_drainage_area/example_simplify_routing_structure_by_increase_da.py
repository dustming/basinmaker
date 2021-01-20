import os
import numpy as np
import pytest
import tempfile
from basinmaker import basinmaker

#############################################
# define working folder, output folder amd data folder
#############################################
num = str(np.random.randint(1, 10000 + 1))
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


#############################################
# initialize basinmaker with working folder
#############################################
basinmaker = basinmaker(path_working_folder=path_working_folder)


#############################################
# obtain simplified catchments and river network before merging lakes
#############################################
basinmaker.simplify_routing_structure_by_drainage_area_method(
    OutputFolder=path_output_folder,
    Path_final_riv_ply=os.path.join(datafolder, "catchment_without_merging_lakes.shp"),
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
