import os

import pytest

from basinmaker import BasinMakerQGIS

#
# BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

###Floder where store the inputs for tests function
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")
path_output_folder = os.path.join("../../tests/testdata", "test3", "in_da")
path_working_folder = os.path.join("../../tests/testdata", "test3")

basinmaker = BasinMakerQGIS(
    path_output_folder=path_output_folder, path_working_folder=path_working_folder
)

basinmaker.simplify_routing_structure_by_drainage_area_method(
    OutputFolder=path_output_folder,
    Path_final_riv_ply=os.path.join(
        path_working_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_working_folder, "river_without_merging_lakes.shp"),
    Path_Con_Lake_ply=os.path.join(path_working_folder, "sl_connected_lake.shp"),
    Path_NonCon_Lake_ply=os.path.join(path_working_folder, "sl_non_connected_lake.shp"),
    Area_Min=100,
    gis_platform="qgis",
)

basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)
