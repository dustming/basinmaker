import os

import pytest

from basinmaker import BasinMakerQGIS

#
# BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

###Floder where store the inputs for tests function
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")
path_output_folder = os.path.join("../../tests/testdata", "test3")
path_working_folder = os.path.join("../../tests/testdata", "test4")

basinmaker = BasinMakerQGIS(
    path_output_folder=path_output_folder, path_working_folder=path_working_folder
)

basinmaker.define_project_extent_method(
    mode="using_hybasin",
    path_dem_in=os.path.join(datafolder, "DEM_big_merit.tif"),
    buffer_distance=0.1,
    hybasin_ply="C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Geo_Data_Base/Shapefiles/HydroBASINS/hybas_na_lev07_v1c.shp",
    down_hybasin_id=7070317700,
)
