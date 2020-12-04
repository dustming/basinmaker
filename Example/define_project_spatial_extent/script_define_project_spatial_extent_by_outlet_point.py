import os

import pytest

from basinmaker import BasinMakerQGIS

#
# BasinMaker_Folder = "C:/Users/dustm/Documents/GitHub/RoutingTool"

###Floder where store the inputs for tests function
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")
path_output_folder = os.path.join("../../tests/testdata", "test1")
path_working_folder = os.path.join("../../tests/testdata", "test2")

basinmaker = BasinMakerQGIS(
    path_output_folder=path_output_folder, path_working_folder=path_working_folder
)

basinmaker.define_project_extent_method(
    mode="using_outlet_pt",
    path_dem_in=os.path.join(datafolder, "DEM_big_merit.tif"),
    outlet_pt=[-92.387, 49.09],
    buffer_distance=0.1,
)
