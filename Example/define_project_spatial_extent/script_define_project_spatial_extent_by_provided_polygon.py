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
    mode="using_provided_ply",
    #    path_dem_in = os.path.join(datafolder, "Sub_Reg_dem.pack"),
    path_dem_in=os.path.join(datafolder, "DEM_big_merit.tif"),
    path_extent_ply=os.path.join(datafolder, "HyMask_region_80011_nobuffer.shp"),
)
