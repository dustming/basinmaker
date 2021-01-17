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

# basinmaker.define_project_extent_method(
#     mode="using_dem", path_dem_in=os.path.join(datafolder, "DEM_samll_merit.tif")
# )

basinmaker.preprocessing_inputs_method(
    path_lakefile_in=os.path.join(datafolder, "HyLake.shp"),
    lake_attributes=["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    gis_platform="qgis",
)
