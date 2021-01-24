import os
import numpy as np
import pytest
import tempfile
from basinmaker import basinmaker

#############################################
# define working folder, output folder amd data folder
#############################################
num = str(np.random.randint(1, 10000 + 1))
num = "101"
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
# generate grid weights from netcdf file 
#############################################

basinmaker.obtain_grids_polygon_from_netcdf_file_method(
    netcdf_path=os.path.join(HRU_Folder,'VIC_runoff.nc'),
    output_folder=path_output_folder,
    coor_x_nm='lon',
    coor_y_nm='lat',
    is_rotated_grid=1,
    r_coor_x_nm='rlon',
    r_coor_y_nm='rlat',
    spatial_ref='EPSG:4326',
    x_add=-360,
    y_add=0,
    gis_platform='qgis',
)