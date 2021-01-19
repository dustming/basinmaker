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
demfolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")

HRU_Folder = os.path.join("../../tests/testdata/", "HRU")
print(path_output_folder)

basinmaker = basinmaker(
    path_output_folder=path_output_folder, path_working_folder=path_working_folder
)

basinmaker.generate_hrus_methods(
    OutputFolder=path_output_folder,
    Path_Subbasin_Ply=os.path.join(
        datafolder, "finalcat_info.shp"
    ),
    Path_Connect_Lake_ply=os.path.join(
        datafolder, "sl_connected_lake.shp"
    ),
    Path_Non_Connect_Lake_ply=os.path.join(
        datafolder, "sl_non_connected_lake.shp"
    ),
    Lake_Id="Hylak_id",    
    Path_Landuse_Ply="#",
    Landuse_ID="Landuse_ID",
    Path_Soil_Ply="#",
    Soil_ID="Soil_ID",
    Path_Veg_Ply="#",
    Veg_ID="Veg_ID",
    Path_Other_Ply_1="#",
    Other_Ply_ID_1="O_ID_1",
    Path_Other_Ply_2="#",
    Other_Ply_ID_2="O_ID_2",
    Landuse_info=os.path.join(HRU_Folder, "landuse_info.csv"),
    Soil_info=os.path.join(HRU_Folder, "soil_info.csv"),
    Veg_info=os.path.join(HRU_Folder, "veg_info.csv"),
    DEM=os.path.join(demfolder, "oih_30_dem.tif"),
    gis_platform="qgis",
)

basinmaker.generate_raven_model_inputs_method(
    path_final_hru_info=os.path.join(path_output_folder,'finalcat_hru_info.shp'),
    startyear=2010,
    endYear=2014,
    CA_HYDAT="#",
    warmup=1,
    template_folder="#",
    lake_as_gauge=False,
    writeobsrvt=False,
    downloadobsdata=False,
    model_name="test",
    subbasingroup_nm_channel=["Allsubbasins"],
    subbasingroup_length_channel=[-1],
    subbasingroup_nm_lake=["AllLakesubbasins"],
    subbasingroup_area_lake=[-1],
    outputfolder=path_output_folder,
    forcing_input_file="#",
)