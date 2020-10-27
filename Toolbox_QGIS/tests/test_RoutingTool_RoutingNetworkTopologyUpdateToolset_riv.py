import pytest
from ToolboxClass import LRRT
import os
import pandas as pd 
from simpledbf import Dbf5
import shutil 
import grass.script as grass
from grass.script import array as garray
import grass.script.setup as gsetup
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import raster as r
from grass.pygrass.modules import Module
from grass_session import Session
import copy 
import numpy as np 

os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='1'))

def Dbf_To_Dataframe(file_path):
    tempinfo = Dbf5(file_path[:-3] + "dbf")
    dataframe = tempinfo.to_dataframe()
    return dataframe
    
def Return_Raster_As_Array(grassdb,grass_location,raster_mn):    
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location,create_opts='')
    Array = copy.deepcopy(garray.array(mapname=raster_mn))
    PERMANENT.close()
    return Array
        
def test_RoutingNetworkTopologyUpdateToolset_riv():

    ###The second version of routing product 
    Data_Folder  = './testdata/Required_data_to_start_from_dem/'
    
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Result   = os.path.join('./testdata','Temporary_output_folder','testout')
    Final_Result_Folder_Result       = os.path.join('./testdata','Final_output_folder','testout')
    
    shutil.rmtree(Temporary_Result_Folder_Result,ignore_errors=True)
    ###The pathes for all inputs 
    Path_DEM_big           = os.path.join(Data_Folder, 'DEM_big_merit.tif')
    Path_DEM_small         = os.path.join(Data_Folder, 'DEM_samll_merit.tif')

    Path_Lake_ply          = os.path.join(Data_Folder, 'HyLake.shp')
    Path_bkf_wd            = os.path.join(Data_Folder, 'Bkfullwidth_depth.shp')
    Path_Landuse           = os.path.join(Data_Folder, 'landuse.tif')
    Path_Roughness_landuse = os.path.join(Data_Folder, 'Landuse.csv')
    ## run generate mask region using input dem  
    
    RTtool=LRRT(dem_in = Path_DEM_small,WidDep = Path_bkf_wd,
               Lakefile = Path_Lake_ply,Landuse = Path_Landuse,
               Landuseinfo = Path_Roughness_landuse,
               OutputFolder = Final_Result_Folder_Result,
               TempOutFolder = Temporary_Result_Folder_Result,
               )
    ### test using extent of input dem as processing extent 
    RTtool.Generatmaskregion()
    RTtool.Generateinputdata()
    RTtool.WatershedDiscretizationToolset(accthresold = 500)
    RTtool.AutomatedWatershedsandLakesFilterToolset(Thre_Lake_Area_Connect = 0,
                                                    Thre_Lake_Area_nonConnect = 0)
    RTtool.RoutingNetworkTopologyUpdateToolset_riv(projection = 'EPSG:3573')
#    RTtool.Define_Final_Catchment(OutputFolder = Final_Result_Folder_Expected,
#                                  Path_final_rivply = os.path.join(Final_Result_Folder_Expected,'finalriv_info_ply.shp'),
#                                  Path_final_riv    = os.path.join(Final_Result_Folder_Expected,'finalriv_info.shp'))
    
    Result_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Result,'finalriv_info_ply.shp')).sort_values(by=['SubId'])
    Result_Finalriv_info     = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Result,'finalriv_info.shp')).sort_values(by=['SubId'])
    Result_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Result,'Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    Result_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Result,'Non_Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    
    Expect_Finalriv_info_ply = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Expected,'finalriv_info_ply.shp')).sort_values(by=['SubId'])
    Expect_Finalriv_info     = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Expected,'finalriv_info.shp')).sort_values(by=['SubId'])
    Expect_Con_Lake_Ply      = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Expected,'Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    Expect_Non_Con_Lake_Ply  = Dbf_To_Dataframe(os.path.join(Final_Result_Folder_Expected,'Non_Con_Lake_Ply.shp')).sort_values(by=['Hylak_id'])
    

    assert len(Result_Finalriv_info_ply['SubId'].values) == len(Expect_Finalriv_info_ply['SubId'].values)
    assert np.allclose(Result_Finalriv_info_ply['BasArea'].values,Expect_Finalriv_info_ply['BasArea'].values,atol = 0.001)
    assert Result_Con_Lake_Ply.equals(Expect_Con_Lake_Ply)
    assert Result_Non_Con_Lake_Ply.equals(Expect_Non_Con_Lake_Ply)
    RTtool.Output_Clean()
    