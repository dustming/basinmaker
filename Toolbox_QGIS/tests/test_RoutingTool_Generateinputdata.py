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
        
def test_Generatmaskregion():

    ###The second version of routing product 
    Data_Folder  = './testdata/Required_data_to_start_from_dem/'
    
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Result   = os.path.join('./testdata','Temporary_output_folder','testout')
    Final_Result_Folder_Result       = os.path.join('./testdata','Final_output_folder','testout')
    
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
    
    Expected_Mask_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'alllake')
    Result_Mask_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'alllake')
    assert (Expected_Mask_Array == Result_Mask_Array).all()
    
    
    Expected_Mask_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'Lake_Bound')
    Result_Mask_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'Lake_Bound')
    assert (Expected_Mask_Array == Result_Mask_Array).all()
    
    
    Expected_Mask_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'acc_grass')
    Result_Mask_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'acc_grass')
    assert (Expected_Mask_Array == Result_Mask_Array).all()
    
        