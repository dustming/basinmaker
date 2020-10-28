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
    """Transfer an input dbf file to dataframe
    
    Parameters
    ---------- 
    file_path   : string
    Full path to a shapefile 
    
    Returns:
    -------
    dataframe   : datafame 
    a pandas dataframe of attribute table of input shapefile    
    """
    tempinfo = Dbf5(file_path[:-3] + "dbf")
    dataframe = tempinfo.to_dataframe()
    return dataframe
    
def Return_Raster_As_Array(grassdb,grass_location,raster_mn):
    """Transfer an rater in grass database into np array
    Parameters
    ---------- 
    grassdb         : string
    Full path to a grass database 
    grass_location  : string
    location name in that grass database   
    raster_mn       : string
    raster name 
        
    Returns:
    -------
    Array            : array  
    np array of the raster. 
       
    """     
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location,create_opts='')
    Array = copy.deepcopy(garray.array(mapname=raster_mn))
    PERMANENT.close()
    return Array
        
def test_Generateinputdata():
    """test function that will: 
    Preprocessing input dataset, Function that used to project and clip 
    input dataset such as DEM, Land use, Lake polygon etc with defined 
    processing extent by function Generatmaskregion. And then it will
    rasterize these vector files.
        
    """
    ###Floder where store the inputs for tests function
    Data_Folder  = './testdata/Required_data_to_start_from_dem/'
    
    ###Folder where store the expected resuts
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_InDEM')
    
    ###Folder where the output will be generated 
    Temporary_Result_Folder_Result   = os.path.join('./testdata','Temporary_output_folder','testout2')
    Final_Result_Folder_Result       = os.path.join('./testdata','Final_output_folder','testout2')
    shutil.rmtree(Temporary_Result_Folder_Result,ignore_errors=True)
    
    ###Define path of input dataset
    Path_DEM_big           = os.path.join(Data_Folder, 'DEM_big_merit.tif')
    Path_DEM_small         = os.path.join(Data_Folder, 'DEM_samll_merit.tif')
    Path_Lake_ply          = os.path.join(Data_Folder, 'HyLake.shp')
    Path_bkf_wd            = os.path.join(Data_Folder, 'Bkfullwidth_depth.shp')
    Path_Landuse           = os.path.join(Data_Folder, 'landuse.tif')
    Path_Roughness_landuse = os.path.join(Data_Folder, 'Landuse.csv')

    ###Generate test resuts    
    RTtool=LRRT(dem_in = Path_DEM_small,WidDep = Path_bkf_wd,
               Lakefile = Path_Lake_ply,Landuse = Path_Landuse,
               Landuseinfo = Path_Roughness_landuse,
               OutputFolder = Final_Result_Folder_Result,
               TempOutFolder = Temporary_Result_Folder_Result,
               )
    RTtool.Generatmaskregion()
    RTtool.Generateinputdata()

    """Evaluate raster alllake 
       alllake is a lake raster which include all lakes
    """ 
    ### transfer expected raster alllake into np array Expected_alllake_Array   
    Expected_alllake_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'alllake')
    ### transfer resulted raster alllake into np array Result_alllake_Array
    Result_alllake_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'alllake')
    ### compare two Expected_alllake_Array and Result_alllake_Array 
    assert (Expected_alllake_Array == Result_alllake_Array).all()
    
    """Evaluate raster Lake_Bound 
       Lake_Bound is a raster represent the lake boundary grids  
    """     
    ### transfer expected raster Lake_Bound into np array Expected_LBound_Array       
    Expected_LBound_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'Lake_Bound')
    ### transfer resulted raster Lake_Bound into np array Result_LBound_Array
    Result_LBound_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'Lake_Bound')
    ### compare two Expected_LBound_Array and Result_LBound_Array                                                  
    assert (Expected_LBound_Array == Result_LBound_Array).all()
    
    """Evaluate raster acc_grass 
       acc_grass is the flow accumulation raster generated by 'r.watershed'  
    """   
    ### transfer expected raster acc_grass into np array Expected_acc_Array               
    Expected_acc_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'acc_grass')
    ### transfer resulted raster acc_grass into np array Result_acc_Array                                        
    Result_acc_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'acc_grass')
    ### compare two Expected_acc_Array and Result_acc_Array                                                  
    assert (Expected_acc_Array == Result_acc_Array).all()
    
    ### clean test output folder
    RTtool.Output_Clean()

        