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
        
def test_WatershedDiscretizationToolset():
    """test function that will: 
    Function that used to Generate a subbasin delineation and river 
    network using user provied flow accumulation thresthold 
    without considering lake. 
    """
    
    ###Floder where store the inputs for tests function
    Data_Folder  = './testdata/Required_data_to_start_from_dem/'
    
    ###Folder where store the expected resuts        
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_InDEM')

    ###Folder where the output will be generated     
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

    ###Generate test resuts    
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

    """Evaluate raster cat1 
    cat1 is a raster represent the delineated subbasins without
    considering lakes     
    """ 
    ### transfer expected raster cat1 into np array Expected_cat1_Array       
    Expected_cat1_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'cat1')
    ### transfer resulted raster cat1 into np array Result_cat1_Array
    Result_cat1_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'cat1')
    ### compare two Expected_Mask_Array and Result_Mask_Array 
    assert (Expected_cat1_Array == Result_cat1_Array).all()

    """Evaluate raster str_grass_r 
    str_grass_r is a river network in raster format     
    """ 
    ### transfer expected raster str_grass_r into np array Expected_str_Array       
    Expected_str_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'str_grass_r')
    ### transfer resulted raster str_grass_r into np array Result_str_Array
    Result_str_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'str_grass_r')
    ### compare two Expected_str_Array and Result_str_Array 
    assert (Expected_str_Array == Result_str_Array).all()

    """Evaluate raster Connect_Lake 
    Connect_Lake is the lake raster only contain lakes connected by 
    str_grass_r 
    """ 
    ### transfer expected raster Connect_Lake into np array Expected_CLake_Array           
    Expected_CLake_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'Connect_Lake')
    ### transfer resulted raster Connect_Lake into np array Result_CLake_Array
    Result_CLake_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'Connect_Lake')
    ### compare two Expected_CLake_Array and Result_CLake_Array 
    assert (Expected_CLake_Array == Result_CLake_Array).all()
    

    """Evaluate raster Nonconnect_Lake 
    Nonconnect_Lake is the lake raster only contain lakes not connected by 
    str_grass_r 
    """ 
    ### transfer expected raster Nonconnect_Lake into np array Expected_NCLake_Array                
    Expected_NCLake_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'Nonconnect_Lake')
    ### transfer resulted raster Nonconnect_Lake into np array Result_NCLake_Array
    Result_NCLake_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'Nonconnect_Lake')
    ### compare two Expected_NCLake_Array and Result_NCLake_Array 
    assert (Expected_NCLake_Array == Result_NCLake_Array).all()
    RTtool.Output_Clean()
     