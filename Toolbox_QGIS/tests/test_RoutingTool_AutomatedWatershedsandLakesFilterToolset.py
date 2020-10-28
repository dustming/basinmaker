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
        
def test_AutomatedWatershedsandLakesFilterToolset():
    """test function that will: 
    Add lake inflow and outflow points as new subabsin outlet 
        
    """
    ###Floder where store the inputs for tests functions 
    Data_Folder  = './testdata/Required_data_to_start_from_dem/'
    
    ###Folder where store the expected resuts  
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_InDEM')
     
    ###Folder where the output will be generated 
    Temporary_Result_Folder_Result   = os.path.join('./testdata','Temporary_output_folder','testout1')
    Final_Result_Folder_Result       = os.path.join('./testdata','Final_output_folder','testout1')
    shutil.rmtree(Temporary_Result_Folder_Result,ignore_errors=True)
    
    ###Define path of input dataset
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
    RTtool.WatershedDiscretizationToolset(accthresold = 500)
    RTtool.AutomatedWatershedsandLakesFilterToolset(Thre_Lake_Area_Connect = 0,
                                                    Thre_Lake_Area_nonConnect = 0)
    
    """Evaluate raster Net_cat 
       Net_cat is a subbsin raster after adding lake inlet and outlet as 
       additional subbasin outlet.      
    """
    ### transfer expected raster Net_cat into np array Expected_Net_cat_Array
    Expected_Net_cat_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'Net_cat')
    
    ### transfer test raster Net_cat into np array Result_Net_cat_Array
    Result_Net_cat_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'Net_cat')
    ### compare two Expected_Net_cat_Array and Result_Net_cat_Array 
    assert (Expected_Net_cat_Array == Result_Net_cat_Array).all()

    """Evaluate raster ndir_grass 
       ndir_grass is a raster represent modified flow direction 
    """
    ### transfer expected raster ndir_grass into np array Expected_ndir_Array
    Expected_ndir_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'ndir_grass')
    ### transfer test raster ndir_grass into np array Result_ndir_Array
    Result_ndir_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                  grass_location = 'Geographic',
                                                  raster_mn = 'ndir_grass')                                          
    ### compare two Expected_ndir_Array and Result_ndir_Array 
    assert (Expected_ndir_Array == Result_ndir_Array).all()
    
    ### clean test output folder 
    RTtool.Output_Clean()