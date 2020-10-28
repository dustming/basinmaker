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
        
def test_Generatmaskregion():
    """test function that will: 
    Define the processing extent. Function that used to define 
    processing spatial extent (PSE). The processing spatial extent 
    is a region where Toolbox will work in. Toolbox will not 
    process grids or features outside the processing spatial extent.
    """
    ###Floder where store the inputs for tests function
    Data_Folder  = './testdata/Required_data_to_start_from_dem/'
    
    ###Folder where store the expected resuts    
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_InDEM')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_InDEM')
    
    ###Folder where the output will be generated 
    Temporary_Result_Folder_Result   = os.path.join('./testdata','Temporary_output_folder','testout3')
    Final_Result_Folder_Result       = os.path.join('./testdata','Final_output_folder','testout3')
    shutil.rmtree(Temporary_Result_Folder_Result,ignore_errors=True)
    
    ###The pathes for all inputs 
    Path_DEM_big           = os.path.join(Data_Folder, 'DEM_big_merit.tif')
    Path_DEM_small         = os.path.join(Data_Folder, 'DEM_samll_merit.tif')
    Path_Lake_ply          = os.path.join(Data_Folder, 'HyLake.shp')
    Path_bkf_wd            = os.path.join(Data_Folder, 'Bkfullwidth_depth.shp')
    Path_Landuse           = os.path.join(Data_Folder, 'landuse.tif')
    Path_Roughness_landuse = os.path.join(Data_Folder, 'Landuse.csv')

    ###Generate test resuts, for option 1 
    RTtool=LRRT(dem_in = Path_DEM_small,WidDep = Path_bkf_wd,
               Lakefile = Path_Lake_ply,Landuse = Path_Landuse,
               Landuseinfo = Path_Roughness_landuse,
               OutputFolder = Final_Result_Folder_Result,
               TempOutFolder = Temporary_Result_Folder_Result,
               )
    RTtool.Generatmaskregion()
    
    """Evaluate raster MASK 
    MASK is a mask raster stored in grass database, which indicate 
    the spatial processing extent.     
    """ 
    ### transfer expected raster alllake into np array Expected_Mask_Array   
    Expected_Mask_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'MASK')
    ### transfer resulted raster alllake into np array Result_Mask_Array
    Result_Mask_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'MASK')
                                                 
    ### compare two Expected_Mask_Array and Result_Mask_Array 
    assert (Expected_Mask_Array == Result_Mask_Array).all()
    shutil.rmtree(Temporary_Result_Folder_Result) 

    ### test another option to define processing extent
    ### Folder where the output will be generated 
    Final_Result_Folder_Expected     = os.path.join('./testdata','Final_output_folder','Expected_Out_Cord')
    Temporary_Result_Folder_Expected = os.path.join('./testdata','Temporary_output_folder','Expected_Out_Cord')
    ###Folder where the output will be generated 
    Temporary_Result_Folder_Result   = os.path.join('./testdata','Temporary_output_folder','testout4')
    Final_Result_Folder_Result       = os.path.join('./testdata','Final_output_folder','testout4')
    shutil.rmtree(Temporary_Result_Folder_Result,ignore_errors=True)
    
    RTtool=LRRT(dem_in = Path_DEM_small,WidDep = Path_bkf_wd,
               Lakefile = Path_Lake_ply,Landuse = Path_Landuse,
               Landuseinfo = Path_Roughness_landuse,
               OutputFolder = Final_Result_Folder_Result,
               TempOutFolder = Temporary_Result_Folder_Result,
               )
    OutletPointxy = [-92.387,49.09]
    RTtool.Generatmaskregion(OutletPoint = OutletPointxy)    
    """Evaluate raster MASK 
    MASK is a mask raster stored in grass database, which indicate 
    the spatial processing extent.     
    """ 
    ### transfer expected raster alllake into np array Expected_Mask_Array   
    Expected_Mask_Array = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Expected,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'MASK')
    ### transfer resulted raster alllake into np array Result_Mask_Array
    Result_Mask_Array   = Return_Raster_As_Array(grassdb = os.path.join(Temporary_Result_Folder_Result,'grassdata_toolbox'),
                                                 grass_location = 'Geographic',
                                                 raster_mn = 'MASK')
                                                 
    ### compare two Expected_Mask_Array and Result_Mask_Array 
    assert (Expected_Mask_Array == Result_Mask_Array).all()
    
    shutil.rmtree(Temporary_Result_Folder_Result) 

    
                   
    
        