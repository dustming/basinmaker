import numpy as np
import sys
import os
import csv
import tempfile 
import copy 
import pandas as pd
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from basinmaker.func.pdtable import *
from osgeo import gdal, ogr

def decode_hru_attri_ids(HRU_temp1,lakehruinfo,Landuse_ID,Soil_ID,Veg_ID,Other_Ply_ID_1,Other_Ply_ID_2):
    #     100000000000000, 
    #          0007422000,
    #          2147483647,
#    'sub':       1000000,
#    'landuse':     10000,
#    'soil':          100,
#    'o1':              1,

    lakehruinfo_id = lakehruinfo[['HRULake_ID','SubId','HRU_IsLake']]
    # remove polygon outside of the subbasin 
    HRU_temp1 = HRU_temp1[HRU_temp1['HRU_ID_T'] != -9999]
    HRU_temp1['HRU_ID_str'] = HRU_temp1['HRU_ID_T'].astype('int64').astype(str)
#    HRU_temp1['HRU_ID_str'] = HRU_temp1['HRU_ID_T'].str.zfill(10)
           
    HRU_temp1['HRULake_ID'] = HRU_temp1['HRU_ID_str'].str[:-6].astype(int)
    HRU_temp1 = pd.merge(HRU_temp1, lakehruinfo_id, how='inner', on = 'HRULake_ID')
    HRU_temp1[Landuse_ID] = HRU_temp1['HRU_ID_str'].str[-6:-4].astype(int)
    HRU_temp1[Veg_ID] = HRU_temp1['HRU_ID_str'].str[-6:-4].astype(int)
    HRU_temp1[Soil_ID] = HRU_temp1['HRU_ID_str'].str[-4:-2].astype(int)
    HRU_temp1[Other_Ply_ID_1] = HRU_temp1['HRU_ID_str'].str[-2:].astype(int)
    HRU_temp1[Other_Ply_ID_2] = HRU_temp1['HRU_ID_str'].str[-2:].astype(int)
    HRU_temp1 = HRU_temp1.loc[(HRU_temp1['HRULake_ID'] > 0) | (HRU_temp1['SubId'] > 0)]
    
    return HRU_temp1
       
def raster_to_vector(output,input,raster_par_list):
    x_min = raster_par_list[0]
    x_max = raster_par_list[1]
    y_min = raster_par_list[2]
    y_max = raster_par_list[3]
    x_res = raster_par_list[4]
    pixel_size = raster_par_list[5]
    y_res = raster_par_list[6]
    crs = raster_par_list[7]
    
    src_ds = gdal.Open(input)
    srcband = src_ds.GetRasterBand(1)


    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(output)
    dst_layer = dst_ds.CreateLayer(output, srs = crs )

    fd = ogr.FieldDefn("HRU_ID_T", ogr.OFTInteger)
    dst_layer.CreateField(fd)
    dst_field = dst_layer.GetLayerDefn().GetFieldIndex("HRU_ID_T")

    gdal.Polygonize( srcband, None, dst_layer, dst_field, [], callback=None)
        
def vector_to_raster(output,input,attribute,raster_par_list,touch = "False",type = gdal.GDT_Float32):
    
    x_min = raster_par_list[0]
    x_max = raster_par_list[1]
    y_min = raster_par_list[2]
    y_max = raster_par_list[3]
    x_res = raster_par_list[4]
    pixel_size = raster_par_list[5]
    y_res = raster_par_list[6]
    crs = raster_par_list[7]
    
    source_ds = ogr.Open(input)
    source_layer = source_ds.GetLayer()

    target_ds = gdal.GetDriverByName('GTiff').Create(output, x_res, y_res, 1, type)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(0)

    target_ds.SetProjection(crs.ExportToWkt())

    gdal.RasterizeLayer(target_ds, 
                        [1], source_layer, 
                        options = ["ALL_TOUCHED=%s"%touch, "ATTRIBUTE=%s" % attribute])
    
    return 

def array_to_raster(output,input,raster_par_list,touch = "False",type = gdal.GDT_Int32):
    
    x_min = raster_par_list[0]
    x_max = raster_par_list[1]
    y_min = raster_par_list[2]
    y_max = raster_par_list[3]
    x_res = raster_par_list[4]
    pixel_size = raster_par_list[5]
    y_res = raster_par_list[6]
    crs = raster_par_list[7]    

    target_ds = gdal.GetDriverByName('GTiff').Create(output, x_res, y_res, 1, type)
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    band = target_ds.GetRasterBand(1)
    band.SetNoDataValue(-9999)

    target_ds.SetProjection(crs.ExportToWkt())
    target_ds.GetRasterBand(1).WriteArray(input)
    target_ds.FlushCache()
    return 

def gdal_slope_raster(input,output):
    gdal.DEMProcessing(output,input, 'slope')

def gdal_aspect_raster(input,output):
    
    #This command outputs a 32-bit float raster with values 
    #between 0° and 360° representing the azimuth that slopes are facing. 
    #The definition of the azimuth is such that : 
    #0° means that the slope is facing the North, 
    #90° it’s facing the East, 180° it’s facing the South and 
    #270° it’s facing the West (provided that the top of your input raster is north oriented). 
    #The aspect value -9999 is used as the nodata value to indicate undefined 
    #aspect in flat areas with slope=0.
    
    gdal.DEMProcessing(output,input, 'aspect')
                
def proj_clip_raster(input,output,dst_src,cutlineDSName= '#'):
    
    input_raster = gdal.Open(input)
    output_raster = output
    if cutlineDSName == '#':
        warp = gdal.Warp(output_raster,input_raster,dstSRS=dst_src)    
    else:
        warp = gdal.Warp(output_raster, 
                    input_raster, 
                    cutlineDSName=cutlineDSName,
                    cropToCutline=True,
                    dstNodata = -9999)



def RasterHRUUnionInt32(OutputFolder,tempfolder,Merge_layer_shp_list,
                        Merge_ID_list,Sub_Lake_HRU_Layer,Sub_ID,
                        Landuse_ID,Soil_ID,Veg_ID,Other_Ply_ID_1,
                        Other_Ply_ID_2,pixel_size,lakehruinfo):
                        
    path_to_sub_lake_raster = os.path.join(tempfolder,'sub_lake.tif')
    path_to_landuse_raster = os.path.join(tempfolder,'landuse.tif')
    path_to_soil_raster = os.path.join(tempfolder,'soil.tif')
    path_to_veg_raster = os.path.join(tempfolder,'veg.tif')
    path_to_other1_raster = os.path.join(tempfolder,'o1.tif')
    path_to_other2_raster = os.path.join(tempfolder,'o2.tif')
    path_to_hru_temp_raster = os.path.join(tempfolder,'HRU.tif')
    path_to_lake_mask_raster = os.path.join(tempfolder,'lake_hru.tif')
    path_to_land_mask_raster = os.path.join(tempfolder,'land_hru.tif')
    path_to_landhru_shp= os.path.join(tempfolder,'land_hru.shp')
    path_to_lakehru_shp= os.path.join(tempfolder,'lake_hru.shp')
    path_to_landuse_shp = os.path.join(tempfolder,'landuse.shp')
    path_to_soil_shp = os.path.join(tempfolder,'soil.shp')
    path_to_veg_shp = os.path.join(tempfolder,'veg.shp')
    path_to_other1_shp = os.path.join(tempfolder,'o1.shp')
    path_to_other2_shp = os.path.join(tempfolder,'o2.shp')
    path_to_hru_temp_shp = os.path.join(tempfolder,'HRU.shp')

    raster_value_multi = {
    #                   2147483647,
    Sub_ID:                1000000,
    Landuse_ID:              10000,
    Soil_ID:                   100,
    Other_Ply_ID_1:              1,
    Veg_ID:                      1,
    Other_Ply_ID_2:              1,
    }
 
     #          12345678910
     #          0007422000,
     #          2147483647,
 #    'sub':       1000000,
 #    'landuse':     10000,
 #    'soil':          100,
 #    'o1':              1,
 
 
    # determine raster parameters        
    vector_fn = Sub_Lake_HRU_Layer
    source_ds = ogr.Open(vector_fn)
    source_layer = source_ds.GetLayer()
    x_min, x_max, y_min, y_max = source_layer.GetExtent()
    ras_crs = source_layer.GetSpatialRef()
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)
    raster_par_list = [x_min, x_max, y_min, y_max,x_res,pixel_size,y_res,ras_crs]

    # rasterize 
    vector_to_raster(path_to_sub_lake_raster,Sub_Lake_HRU_Layer,'HRULake_ID',raster_par_list,touch='False')
               
    
    dissolve_filedname_list = ["HRULake_ID"]
    
    sub_raster = gdal.Open(path_to_sub_lake_raster)
    sub_raster_a = np.array(sub_raster.GetRasterBand(1).ReadAsArray())
    
    #check if some subbasin is missing 
    raster_unique_id = np.unique(sub_raster_a)
    raster_unique_id = raster_unique_id[raster_unique_id > 0]
    vector_unique_id = np.unique(lakehruinfo['HRULake_ID'].values)
    vector_unique_id = vector_unique_id[vector_unique_id > 0]
    mask = np.isin(vector_unique_id,raster_unique_id)
    
    if len(vector_unique_id[~mask]) > 0:
        print("Please increase the pixel_size parameter")
        print("The following element in path_output_folder/finalcat_hru_lake_info.shp is lost ")
        print(vector_unique_id[~mask])
        sys.exit()
        
    
    HRU_shape = sub_raster_a.shape
    HRU_a = sub_raster_a * raster_value_multi[Sub_ID]
    for i in range(0,len(Merge_layer_shp_list)):
        merge_shpfile = Merge_layer_shp_list[i]
        ID_shpfile = Merge_ID_list[i]
        out_raster = os.path.join(tempfolder,"%s_ras_layer.tif"%(ID_shpfile))
        vector_to_raster(out_raster,merge_shpfile,ID_shpfile,raster_par_list,touch='False')

        raster_layer_i = gdal.Open(out_raster)
        raster_layer_i_a = np.array(raster_layer_i.GetRasterBand(1).ReadAsArray())
        if raster_layer_i_a.shape[0] == HRU_shape[0] and raster_layer_i_a.shape[1] == HRU_shape[1]:
            HRU_a = HRU_a + raster_layer_i_a * raster_value_multi[ID_shpfile]
        else:
            print(ID_shpfile,"polygon may not fully cover the subbasin domain ",land_raster_a.shape)
            sys.exit()
        

    HRU_a[HRU_a < raster_value_multi[Sub_ID]] = -9999
    array_to_raster(path_to_hru_temp_raster,HRU_a,raster_par_list,touch = "False",type = gdal.GDT_Int32)
    
    raster_to_vector(path_to_hru_temp_shp,path_to_hru_temp_raster,raster_par_list)

    print("union done")

    return path_to_hru_temp_shp