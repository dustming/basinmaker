from arcgis import GIS
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import os
import pandas as pd 
import tempfile
import arcpy.cartography as CA
from json import load, JSONEncoder
import json
import shutil
#####
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

def select_feature_by_attributes_arcgis(input,Attri_NM,Attri_v,output):
    where_clause = '"%s" IN' % (Attri_NM) 
    where_clause = where_clause + " (" 
    for i in range(0,len(Attri_v)):
        if i == 0:
            where_clause = where_clause + str(Attri_v[i])
        else:
            where_clause = where_clause + "," + str(Attri_v[i])
    where_clause = where_clause + ")"
    
    arcpy.Select_analysis(input, output, where_clause)
    return
##################


def Remove_Unselected_Lake_Attribute_In_Finalcatinfo_Arcgis(finalcat_ply, Conn_Lake_Ids):
    """Functions will set lake id not in Conn_Lake_Ids to -1.2345 in attribute
        table of Path_Finalcatinfo
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """
    
    mask1 = np.logical_not(finalcat_ply['HyLakeId'].isin(Conn_Lake_Ids))
    mask2 = finalcat_ply['Lake_Cat'] != 2
    mask = np.logical_and(mask1,mask2)
    
    finalcat_ply.loc[mask,'HyLakeId'] = 0
    finalcat_ply.loc[mask,'LakeVol'] = 0
    finalcat_ply.loc[mask,'LakeArea'] = 0
    finalcat_ply.loc[mask,'LakeDepth'] = 0
    finalcat_ply.loc[mask,'Laketype'] =0
    finalcat_ply.loc[mask,'Lake_Cat'] = 0
    
    return finalcat_ply


def save_modified_attributes_to_outputs(mapoldnew_info,tempfolder,OutputFolder,cat_name,riv_name,Path_final_riv,dis_col_name='SubId'):

    mapoldnew_info.spatial.to_featureclass(location=os.path.join(tempfolder,'updateattri.shp'),overwrite=True,sanitize_columns=False)
    arcpy.Dissolve_management(os.path.join(tempfolder,'updateattri.shp'), os.path.join(OutputFolder,cat_name), [dis_col_name])
    arcpy.JoinField_management(os.path.join(OutputFolder,cat_name), dis_col_name, os.path.join(tempfolder,'updateattri.shp'), dis_col_name)
    arcpy.DeleteField_management(os.path.join(OutputFolder,cat_name), 
        ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters"]
    )
    
    
    if riv_name != '#':
        arcpy.CalculateGeometryAttributes_management(os.path.join(OutputFolder, cat_name), [["centroid_x", "CENTROID_X"], ["centroid_y", "CENTROID_Y"]],coordinate_system = arcpy.SpatialReference(4326))
        cat_colnms = mapoldnew_info.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(["SHAPE","SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId"])]
        cat_pd = mapoldnew_info.drop(columns=drop_cat_colnms)

        riv_pd = pd.DataFrame.spatial.from_featureclass(Path_final_riv)
        riv_pd['Old_SubId'] = riv_pd['SubId']
        # remove all columns 
        riv_pd = riv_pd[['SHAPE','Old_SubId']]        
        riv_pd = pd.merge(riv_pd, cat_pd, on='Old_SubId', how='left')
        
        riv_pd = riv_pd.drop(columns=['Old_SubId'])
        mask = np.logical_or(riv_pd['RivLength'] > 0,riv_pd['Lake_Cat'] > 0)
        riv_pd = riv_pd[mask]
        riv_pd.drop(columns=['centroid_x','centroid_y'])
        riv_pd.spatial.to_featureclass(location=os.path.join(tempfolder,'riv_attri.shp'),overwrite=True,sanitize_columns=False)
        arcpy.Dissolve_management(os.path.join(tempfolder,'riv_attri.shp'), os.path.join(OutputFolder,riv_name), ["SubId"])
        arcpy.JoinField_management(os.path.join(OutputFolder,riv_name), "SubId", os.path.join(tempfolder,'riv_attri.shp'), "SubId")
        arcpy.DeleteField_management(os.path.join(OutputFolder,riv_name), 
            ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters",'centroid_x','centroid_y']
        )
    
    # if "finalcat_info" in cat_name:
    #     create_geo_jason_file(os.path.join(OutputFolder,cat_name))
    
def clean_attribute_name_arcgis(table,names):
    remove_column_names = table.columns[np.logical_not(np.isin(table.columns,names))]
    table = table.drop(columns=remove_column_names)
    return table 


####
def create_geo_jason_file(Input_Polygon_path):

    arcpy.env.overwriteOutput = True
    product_dir = os.path.dirname(Input_Polygon_path)
    Names_in = os.path.basename(Input_Polygon_path).split('_')
    n_charc = len(Names_in)
    version  = Names_in[n_charc - 1][0:4]
    TOLERANCEs = [0.0001,0.0005,0.001,0.005,0.01,0.05]


    head_name_cat = "finalcat_info"
    head_name_riv = "finalcat_info_riv"
    head_name_slake = "sl_connected_lake"
    head_name_nlake = "sl_non_connected_lake"

    Input_file_name = []
    Output_file_name = []                            
    if 'v' in version:
        Input_file_name = [
                           head_name_cat + "_"+version+'.shp',
                           head_name_riv + "_"+version+'.shp',
                           head_name_slake + "_"+version+'.shp',
                           head_name_nlake + "_"+version+'.shp',
                           ]
        Output_file_name = [
                           head_name_cat + "_"+version+'.geojson',
                           head_name_riv + "_"+version+'.geojson',
                           head_name_slake + "_"+version+'.geojson',
                           head_name_nlake + "_"+version+'.geojson',
                           ]
    else:
        Input_file_name = [
                           head_name_cat +'.shp',
                           head_name_riv +'.shp',
                           head_name_slake +'.shp',
                           head_name_nlake +'.shp',
                           ]
        Output_file_name = [
                           head_name_cat +'.geojson',
                           head_name_riv +'.geojson',
                           head_name_slake +'.geojson',
                           head_name_nlake +'.geojson',
                           ]
    created_jason_files = []   
    created_jason_files_lake_riv = [] 
    
    for i  in  range(0,len(Input_file_name)):
        input_path = os.path.join(product_dir,Input_file_name[i])
        output_jason_path = os.path.join(product_dir,Output_file_name[i])
        if not os.path.exists(input_path):
            continue 
        created_jason_files.append(output_jason_path) 
        
        if 'finalcat_info_riv' in Input_file_name[i] or 'connected_lake' in Input_file_name[i]:
            created_jason_files_lake_riv.append(output_jason_path)  
                            
        # reproject to WGS84
        input_wgs_84 = os.path.join(tempfile.gettempdir(),"input_wgs_84.shp")
        arcpy.Project_management(input_path, input_wgs_84, arcpy.SpatialReference(int(4326)))


        if 'finalcat_info' in Input_file_name[i] or "finalcat_info_riv" in Input_file_name[i]:
            arcpy.AddField_management(input_wgs_84, 'rvhName', "TEXT")
            arcpy.CalculateField_management(input_wgs_84, 'rvhName', "'sub' + str(int(\"!SubId!\"))", "PYTHON3")

        arcpy.RepairGeometry_management(input_wgs_84)

        for TOLERANCE in TOLERANCEs:
            input_wgs_84_simplify = os.path.join(tempfile.gettempdir(),"input_wgs_84_simplify.shp")
            if arcpy.Exists(input_wgs_84_simplify):
                arcpy.Delete_management(input_wgs_84_simplify)
                
                
            if "finalcat_info_riv" not in Input_file_name[i]:
                CA.SimplifyPolygon(input_wgs_84, input_wgs_84_simplify, "POINT_REMOVE", TOLERANCE)
            else:
                CA.SimplifyLine(input_wgs_84, input_wgs_84_simplify, "POINT_REMOVE", TOLERANCE)

            arcpy.FeaturesToJSON_conversion(input_wgs_84_simplify, 
                                            output_jason_path,"FORMATTED","NO_Z_VALUES","NO_M_VALUES","GEOJSON")
                                
            json_file_size = os.stat(output_jason_path).st_size/1024/1024 #to MB
            if json_file_size <= 100:
                break                                

    if len(created_jason_files_lake_riv) > 1 and os.stat(os.path.join(product_dir,Output_file_name[1])).st_size/1024/1024 < 500:
        for i in range(0,len(created_jason_files_lake_riv)):
            injson2 = load(open(created_jason_files_lake_riv[i]))
            if 'finalcat_info_riv' in created_jason_files_lake_riv[i]:
                new_features = []
                for element in injson2["features"]:
                    if element["properties"]["Lake_Cat"] == 0:    
                        new_features.append(element) 
                injson2["features"] = new_features

            if i == 0:
                output_jason_lake_riv = injson2
            else:
                output_jason_lake_riv['features'] += injson2['features']
                
        with open(os.path.join(product_dir,'routing_product_lake_river.geojson'), 'w', encoding='utf-8') as f:
            json.dump(output_jason_lake_riv, f, ensure_ascii=False, indent=4)
    else:
        shutil.copy(created_jason_files_lake_riv[0], os.path.join(product_dir,'routing_product_lake_river.geojson'))             


                                 
    return 


    
    
