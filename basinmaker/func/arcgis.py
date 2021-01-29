from arcgis import GIS
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import os
import pandas as pd 

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
    mask2 = finalcat_ply['IsLake'] != 2
    mask = np.logical_and(mask1,mask2)
    
    finalcat_ply.loc[mask,'HyLakeId'] = -9999
    finalcat_ply.loc[mask,'LakeVol'] = -9999
    finalcat_ply.loc[mask,'LakeArea'] = -9999
    finalcat_ply.loc[mask,'LakeDepth'] = -9999
    finalcat_ply.loc[mask,'Laketype'] = -9999
    finalcat_ply.loc[mask,'IsLake'] = -9999
    
    return finalcat_ply


def save_modified_attributes_to_outputs(mapoldnew_info,tempfolder,OutputFolder,cat_name,riv_name,Path_final_riv,dis_col_name='SubId'):

    mapoldnew_info.spatial.to_featureclass(location=os.path.join(tempfolder,'updateattri.shp'))
    arcpy.Dissolve_management(os.path.join(tempfolder,'updateattri.shp'), os.path.join(OutputFolder,cat_name), [dis_col_name])
    arcpy.JoinField_management(os.path.join(OutputFolder,cat_name), dis_col_name, os.path.join(tempfolder,'updateattri.shp'), dis_col_name)
    arcpy.DeleteField_management(os.path.join(OutputFolder,cat_name), 
        ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters"]
    )
    
    
    if riv_name != '#':
        arcpy.CalculateGeometryAttributes_management(os.path.join(OutputFolder, cat_name), [["centroid_x", "CENTROID_X"], ["centroid_y", "CENTROID_Y"]])
        riv_pd = pd.DataFrame.spatial.from_featureclass(Path_final_riv)
        riv_pd['SubID_Oldriv'] = riv_pd['SubId']
        # remove all columns 
        riv_pd = riv_pd[['SHAPE','SubID_Oldriv']]
        riv_pd.spatial.to_featureclass(location=os.path.join(tempfolder,'riv_no_attri.shp'))
        arcpy.SpatialJoin_analysis(os.path.join(tempfolder,'riv_no_attri.shp'), os.path.join(OutputFolder, cat_name), os.path.join(tempfolder,'riv_attri.shp'),match_option='WITHIN')
        arcpy.Dissolve_management(os.path.join(tempfolder,'riv_attri.shp'), os.path.join(OutputFolder,riv_name), ["SubId"])
        arcpy.JoinField_management(os.path.join(OutputFolder,riv_name), "SubId", os.path.join(tempfolder,'riv_attri.shp'), "SubId")
        arcpy.DeleteField_management(os.path.join(OutputFolder,cat_name), 
            ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters"]
        )

def clean_attribute_name_arcgis(table,names):
    remove_column_names = table.columns[np.logical_not(np.isin(table.columns,names))]
    table = table.drop(columns=remove_column_names)
    return table 
    
    
