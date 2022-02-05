import geopandas 
import numpy as np
import os
import pandas as pd 


def save_modified_attributes_to_outputs(mapoldnew_info,tempfolder,OutputFolder,cat_name,riv_name,Path_final_riv,dis_col_name='SubId'):


    print(tempfolder)
    NEED_TO_REMOVE_IDS = ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId"]

    
    
    if riv_name != '#':

        riv_pd = geopandas.read_file(Path_final_riv)
        riv_pd['Old_SubId'] = riv_pd['SubId']
        
        cat_pd = mapoldnew_info.drop(columns = 'geometry').copy(deep=True)
        # remove all columns 
        riv_pd = riv_pd[['geometry','Old_SubId']]        
        riv_pd = pd.merge(riv_pd, cat_pd, on='Old_SubId', how='left')
        
        riv_pd_nncls_routing_info = riv_pd[riv_pd['Lake_Cat'] != 2][['SubId','DowSubId']].copy(deep=True)
        remove_channel = []
        for subid in riv_pd_nncls_routing_info['SubId'].values:
            if subid not in riv_pd_nncls_routing_info['DowSubId'].values:
                remove_channel.append(subid)                
        riv_pd = riv_pd[~riv_pd['SubId'].isin(remove_channel)]                
        riv_pd = riv_pd.dissolve(by=dis_col_name, aggfunc='first')
        
        
        mapoldnew_info = mapoldnew_info.dissolve(by=dis_col_name, aggfunc='first')
        mapoldnew_info["centroid_y"] = mapoldnew_info.geometry.centroid.y
        mapoldnew_info["centroid_x"] = mapoldnew_info.geometry.centroid.x
        
        cat_c_x_y = mapoldnew_info[["centroid_y","centroid_x"]].copy(deep=True)
        riv_pd = riv_pd.drop(columns = ["centroid_y","centroid_x"])
        riv_pd = riv_pd.join(cat_c_x_y) 
        
        cat_colnms = riv_pd.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(NEED_TO_REMOVE_IDS)]
        riv_pd = riv_pd.drop(columns=drop_cat_colnms)
        riv_pd.to_file(os.path.join(OutputFolder,riv_name))
 
        cat_colnms = mapoldnew_info.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(NEED_TO_REMOVE_IDS)]
        mapoldnew_info = mapoldnew_info.drop(columns=drop_cat_colnms)
        mapoldnew_info.to_file(os.path.join(OutputFolder,cat_name))
  
    
    else: 

        mapoldnew_info = mapoldnew_info.dissolve(by=dis_col_name, aggfunc='first')
        mapoldnew_info["centroid_y"] = mapoldnew_info.geometry.centroid.y
        mapoldnew_info["centroid_x"] = mapoldnew_info.geometry.centroid.x
        
        cat_colnms = mapoldnew_info.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(["SHAPE","SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId"])]
        mapoldnew_info = mapoldnew_info.drop(columns=drop_cat_colnms)
        mapoldnew_info = mapoldnew_info.drop(columns=drop_column)
        mapoldnew_info.to_file(os.path.join(OutputFolder,cat_name))


def clean_attribute_name_arcgis(table,names):
    remove_column_names = table.columns[np.logical_not(np.isin(table.columns,names))]
    table = table.drop(columns=remove_column_names)
    return table 