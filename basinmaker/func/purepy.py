import geopandas 
import numpy as np
import os
import pandas as pd 


def save_modified_attributes_to_outputs(mapoldnew_info,tempfolder,OutputFolder,cat_name,riv_name,Path_final_riv,dis_col_name='SubId'):


    NEED_TO_REMOVE_IDS = ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId","SubIdt2"]

    
    
    if riv_name != '#':

        riv_pd = geopandas.read_file(Path_final_riv)
        riv_pd['Old_SubId'] = riv_pd['SubId']
        
        cat_pd = mapoldnew_info.drop(columns = 'geometry').copy(deep=True)
        # remove all columns 
        riv_pd = riv_pd[['geometry','Old_SubId']]        
        riv_pd = pd.merge(riv_pd, cat_pd, on='Old_SubId', how='left')             
        riv_pd = riv_pd.dissolve(by=dis_col_name, aggfunc='first',as_index=False)
        
        
        mapoldnew_info = mapoldnew_info.dissolve(by=dis_col_name, aggfunc='first',as_index=False)
        mapoldnew_info = add_centroid_in_wgs84(mapoldnew_info,"centroid_x","centroid_y")
        
        cat_c_x_y = mapoldnew_info[["centroid_y","centroid_x"]].copy(deep=True)
        riv_pd = riv_pd.drop(columns = ["centroid_y","centroid_x"])
        riv_pd = riv_pd.join(cat_c_x_y) 

        riv_pd_nncls_routing_info = mapoldnew_info[mapoldnew_info['Lake_Cat'] != 2][['SubId','DowSubId']].copy(deep=True)
        remove_channel = []
        for subid in riv_pd_nncls_routing_info['SubId'].values:
            if subid not in riv_pd_nncls_routing_info['DowSubId'].values:
                remove_channel.append(subid)                
        riv_pd = riv_pd[~riv_pd.SubId.isin(remove_channel)]   
        cat_colnms = riv_pd.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(NEED_TO_REMOVE_IDS)]
        riv_pd = riv_pd.drop(columns=drop_cat_colnms)
        if len(riv_pd) > 0:
            riv_pd.to_file(os.path.join(OutputFolder,riv_name))
        
        mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'RivSlope'] = -1.2345
        mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'RivLength'] = -1.2345
        mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'FloodP_n'] = -1.2345
        mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'Ch_n'] = -1.2345
        mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'Max_DEM'] = -1.2345
        mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'Min_DEM'] = -1.2345
        
        cat_colnms = mapoldnew_info.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(NEED_TO_REMOVE_IDS)]
        mapoldnew_info = mapoldnew_info.drop(columns=drop_cat_colnms)
        mapoldnew_info.to_file(os.path.join(OutputFolder,cat_name))
        
        create_geo_jason_file(os.path.join(OutputFolder,cat_name))
  
    else: 

        mapoldnew_info = mapoldnew_info.dissolve(by=dis_col_name, aggfunc='first',as_index=False)
    
        if "centroid_y" in mapoldnew_info.columns:

            mapoldnew_info = add_centroid_in_wgs84(mapoldnew_info,"centroid_x","centroid_y")
            mapoldnew_info["SubId"] = mapoldnew_info.index
            riv_pd_nncls_routing_info = mapoldnew_info[mapoldnew_info['Lake_Cat'] != 2][['SubId','DowSubId']].copy(deep=True)
            remove_channel = []
            for subid in riv_pd_nncls_routing_info['SubId'].values:
                if subid not in riv_pd_nncls_routing_info['DowSubId'].values:
                    remove_channel.append(subid)                
                            
            mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'RivSlope'] = -1.2345
            mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'RivLength'] = -1.2345
            mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'FloodP_n'] = -1.2345
            mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'Ch_n'] = -1.2345
            mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'Max_DEM'] = -1.2345
            mapoldnew_info.loc[mapoldnew_info.SubId.isin(remove_channel),'Min_DEM'] = -1.2345
            
                
        cat_colnms = mapoldnew_info.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(["SHAPE","SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId"])]
        mapoldnew_info = mapoldnew_info.drop(columns=drop_cat_colnms)
        mapoldnew_info.to_file(os.path.join(OutputFolder,cat_name))
        return mapoldnew_info


def Remove_Unselected_Lake_Attribute_In_Finalcatinfo_purepy(finalcat_ply, Conn_Lake_Ids):
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
    
def clean_attribute_name_purepy(table,names):
    remove_column_names = table.columns[np.logical_not(np.isin(table.columns,names))]
    table = table.drop(columns=remove_column_names)
    return table 
    
def clean_geometry_purepy(data):
    narow = ~data['geometry'].isna()
    emrow = ~data.is_empty
    arearow = data.area > 0.00000001

    row1 = np.logical_and(narow,emrow)
    rowselect = np.logical_and(arearow,row1)
    data = data.loc[rowselect]

    return data   
    
def add_area_in_m2(data,prj_crs,area_col):
    src_src = data.crs
    tost = data.copy()

    tost= data.to_crs(prj_crs)
    tost[area_col] = tost.area
    
    out= tost.copy(deep=True).to_crs(src_src)

    return out 

def add_centroid_in_wgs84(data,colx,coly):
    src_src = data.crs
    tost = data.copy()
    
    tost= tost.to_crs('EPSG:4326')

    tost[coly] = tost.geometry.centroid.y
    tost[colx] = tost.geometry.centroid.x
    
    out= tost.copy(deep=True).to_crs(src_src)
    
    return out
    

def create_geo_jason_file(Input_Polygon_path):


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
        input_pd = geopandas.read_file(input_path)
        
        input_wgs_84 = input_pd.to_crs('EPSG:4326') 

        
        
        if 'finalcat_info' in Input_file_name[i] or "finalcat_info_riv" in Input_file_name[i]:
            input_wgs_84['rvhName'] = input_wgs_84['SubId'].astype(int).astype(str)
            
        
        input_tojson = input_wgs_84
        
        for TOLERANCE in TOLERANCEs:                               
            input_tojson['geometry'] = input_tojson.simplify(TOLERANCE)
            input_tojson.to_file(output_jason_path, driver="GeoJSON") 

            json_file_size = os.stat(output_jason_path).st_size/1024/1024 #to MB
            if json_file_size <= 100:
                break
                                
    return     
    
    
    