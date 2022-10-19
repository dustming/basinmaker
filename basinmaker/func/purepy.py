import geopandas
import numpy as np
import os
import pandas as pd
from json import load, JSONEncoder
import json
import requests
import sys
#import shapely.wkt
#import pygeos as pg
from osgeo import gdal, ogr
from rasterstats import zonal_stats

def zonal_stats_pd(shp_gpd,raster,stats,key):

    result = zonal_stats(shp_gpd, raster, stats = stats,geojson_out=True,all_touched=True)

    reault_list_dic = list(map(lambda x : {key:x['properties'][key],'mean':x['properties']['mean']}, result))

    reault_pd = pd.DataFrame(reault_list_dic)

    return reault_pd



def ZonalStats(shp_gpd, raster, stats,key,col):
    # shape - shapefile path
    # raster - raster path
    # stats - stats as list, f.e. 'min mean max' ; 'min'
    # the result is final_gdf as GeoDataFrame

    nodata_pd = shp_gpd.copy(deep=True)
    i = 0
    while len(nodata_pd) >0 and i <=5:
        temp_pd = zonal_stats_pd(nodata_pd,raster,stats,key)
        bad_hru_values = temp_pd[~(temp_pd['mean'] > 0)]['HRU_ID_New']
        nodata_pd = nodata_pd[nodata_pd['HRU_ID_New'].isin(bad_hru_values)].copy(deep=True)
        if i == 0:
            reault_pd = temp_pd[temp_pd['mean'] > 0]
            i = i + 1
        else:
            i = i + 1
            good_pd = temp_pd[temp_pd['mean'] > 0]
            reault_pd = reault_pd.append(good_pd)

    if len(nodata_pd) > 0:
        nodata_pd['mean'] = nodata_pd[col]
        reault_pd = reault_pd.append(nodata_pd)
#    print(len(shp_gpd),len(nodata_pd),len(reault_pd))
    return reault_pd


def Reproj_Clip_Dissolve_Simplify_Polygon_purepy(
    layer_path, Class_Col, tempfolder,mask_layer,Class_NM_Col = '#',info_table = '#'
):
    """Preprocess user provided polygons

    Function that will reproject clip input polygon with subbasin polygon
    and will dissolve the input polygon based on their ID, such as landuse id
    or soil id.

    Parameters
    ----------
    processing                        : qgis object
    context                           : qgis object
    layer_path                        : string
        The path to a specific polygon, for example path to landuse layer
    Project_crs                       : string
        the EPSG code of a projected coodinate system that will be used to
        calcuate HRU area and slope.
    trg_crs                           : string
        the EPSG code of a  coodinate system that will be used to
        calcuate reproject input polygon
    Class_Col                         : string
        the column name in the input polygon (layer_path) that contains
        their ID, for example land use ID or soil ID.
    Layer_clip                        : qgis object
        A shpfile with extent of the watershed, will be used to clip input
        input polygon
    Notes
    -------
        # TODO: May be add some function to simplify the input polygons
                for example, remove the landuse type with small areas
                or merge small landuse polygon into the surrounding polygon

    Returns:
    -------
        layer_dis                  : qgis object
            it is a polygon after preprocess
    """

    new_crs = mask_layer.crs
    data = geopandas.read_file(layer_path)

    projected = data.to_crs(new_crs)
    projected = projected.explode(ignore_index=True)
    clipped = projected.clip(mask_layer)
#    dissolved = clipped.dissolve(by=[Class_Col], aggfunc='first',as_index=False)
#    print(Class_Col,new_crs,clipped.crs)
    cleaned = clean_geometry_purepy(clipped,1)
    if Class_NM_Col != '#':
        if Class_Col not in clipped.columns:
            print(Class_Col," not in the attribute table of provided shapefile")
            sys.exit()
        if Class_Col not in info_table.columns:
            print(Class_Col," not in the attribute table of provided info table")
            sys.exit()
        if Class_NM_Col not in info_table.columns:
            print(Class_NM_Col," not in the attribute table of provided info table")
            sys.exit()

        info_table_copy = info_table.copy(deep=True)
        info_table_copy = info_table_copy.drop_duplicates(subset=[Class_NM_Col], keep='last',ignore_index=True)
        info_table_copy['New_ID'] = info_table_copy[Class_Col]
        info_table_copy = info_table_copy[[Class_NM_Col,'New_ID']]
        info_table_copy2 = pd.merge(info_table, info_table_copy, how='inner', on = Class_NM_Col).copy(deep=True)

        cleaned = pd.merge(cleaned, info_table_copy2, how='inner', on = Class_Col)
        cleaned[Class_Col] = cleaned['New_ID']
        #update
    #    clipped Class_Col based on the Class_NM_Col
    #    arcpy.RepairGeometry_management(os.path.join(tempfolder,Class_Col+"_dislve.shp"))

    #    arcpy.AddSpatialIndex_management(os.path.join(tempfolder,Class_Col+"_dislve.shp"))
    return cleaned


def save_modified_attributes_to_outputs(mapoldnew_info,tempfolder,OutputFolder,cat_name,riv_name,Path_final_riv,dis_col_name='SubId'):


    NEED_TO_REMOVE_IDS = ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId","SubIdt2"]

    #obtain readme file
    if "DA_Chn_L" in  mapoldnew_info.columns:
        url = 'https://github.com/dustming/RoutingTool/wiki/Files/README_OIH.pdf'
    else:
        url = 'https://github.com/dustming/RoutingTool/wiki/Files/README_NA.pdf'

    response = requests.get(url)
    with open(os.path.join(OutputFolder,"README.pdf"), 'wb') as f:
        f.write(response.content)

    if riv_name != '#':

        if Path_final_riv != '#':
            riv_pd = geopandas.read_file(Path_final_riv)
            riv_pd['Old_SubId'] = riv_pd['SubId']

        cat_pd = mapoldnew_info.drop(columns = 'geometry').copy(deep=True)
        # remove all columns
        if Path_final_riv != '#':
            riv_pd = riv_pd[['geometry','Old_SubId']]
            riv_pd = pd.merge(riv_pd, cat_pd, on='Old_SubId', how='left')
            riv_pd = riv_pd.dissolve(by=dis_col_name, aggfunc='first',as_index=False)


        mapoldnew_info = mapoldnew_info.dissolve(by=dis_col_name, aggfunc='first',as_index=False)
        mapoldnew_info = add_centroid_in_wgs84(mapoldnew_info,"centroid_x","centroid_y")

        cat_c_x_y = mapoldnew_info[["centroid_y","centroid_x"]].copy(deep=True)
        if Path_final_riv != '#':
            riv_pd = riv_pd.drop(columns = ["centroid_y","centroid_x"])
            riv_pd = riv_pd.join(cat_c_x_y)

        riv_pd_nncls_routing_info = mapoldnew_info[mapoldnew_info['Lake_Cat'] != 2][['SubId','DowSubId']].copy(deep=True)
        remove_channel = []
        for subid in riv_pd_nncls_routing_info['SubId'].values:
            if subid not in riv_pd_nncls_routing_info['DowSubId'].values:
                remove_channel.append(subid)
        if Path_final_riv != '#':
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

def clean_geometry_purepy(data,set_precision = -1):

#    data["geometry"] = data["geometry"].apply(lambda x: shapely.wkt.loads(shapely.wkt.dumps(x, rounding_precision=4)))
    if set_precision > 0:
        data['geometry'] = data['geometry'].buffer(0.00000000001)

#        data.geometry = pg.set_precision(data.geometry.values.data, 1e-6)
#        data["geometry"] = data["geometry"].apply(lambda x: shapely.wkt.loads(shapely.wkt.dumps(x, rounding_precision=4)))
#        data['geometry'] = data['geometry'].buffer(0)
#    print("aaaa")
    narow = ~data['geometry'].isna()
#    print("a1",len(data.loc[narow]))
    emrow = ~data.is_empty
#    print("a2",len(data.loc[emrow]))
    arearow = data.area > 0
#    print("a3",len(data.loc[arearow]))
#    arevalid = data.is_valid
    row1 = np.logical_and(narow,emrow)
    rowselect = np.logical_and(arearow,row1)
#    print("a",len(data))
#    rowselect = np.logical_and(rowselect,arevalid)
    data = data.loc[rowselect]
#    print("b",len(data))
    data = data.loc[data.geom_type != 'Point']
    if len(data.loc[data.geom_type == 'GeometryCollection']) > 0:
        print("###########################")
        print("check the following features")
        print(data.loc[data.geom_type == 'GeometryCollection'])
        print("###########################")

    data = data.loc[data.geom_type != 'GeometryCollection']
    data.sindex
#    print("c",len(data))
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

    if "finalcat_info" not in Input_Polygon_path:
#        print(Input_Polygon_path)
        return

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
            for idx in input_wgs_84.index:
                input_wgs_84.loc[idx,'rvhName'] = 'sub' + input_wgs_84.loc[idx,'rvhName']


        input_tojson = input_wgs_84

        for TOLERANCE in TOLERANCEs:
            input_tojson['geometry'] = input_tojson.simplify(TOLERANCE)
            input_tojson.to_file(output_jason_path, driver="GeoJSON")

            json_file_size = os.stat(output_jason_path).st_size/1024/1024 #to MB
            if json_file_size <= 100:
                break

    if len(created_jason_files_lake_riv) > 1 and os.stat(os.path.join(product_dir,Output_file_name[0])).st_size/1024/1024 < 500:
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

    return
