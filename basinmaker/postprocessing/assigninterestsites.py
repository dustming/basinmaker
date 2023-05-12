import numpy as np
import sys
import os
import csv
import tempfile
import copy
import pandas as pd
import shutil
import geopandas

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from basinmaker.func.purepy import *
from basinmaker.func.pdtable import *


def define_interest_sites(
    routing_product_folder= '#',
    path_to_points_of_interest_points = '#',
    clean_exist_pois = True,
    path_output_folder = "#",
):
    sub_colnm = "SubId"
    down_colnm = "DowSubId"
    DA_colnm = "DrainArea"
    SegID_colnm = "Seg_ID"
    OutputFolder =  path_output_folder

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

    Path_Catchment_Polygon="#"
    Path_River_Polyline="#"
    Path_Con_Lake_ply="#"
    Path_NonCon_Lake_ply="#"
    Path_obs_gauge_point="#"
    Path_final_cat_ply="#"
    Path_final_cat_riv="#"

    ##define input files from routing prodcut
    for file in os.listdir(routing_product_folder):
        if file.endswith(".shp"):
            if 'catchment_without_merging_lakes' in file:
                Path_Catchment_Polygon = os.path.join(routing_product_folder, file)
                catname = file
            if 'river_without_merging_lakes' in file:
                Path_River_Polyline = os.path.join(routing_product_folder, file)
                rivname = file
            if 'sl_connected_lake' in file:
                Path_Con_Lake_ply = os.path.join(routing_product_folder, file)
            if 'sl_non_connected_lake' in file:
                Path_NonCon_Lake_ply = os.path.join(routing_product_folder, file)
            if 'obs_gauges' in file or 'poi' in file:
                Path_obs_gauge_point = os.path.join(routing_product_folder, file)

    if Path_Catchment_Polygon == '#' or  Path_River_Polyline =='#':
        print("Invalid routing product folder ")

    Path_final_riv_ply = Path_Catchment_Polygon
    Path_final_riv = Path_River_Polyline

    # read attribute table, and
    Path_final_rviply = Path_final_riv_ply
    Path_final_riv = Path_final_riv
    Path_Conl_ply = Path_Con_Lake_ply
    Path_Non_ConL_ply = Path_NonCon_Lake_ply

    finalriv_infoply = geopandas.read_file(Path_final_rviply)

    # For ECCC and other product we still use DA_error
    DA_ERR_COL = 'DA_error'
    if DA_ERR_COL not in finalriv_infoply.columns:
        DA_ERR_COL = 'DA_Diff'

    # Remove DA_Obs and DA_Diff columns from subbasin layer
    for col in ['DA_Obs', DA_ERR_COL]:
        if col in finalriv_infoply.columns:
            finalriv_infoply = finalriv_infoply.drop(columns=[col])
        else:
            pass

    cat_ply = finalriv_infoply.drop(columns=['Obs_NM', 'SRC_obs']).copy(deep=True)

    if clean_exist_pois:
        finalriv_infoply["Has_POI"] = 0
        finalriv_infoply["SRC_obs"] = "nan"
        finalriv_infoply["Obs_NM"]  = "nan"
        # finalriv_infoply["DA_Obs"] = -1.2345
        # finalriv_infoply[DA_ERR_COL] = -999

    interest_site = geopandas.read_file(path_to_points_of_interest_points)
    interest_site = interest_site.to_crs(cat_ply.crs)
    interest_site = interest_site.sjoin(cat_ply, how="left")
    non_lake = interest_site[interest_site["Type"] == "River"].copy(deep=True)
    lake     = interest_site[interest_site["Type"] == "Lake"].copy(deep=True)

    # print(interest_site.columns)

    if len(non_lake) > 0:
        if len(non_lake[non_lake["Lake_Cat"] == 1]) > 0:
            print("Folllowing river gauges are located in the lake, they will be ignored")
            print(non_lake[non_lake["Lake_Cat"] == 1]["Obs_NM"].values)
            non_lake = non_lake[non_lake["Lake_Cat"] == 0].copy(deep=True)

    if len(lake) > 0:
        if len(lake[lake["Lake_Cat"] == 0]) > 0:
            print("Folllowing lake gauges are not located in the lake, they will be ignored")
            print(lake[lake["Lake_Cat"] == 0]["Obs_NM"].values)
            lake = lake[lake["Lake_Cat"] == 1].copy(deep=True)

    # obtain subid of lake outlet
    lake_subid = np.unique(lake['HyLakeId'].values)
    sub_of_interested_lake = finalriv_infoply[finalriv_infoply["HyLakeId"].isin(lake_subid)].copy(deep=True)
    sub_of_interested_lake = sub_of_interested_lake.sort_values(by='DrainArea', ascending=False)
    sub_of_interested_lake_outlet = sub_of_interested_lake.drop_duplicates(subset=['HyLakeId'],keep='first')
    interest_subids_lake = sub_of_interested_lake_outlet[["SubId","HyLakeId"]]
    interest_subids_lake["Lake_Outlet"] = interest_subids_lake["SubId"]
    interest_subids_lake = interest_subids_lake.drop(columns='SubId')
    lake = lake.merge(interest_subids_lake,on='HyLakeId',how='left')
    lake["SubId"] = lake["Lake_Outlet"]
    lake = lake.drop(columns='Lake_Outlet')
    poi_subs = pd.concat([non_lake, lake])
    poi_subs = poi_subs.reset_index()

    if Path_obs_gauge_point != "#":
        exist_poi = geopandas.read_file(Path_obs_gauge_point)
        exist_poi = exist_poi[["Obs_NM","geometry",'DA_Obs','SRC_obs','DrainArea']].copy(deep=True)
        cat_table = finalriv_infoply[['Obs_NM','SubId']].copy(deep=True)
        exist_poi = exist_poi.merge(cat_table,on='Obs_NM',how='left')
    else:
        exist_poi = poi_subs[["Obs_NM","geometry",'DA_Obs','SRC_obs','DrainArea','SubId']].copy(deep=True)

    for index in poi_subs.index:

        trg_SubId = poi_subs.loc[index,"SubId"]
        SRC_obs_in = poi_subs.loc[index,"SRC_obs"]
        Obs_NM_in =  poi_subs.loc[index,"Obs_NM"]
        # DA_Obs_in =  poi_subs.loc[index,"DA_Obs"]
        row = poi_subs[poi_subs.index == index]
        if np.isnan(trg_SubId):
            print("POI: ",Obs_NM_in,"   is not located in any subbasin, please check the location of this POI")
            continue
        trg_sub_info = finalriv_infoply[finalriv_infoply["SubId"] == trg_SubId].copy(deep=True)
        Cur_SRC_obs  =  trg_sub_info.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"].values[0]
        Cur_Obs_NM   =  trg_sub_info.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"].values[0]

        # Initialze this variable
        # DAerror_in   = -999    
        # no point is here
        if len(exist_poi) < 0:
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = SRC_obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"] = Obs_NM_in
            # finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DA_Obs"] = DA_Obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            # if DA_Obs_in > 0 :
            #     da_sim = finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DrainArea"].values[0]
            #     if DA_ERR_COL == 'DA_error':
            #         DAerror_in = da_sim/1000/1000/DA_Obs_in
            #     else:
            #         DAerror_in = (da_sim/1000/1000 - DA_Obs_in)/DA_Obs_in
            #     finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,DA_ERR_COL] = DAerror_in
            # continue

        # the subbasin already contain the ObsNM
        if Obs_NM_in in Cur_Obs_NM:
            print("The POI ",Obs_NM_in,"Already exists")
            continue

        # Obs_NM not in current subbasin but in other subbasins
        if Obs_NM_in in finalriv_infoply["Obs_NM"].values:
            # delelte the old point of interest

            mask_cat  = finalriv_infoply["Obs_NM"] == Obs_NM_in

            exist_poi = exist_poi.reset_index(drop=True)

            mask_point = exist_poi["Obs_NM"] == Obs_NM_in

            finalriv_infoply.loc[mask_cat,"Obs_NM"] = "nan"
            finalriv_infoply.loc[mask_cat,"SRC_obs"] = "nan"
            finalriv_infoply.loc[mask_cat,"Has_POI"] = 0
            # finalriv_infoply.loc[mask_cat,"DA_Obs"] = -1.2345
            # finalriv_infoply.loc[mask_cat,DA_ERR_COL] = -999
            exist_poi=exist_poi.drop(exist_poi.index[mask_point])
            print("The existing : ",Obs_NM_in," in the product is removed")

        if trg_sub_info["Has_POI"].values[0] <= 0:
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = SRC_obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"] = Obs_NM_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            
            # finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DA_Obs"] = DA_Obs_in
            # finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            # if DA_Obs_in > 0 :
            #     da_sim = finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DrainArea"].values[0]
            #     if DA_ERR_COL == 'DA_error':
            #         DAerror_in = da_sim/1000/1000/DA_Obs_in
            #     else:
            #         DAerror_in = (da_sim/1000/1000 - DA_Obs_in)/DA_Obs_in
            #     finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,DA_ERR_COL] = DAerror_in


            exist_poi = pd.concat([exist_poi, row[["Obs_NM","geometry",'DA_Obs','SRC_obs','DrainArea','SubId']]])
            print("Add the POI : ",Obs_NM_in," into the routing product ")

        else:
            ## POI in this subbasin

            SRC_obs_in  =  Cur_SRC_obs+ "&" + SRC_obs_in
            Obs_NM_in   =   Cur_Obs_NM + "&" + Obs_NM_in


            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = SRC_obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"] = Obs_NM_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1

            # finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DA_Obs"] = DA_Obs_in
            # finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            # if DA_Obs_in > 0 :
            #     da_sim = finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DrainArea"].values[0]
            #     if DA_ERR_COL == 'DA_error':
            #         DAerror_in = da_sim/1000/1000/DA_Obs_in
            #     else:
            #         DAerror_in = (da_sim/1000/1000 - DA_Obs_in)/DA_Obs_in

            #     finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,DA_ERR_COL] = DAerror_in

            # else:
            #     finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,DA_ERR_COL] = -999

            # row["Obs_NM"] = Obs_NM_in
            # exist_poi.loc[exist_poi["Obs_NM"] == Cur_Obs_NM ,"Obs_NM"] = Obs_NM_in

            exist_poi = pd.concat([exist_poi, row[["Obs_NM","geometry",'DA_Obs','SRC_obs','DrainArea','SubId']]])

            print("Append POI : ", Obs_NM_in , ',' ,  SRC_obs_in )
    
    # Output subbasin layer
    finalriv_infoply.to_file(os.path.join(OutputFolder,os.path.basename(Path_final_riv_ply)))

    # cat_table = finalriv_infoply.drop(columns=["geometry",'DA_Obs','SRC_obs','DrainArea'])
    # exist_poi = exist_poi.merge(cat_table,on='Obs_NM',how='left')
    exist_poi['Use_region'] = 1
    exist_poi[DA_ERR_COL]   = -999
    mask = exist_poi["DA_Obs"] > 0 

    if DA_ERR_COL == 'DA_error':
        exist_poi.loc[mask,DA_ERR_COL] = exist_poi.loc[mask,'DrainArea']/1000/1000/exist_poi.loc[mask,'DA_Obs']                                                                           
    else:
        exist_poi.loc[mask,DA_ERR_COL] = (exist_poi.loc[mask,'DrainArea']/1000/1000 - exist_poi.loc[mask,'DA_Obs'])/exist_poi.loc[mask,'DA_Obs']                                                                          

    exist_poi = exist_poi[['SubId','Obs_NM','DA_Obs','DrainArea',DA_ERR_COL,'SRC_obs','Use_region','geometry']].copy(deep=True)

    # Revise DA_Diff column type and value format
    exist_poi.loc[exist_poi['DA_Obs'] <= 0,DA_ERR_COL] = -999
    # non-null values are shown as percentage
    if DA_ERR_COL == 'DA_Diff':
        exist_poi[DA_ERR_COL] = exist_poi[DA_ERR_COL].apply(lambda x: f"{x*100:.3f}%" if x != -999 else "<NA>")
    exist_poi.loc[exist_poi['DA_Obs'] <= 0, 'DA_Obs'] = -1.2345

    if Path_obs_gauge_point != "#":
        exist_poi.to_file(os.path.join(OutputFolder,os.path.basename(Path_obs_gauge_point)))
    else:
        exist_poi.to_file(os.path.join(OutputFolder,"poi"+"_v1-0"+".shp"))


    if Path_Con_Lake_ply != "#":
        con_lake = geopandas.read_file(Path_Con_Lake_ply)
        con_lake.to_file(os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply)))

    if Path_NonCon_Lake_ply != "#":
        noncon_lake = geopandas.read_file(Path_NonCon_Lake_ply)
        noncon_lake.to_file(os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)))

    if Path_River_Polyline != "#":
        riv = geopandas.read_file(Path_River_Polyline)
        riv = riv[["SubId","geometry"]]
        cat_table = finalriv_infoply.drop(columns=['geometry'])
        riv = riv.merge(cat_table,on='SubId',how='left')
        riv.to_file(os.path.join(OutputFolder,os.path.basename(Path_River_Polyline)))

    return
