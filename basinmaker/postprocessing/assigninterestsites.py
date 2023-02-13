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
            if 'obs_gauges' in file:
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
    cat_ply = finalriv_infoply.drop(columns=['DA_Obs','SRC_obs','Obs_NM']).copy(deep=True)

    if clean_exist_pois:
        finalriv_infoply["Has_POI"] = 0
        finalriv_infoply["SRC_obs"] = "nan"
        finalriv_infoply["Obs_NM"]  = "nan"
        finalriv_infoply["DA_Obs"] = 0
        finalriv_infoply["DA_error"] = -1.2345

    interest_site = geopandas.read_file(path_to_points_of_interest_points)
    interest_site = interest_site.to_crs(cat_ply.crs)
    interest_site = interest_site.sjoin(cat_ply, how="left")
    non_lake = interest_site[interest_site["Type"] == "River"].copy(deep=True)
    lake     = interest_site[interest_site["Type"] == "Lake"].copy(deep=True)

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
    poi_subs = non_lake.append(lake)
    poi_subs = poi_subs.reset_index()

    if Path_obs_gauge_point != "#":
        exist_poi = geopandas.read_file(Path_obs_gauge_point)
        exist_poi = exist_poi[["Obs_NM","geometry"]].copy(deep=True)
    else:
        exist_poi = poi_subs[["Obs_NM","geometry"]].copy(deep=True)

    for index in poi_subs.index:

        trg_SubId = poi_subs.loc[index,"SubId"]
        SRC_obs_in = poi_subs.loc[index,"SRC_obs"]
        Obs_NM_in =  poi_subs.loc[index,"Obs_NM"]
        DA_Obs_in =  poi_subs.loc[index,"DA_Obs"]
        row = poi_subs[poi_subs.index == index]

        trg_sub_info = finalriv_infoply[finalriv_infoply["SubId"] == trg_SubId].copy(deep=True)
        Cur_SRC_obs  =  trg_sub_info.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"].values[0]
        Cur_Obs_NM   =  trg_sub_info.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"].values[0]

        # no point is here
        if len(exist_poi) < 0:
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = SRC_obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"] = Obs_NM_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = DA_Obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            if DA_Obs_in > 0 :
                DAerror_in = finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DrainArea"]/1000/1000/DA_Obs_in
                finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DA_error"] = DAerror_in
            continue

        # the subbasin already contain the ObsNM
        if Obs_NM_in in Cur_Obs_NM:
            print("The POI ",Obs_NM_in,"Already exists")
            continue

        # Obs_NM not in current subbasin but in other subbasins
        if Obs_NM_in in finalriv_infoply["Obs_NM"].values:
            # delelte the old point of interest

            mask_cat  = finalriv_infoply["Obs_NM"] == Obs_NM_in

            exist_poi = exist_poi.reset_index()

            mask_point = exist_poi["Obs_NM"] == Obs_NM_in

            finalriv_infoply.loc[mask_cat,"Obs_NM"] = "nan"
            finalriv_infoply.loc[mask_cat,"SRC_obs"] = "nan"
            finalriv_infoply.loc[mask_cat,"Has_POI"] = 0
            finalriv_infoply.loc[mask_cat,"DA_Obs"] = 0
            finalriv_infoply.loc[mask_cat,"DA_error"] = -1.2345
            exist_poi=exist_poi.drop(exist_poi.index[mask_point])
            print("The existing : ",Obs_NM_in," in the product is removed")

        if trg_sub_info["Has_POI"].values[0] <= 0:
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = SRC_obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"] = Obs_NM_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DA_Obs"] = DA_Obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            if DA_Obs_in > 0 :
                DAerror_in = finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DrainArea"]/1000/1000/DA_Obs_in
                finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"DA_error"] = DAerror_in
            exist_poi = exist_poi.append(row[["Obs_NM","geometry"]])
            print("Add the POI : ",Obs_NM_in," into the routing product ")

        else:

            if SRC_obs_in != "nan":
                SRC_obs_in  =  Cur_SRC_obs+"&" + SRC_obs_in
            else:
                SRC_obs_in  = Cur_SRC_obs

            Obs_NM_in   =   Cur_Obs_NM+"&" + Obs_NM_in

            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"SRC_obs"] = SRC_obs_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Obs_NM"] = Obs_NM_in
            finalriv_infoply.loc[finalriv_infoply["SubId"] == trg_SubId,"Has_POI"] = 1
            row["Obs_NM"] = Obs_NM_in
            exist_poi.loc[exist_poi["Obs_NM"] == Cur_Obs_NM ,"Obs_NM"] = Obs_NM_in
            exist_poi = exist_poi.append(row[["Obs_NM","geometry"]])
            print("Append the POI : ",Obs_NM_in," into the routing product ")

    finalriv_infoply.to_file(os.path.join(OutputFolder,os.path.basename(Path_final_riv_ply)))

    cat_table = finalriv_infoply.drop(columns="geometry")
    exist_poi = exist_poi.merge(cat_table,on='Obs_NM',how='left')
    exist_poi = exist_poi[["geometry",'DA_Obs','SRC_obs','Obs_NM','SubId','DrainArea']].copy(deep=True)

    if Path_obs_gauge_point != "#":
        exist_poi.to_file(os.path.join(OutputFolder,os.path.basename(Path_obs_gauge_point)))
    else:
        exist_poi.to_file(os.path.join(OutputFolder,"obs_gauges"+"_v1-0"+".shp"))


    if Path_Con_Lake_ply != "#":
        con_lake = geopandas.read_file(Path_Con_Lake_ply)
        con_lake.to_file(os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply)))

    if Path_NonCon_Lake_ply != "#":
        noncon_lake = geopandas.read_file(Path_NonCon_Lake_ply)
        noncon_lake.to_file(os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)))

    if Path_River_Polyline != "#":
        riv = geopandas.read_file(Path_River_Polyline)
        riv = riv[["SubId","geometry"]]
        riv = riv.merge(cat_table,on='SubId',how='left')
        riv.to_file(os.path.join(OutputFolder,os.path.basename(Path_River_Polyline)))

    return
