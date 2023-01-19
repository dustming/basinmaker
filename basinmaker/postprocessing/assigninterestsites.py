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
    Routing_Product_Folder= '#',
    path_to_interest_sites = '#',
):
    sub_colnm = "SubId"
    down_colnm = "DowSubId"
    DA_colnm = "DrainArea"
    SegID_colnm = "Seg_ID"
    OutputFolder =  Routing_Product_Folder

    tempfolder = os.path.join(
        tempfile.gettempdir(), "basinmaker_inda_" + str(np.random.randint(1, 10000 + 1))
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)

    Path_Catchment_Polygon="#"
    Path_River_Polyline="#"
    Path_Con_Lake_ply="#"
    Path_NonCon_Lake_ply="#"
    Path_obs_gauge_point="#"
    Path_final_cat_ply="#"
    Path_final_cat_riv="#"

    ##define input files from routing prodcut
    for file in os.listdir(Routing_Product_Folder):
        if file.endswith(".shp"):
            if 'catchment_without_merging_lakes' in file:
                Path_Catchment_Polygon = os.path.join(Routing_Product_Folder, file)
                catname = file
            if 'river_without_merging_lakes' in file:
                Path_River_Polyline = os.path.join(Routing_Product_Folder, file)
                rivname = file
            if 'sl_connected_lake' in file:
                Path_Con_Lake_ply = os.path.join(Routing_Product_Folder, file)
            if 'sl_non_connected_lake' in file:
                Path_NonCon_Lake_ply = os.path.join(Routing_Product_Folder, file)
            if 'obs_gauges' in file:
                Path_obs_gauge_point = os.path.join(Routing_Product_Folder, file)
            if 'finalcat_info' in file:
                Path_final_cat_ply = os.path.join(Routing_Product_Folder, file)
            if 'finalcat_info_riv' in file:
                Path_final_cat_riv = os.path.join(Routing_Product_Folder, file)

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
    finalriv_inforiv = geopandas.read_file(Path_final_riv)
    interest_site = geopandas.read_file(path_to_interest_sites)
    interest_site = interest_site.to_crs(finalriv_infoply.crs)
    interest_site = interest_site.sjoin(finalriv_infoply, how="left")
    interest_site_nolake = interest_site[interest_site["Lake_Cat"] == 0].copy(deep=True)
    interest_subids_nolake = interest_site_nolake['SubId'].values

    interest_site_haslake = interest_site[interest_site["Lake_Cat"] != 0].copy(deep=True)

    lake_subid = np.unique(interest_site_haslake['HyLakeId'].values)
    sub_of_interested_lake = finalriv_infoply[finalriv_infoply["HyLakeId"].isin(lake_subid)].copy(deep=True)
    sub_of_interested_lake = sub_of_interested_lake.sort_values(by='DrainArea', ascending=False)
    sub_of_interested_lake_outlet = sub_of_interested_lake.drop_duplicates(subset=['HyLakeId'],keep='first')
    interest_subids_lake = sub_of_interested_lake_outlet["SubId"].values


    interest_subids = np.concatenate([interest_subids_nolake, interest_subids_lake])

    finalriv_infoply["Has_POI"] = 0
    mask = finalriv_infoply["SubId"].isin(interest_subids)
    finalriv_infoply.loc[mask,"Has_POI"] = 1

    finalriv_inforiv["Has_POI"] = 0
    mask = finalriv_inforiv["SubId"].isin(interest_subids)
    finalriv_inforiv.loc[mask,"Has_POI"] = 1

    finalriv_infoply.to_file(os.path.join(OutputFolder,catname))
    finalriv_inforiv.to_file(os.path.join(OutputFolder,rivname))

    return
