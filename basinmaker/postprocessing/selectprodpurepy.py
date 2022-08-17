
import numpy as np
import geopandas
import sys
import os
import csv
import pandas as pd
import tempfile 
import copy 

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from basinmaker.func.purepy import *
from basinmaker.func.pdtable import *



def Select_Routing_product_based_SubId_purepy(
    OutputFolder,
    Routing_Product_Folder,
    mostdownid=-1,
    mostupstreamid=-1,
):
    """Extract region of interest based on provided Subid

    Function that used to obtain the region of interest from routing
    product based on given SubId

    Parameters
    ----------
    OutputFolder                   : string
        Folder path that stores extracted routing product
    Path_Catchment_Polygon         : string
        Path to the catchment polygon
    Path_River_Polyline            : string (optional)
        Path to the river polyline
    Path_Con_Lake_ply              : string (optional)
        Path to a connected lake polygon. Connected lakes are lakes that
        are connected by Path_final_cat_riv or Path_final_riv.
    Path_NonCon_Lake_ply           : string (optional)
        Path to a non connected lake polygon. Connected lakes are lakes
        that are not connected by Path_final_cat_riv or Path_final_riv.
    mostdownid                     : integer
        It is the most downstream subbasin ID in the region of interest
    mostupstreamid                 : integer (optional)
        It is the most upstream subbasin ID in the region of interest.
        Normally it is -1, indicating all subbasin drainage to mostdownid
        is needed. In some case, if not all subbasin drainage to mostdownid
        is needed, then the most upstream subbsin ID need to be provided
        here.


    Notes
    -------
    This function has no return values, instead following fiels will be
    generated. The output files have are same as inputs expect the extent
    are different.

    os.path.join(OutputFolder,os.path.basename(Path_Catchment_Polygon))
    os.path.join(OutputFolder,os.path.basename(Path_River_Polyline))
    os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply))
    os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply))

    Returns:
    -------
    None

    Examples
    -------

    """
    tempfolder = os.path.join(
        tempfile.gettempdir(),
        "basinmaker_locsubid" + str(101),#np.random.randint(1, 10000 + 1)),
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
            if 'river_without_merging_lakes' in file:
                Path_River_Polyline = os.path.join(Routing_Product_Folder, file)
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
        return()



    sub_colnm = "SubId"
    down_colnm = "DowSubId"

    ##3
     
    cat_ply = geopandas.read_file(Path_Catchment_Polygon)

    hyshdinfo = cat_ply[[sub_colnm, down_colnm]].astype("int32").values

    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in cat_ply.columns:
        Gauge_col_Name = "Has_Gauge"
        
    ### Loop for each downstream id

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

    ## find all subid control by this subid
    for i_down in range(0,len(mostdownid)):
        ### Loop for each downstream id
        OutHyID = mostdownid[i_down]
        OutHyID2 = mostupstreamid[i_down]
            
        ## find all subid control by this subid
        HydroBasins1 = defcat(hyshdinfo, OutHyID)
        if OutHyID2 > 0:
            HydroBasins2 = defcat(hyshdinfo, OutHyID2)
            ###  exculde the Ids in HydroBasins2 from HydroBasins1
            for i in range(len(HydroBasins2)):
                rows = np.argwhere(HydroBasins1 == HydroBasins2[i])
                HydroBasins1 = np.delete(HydroBasins1, rows)
            HydroBasins = HydroBasins1
        else:
            HydroBasins = HydroBasins1

        if i_down == 0:
            HydroBasins_All = HydroBasins
        else:
            HydroBasins_All = np.concatenate((HydroBasins_All, HydroBasins), axis=0)

    Outputfilename_cat = os.path.join(
        OutputFolder, os.path.basename(Path_Catchment_Polygon)
    )

    cat_ply_select = cat_ply.loc[cat_ply['SubId'].isin(HydroBasins_All)]

    cat_ply_select.to_file(Outputfilename_cat) 
     
    Outputfilename_cat_riv = os.path.join(
        OutputFolder, os.path.basename(Path_River_Polyline)
    )

    cat_riv = geopandas.read_file(Path_River_Polyline)
    

    cat_riv_select = cat_riv.loc[cat_riv['SubId'].isin(HydroBasins)]
    
    cat_riv_select.to_file(Outputfilename_cat_riv) 
    
    cat_ply_select = geopandas.read_file(Outputfilename_cat)
    Connect_Lake_info = cat_ply_select.loc[cat_ply_select["Lake_Cat"] == 1]
    Connect_Lakeids = np.unique(Connect_Lake_info["HyLakeId"].values)
    Connect_Lakeids = Connect_Lakeids[Connect_Lakeids > 0]
    
    NConnect_Lake_info = cat_ply_select.loc[cat_ply_select["Lake_Cat"] == 2]
    NonCL_Lakeids = np.unique(NConnect_Lake_info["HyLakeId"].values)
    NonCL_Lakeids = NonCL_Lakeids[NonCL_Lakeids > 0]

    if len(Connect_Lakeids) > 0 and Path_Con_Lake_ply != "#":

        sl_con_lakes = geopandas.read_file(Path_Con_Lake_ply)
        sl_con_lakes = sl_con_lakes.loc[sl_con_lakes['Hylak_id'].isin(Connect_Lakeids)]
        sl_con_lakes.to_file(os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply)))


    if len(NonCL_Lakeids) > 0 and Path_NonCon_Lake_ply != "#":
        sl_non_con_lakes = geopandas.read_file(Path_NonCon_Lake_ply)
        sl_non_con_lakes = sl_non_con_lakes.loc[sl_non_con_lakes['Hylak_id'].isin(NonCL_Lakeids)]
        sl_non_con_lakes.to_file(os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)))
    
    sl_gauge_info = cat_ply_select.loc[cat_ply_select[Gauge_col_Name] > 0]
    sl_gauge_nm = np.unique(sl_gauge_info["Obs_NM"].values)
    sl_gauge_nm = sl_gauge_nm[sl_gauge_nm != 'nan']
    if len(sl_gauge_nm) > 0 and Path_obs_gauge_point !='#':
        all_gauge = geopandas.read_file(Path_obs_gauge_point)
        sl_gauge = all_gauge.loc[all_gauge['Obs_NM'].isin(sl_gauge_nm)]
        sl_gauge.to_file(os.path.join(OutputFolder,os.path.basename(Path_obs_gauge_point)))
    
    return 


################################################################