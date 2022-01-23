
import numpy as np
import arcpy
from arcpy import env
from arcpy.sa import *
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import sys
import os
import csv
import pandas as pd
import tempfile 
import copy 

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")


def Locate_subid_needsbyuser_arcgis(
    Path_Points="#", Gauge_NMS="#", Path_products="#"
):
    """Get Subbasin Ids

    Function that used to obtain subbasin ID of certain gauge.
    or subbasin ID of the polygon that includes the given point
    shapefile.

    Parameters
    ----------
    Path_Points      : string (Optional)
        It is the path of the point shapefile. If the point shapefile is
        provided. The function will return subids of those catchment
        polygons that includes these point in the point shapefile

    Gauge_NMS        : list
        Name of the streamflow gauges, such as ['09PC019'], if the gauge
        name is provided, the subbasin ID that contain this gauge will be
        returned
    Path_products    : string
        The path of the subbasin polygon shapefile.
        The shapefile should at least contains following columns
        ##############Subbasin related attributes###########################
        SubID           - integer, The subbasin Id
        DowSubId        - integer, The downstream subbasin ID of this
                                   subbasin
        Obs_NM          - The streamflow obervation gauge name.

    Notes
    -------
    Path_Points or Gauge_NMS should only provide one each time
    to use this function

    Returns:
    -------
    SubId_Selected  : list
        It is a list contains the selected subid based on provided
        streamflow gauge name or provided point shapefile

    Examples
    -------

    """

    tempfolder = os.path.join(
        tempfile.gettempdir(),
        "basinmaker_locsubid" + str(np.random.randint(1, 10000 + 1)),
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)
    arcpy.env.workspace = tempfolder

    SubId_Selected = -1
    if Gauge_NMS[0] != "#":

        hyshdinfo2 = pd.DataFrame.spatial.from_featureclass(Path_products)
        hyshdinfo2 = hyshdinfo2.loc[hyshdinfo2["Obs_NM"] != "-9999.0"]
        hyshdinfo2 = hyshdinfo2.loc[hyshdinfo2["Obs_NM"].isin(Gauge_NMS)]
        hyshdinfo2 = hyshdinfo2[["Obs_NM", "SubId"]]
        #            hyshdinfo2.to_csv(os.path.join(self.OutputFolder,'SubIds_Selected.csv'),sep=',', index = None)
        SubId_Selected = hyshdinfo2["SubId"].values

    if Path_Points != "#":
        SpRef_in = arcpy.Describe(Path_products).spatialReference

        arcpy.Project_management(Path_products,"Obspoint_project2.shp", out_coordinate_system)



        arcpy.SpatialJoin_analysis("Obspoint_project2.shp", Path_products, 'Sub_Selected_by_Points')

        hyshdinfo2 = Outputfilename_cat(os.path.join(tempfolder, "Sub_Selected_by_Points.shp"))

        SubId_Selected = hyshdinfo2["SubId"].values
        SubId_Selected = SubId_Selected[SubId_Selected > 0]

    return SubId_Selected


def Select_Routing_product_based_SubId_arcgis(
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
    arcpy.env.workspace = tempfolder


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
        arcpy.AddMessage(Path_Catchment_Polygon)
        return()



    sub_colnm = "SubId"
    down_colnm = "DowSubId"

    ##3
    
    cat_ply = pd.DataFrame.spatial.from_featureclass(Path_Catchment_Polygon)

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

    cat_ply_select.spatial.to_featureclass(location=Outputfilename_cat,overwrite=True,sanitize_columns=False) 
    Outputfilename_cat_riv = os.path.join(
        OutputFolder, os.path.basename(Path_River_Polyline)
    )

    cat_riv = pd.DataFrame.spatial.from_featureclass(Path_River_Polyline)
    

    cat_riv_select = cat_riv.loc[cat_riv['SubId'].isin(HydroBasins)]
    
    cat_riv_select.spatial.to_featureclass(location=Outputfilename_cat_riv,overwrite=True,sanitize_columns=False) 
    
    cat_ply_select = pd.DataFrame.spatial.from_featureclass(Outputfilename_cat)
    Connect_Lake_info = cat_ply_select.loc[cat_ply_select["Lake_Cat"] == 1]
    Connect_Lakeids = np.unique(Connect_Lake_info["HyLakeId"].values)
    Connect_Lakeids = Connect_Lakeids[Connect_Lakeids > 0]
    
    NConnect_Lake_info = cat_ply_select.loc[cat_ply_select["Lake_Cat"] == 2]
    NonCL_Lakeids = np.unique(NConnect_Lake_info["HyLakeId"].values)
    NonCL_Lakeids = NonCL_Lakeids[NonCL_Lakeids > 0]

    if len(Connect_Lakeids) > 0 and Path_Con_Lake_ply != "#":

        sl_con_lakes = pd.DataFrame.spatial.from_featureclass(Path_Con_Lake_ply)
        sl_con_lakes = sl_con_lakes.loc[sl_con_lakes['Hylak_id'].isin(Connect_Lakeids)]
        sl_con_lakes.spatial.to_featureclass(location=os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply)),overwrite=True,sanitize_columns=False)


    if len(NonCL_Lakeids) > 0 and Path_NonCon_Lake_ply != "#":
        sl_non_con_lakes = pd.DataFrame.spatial.from_featureclass(Path_NonCon_Lake_ply)
        sl_non_con_lakes = sl_non_con_lakes.loc[sl_non_con_lakes['Hylak_id'].isin(NonCL_Lakeids)]
        sl_non_con_lakes.spatial.to_featureclass(location=os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)),overwrite=True,sanitize_columns=False)
    
    sl_gauge_info = cat_ply_select.loc[cat_ply_select[Gauge_col_Name] > 0]
    sl_gauge_nm = np.unique(sl_gauge_info["Obs_NM"].values)
    sl_gauge_nm = sl_gauge_nm[sl_gauge_nm != 'nan']
    if len(sl_gauge_nm) > 0 and Path_obs_gauge_point !='#':
        all_gauge = pd.DataFrame.spatial.from_featureclass(Path_obs_gauge_point)
        sl_gauge = all_gauge.loc[all_gauge['Obs_NM'].isin(sl_gauge_nm)]
        sl_gauge.spatial.to_featureclass(location=os.path.join(OutputFolder,os.path.basename(Path_obs_gauge_point)),overwrite=True,sanitize_columns=False)
    
    return 


################################################################