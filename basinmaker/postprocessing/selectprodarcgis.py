
import numpy as np
import arcpy
from arcpy import env
from arcpy.sa import *
import sys
import os
import csv
import pandas as pd
from arcgis import GIS
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import tempfile 
import copy 
from ..func.pdtable import *
from ..func.arcgis import *

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
    Path_Catchment_Polygon="#",
    Path_River_Polyline="#",
    Path_Con_Lake_ply="#",
    Path_NonCon_Lake_ply="#",
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
        "basinmaker_locsubid" + str(np.random.randint(1, 10000 + 1)),
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)
    arcpy.env.workspace = tempfolder


    sub_colnm = "SubId"
    down_colnm = "DowSubId"

    ##3
    hyshdinfo2 = pd.DataFrame.spatial.from_featureclass(Path_Catchment_Polygon)

    hyshdinfo = hyshdinfo2[[sub_colnm, down_colnm]].astype("int32").values

    ### Loop for each downstream id
    OutHyID = mostdownid
    OutHyID2 = mostupstreamid

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

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

    Outputfilename_cat = os.path.join(
        OutputFolder, os.path.basename(Path_Catchment_Polygon)
    )

    arcpy.AddMessage(hyshdinfo)
    arcpy.AddMessage(OutHyID)
    arcpy.AddMessage(OutHyID2)
    arcpy.AddMessage(OutHyID in hyshdinfo[:,0])
    arcpy.AddMessage(OutHyID in hyshdinfo[:,1])
    arcpy.AddMessage(HydroBasins)
        
    select_feature_by_attributes_arcgis(
        input = Path_Catchment_Polygon, 
        Attri_NM = "SubId",
        Attri_v = HydroBasins,
        output = Outputfilename_cat,
    )

    if Path_River_Polyline != "#":
        Outputfilename_cat_riv = os.path.join(
            OutputFolder, os.path.basename(Path_River_Polyline)
        )

        select_feature_by_attributes_arcgis(
            input = Path_River_Polyline, 
            Attri_NM = "SubId",
            Attri_v = HydroBasins,
            output = Outputfilename_cat_riv,
        )

    finalcat_info = pd.DataFrame.spatial.from_featureclass(Outputfilename_cat)

    Connect_Lake_info = finalcat_info.loc[finalcat_info["IsLake"] == 1]
    Connect_Lakeids = np.unique(Connect_Lake_info["HyLakeId"].values)
    Connect_Lakeids = Connect_Lakeids[Connect_Lakeids > 0]

    NConnect_Lake_info = finalcat_info.loc[finalcat_info["IsLake"] == 2]
    NonCL_Lakeids = np.unique(NConnect_Lake_info["HyLakeId"].values)
    NonCL_Lakeids = NonCL_Lakeids[NonCL_Lakeids > 0]

    if len(Connect_Lakeids) > 0 and Path_Con_Lake_ply != "#":

        select_feature_by_attributes_arcgis(
            input = Path_Con_Lake_ply, 
            Attri_NM = "Hylak_id",
            Attri_v = Connect_Lakeids,
            output = os.path.join(OutputFolder, os.path.basename(Path_Con_Lake_ply)),
        )

    if len(NonCL_Lakeids) > 0 and Path_NonCon_Lake_ply != "#":
        select_feature_by_attributes_arcgis(
            input = Path_NonCon_Lake_ply, 
            Attri_NM = "Hylak_id",
            Attri_v = NonCL_Lakeids,
            output = os.path.join(OutputFolder, os.path.basename(Path_NonCon_Lake_ply)),
        )
    return 


################################################################


OutputFolder =sys.argv[6]
Path_Catchment_Polygon=sys.argv[1]
Path_River_Polyline=sys.argv[2]
Path_Con_Lake_ply=sys.argv[3]
Path_NonCon_Lake_ply=sys.argv[4]
mostdownid=int(sys.argv[5])
mostupstreamid=int(sys.argv[6])

Select_Routing_product_based_SubId_arcgis(
    OutputFolder = OutputFolder,
    Path_Catchment_Polygon=Path_Catchment_Polygon,
    Path_River_Polyline=Path_River_Polyline,
    Path_Con_Lake_ply=Path_Con_Lake_ply,
    Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
    mostdownid=mostdownid,
    mostupstreamid=mostupstreamid,
)


# gis = GIS()
# item = gis.content.get(Path_Catchment_Polygon)
# flayer = item.layers[0]
# sdf = pd.DataFrame.spatial.from_layer(flayer)
# sdf.head()
