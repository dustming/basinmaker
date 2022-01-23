from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import pandas as pd
import numpy as np
import tempfile


def Locate_subid_needsbyuser_qgis(
    Path_Points="#", Gauge_NMS="#", Path_products="#", qgis_prefix_path="#"
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

    # obtain subbasin ID based on either points or guage names
    QgsApplication.setPrefixPath(qgis_prefix_path, True)
    Qgs = QgsApplication([], False)
    Qgs.initQgis()
    from processing.core.Processing import Processing
    from processing.tools import dataobjects
    from qgis import processing

    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    SubId_Selected = -1
    if Gauge_NMS[0] != "#":
        hyinfocsv = Path_products[:-3] + "dbf"
        tempinfo = Dbf5(hyinfocsv)
        hyshdinfo2 = tempinfo.to_dataframe().drop_duplicates("SubId", keep="first")
        hyshdinfo2 = hyshdinfo2.loc[hyshdinfo2["Obs_NM"] != "-9999.0"]
        hyshdinfo2 = hyshdinfo2.loc[hyshdinfo2["Obs_NM"].isin(Gauge_NMS)]
        hyshdinfo2 = hyshdinfo2[["Obs_NM", "SubId"]]
        #            hyshdinfo2.to_csv(os.path.join(self.OutputFolder,'SubIds_Selected.csv'),sep=',', index = None)
        SubId_Selected = hyshdinfo2["SubId"].values

    if Path_Points != "#":
        vector_layer = qgis_vector_read_vector(processing, context, Path_products)
        SpRef_in = qgis_vector_return_crs_id(processing, context, vector_layer)

        qgis_vector_reproject_layers(
            processing,
            context,
            INPUT=Path_Points,
            TARGET_CRS=SpRef_in,
            OUTPUT=os.path.join(tempfolder, "Obspoint_project2.shp"),
        )

        qgis_vector_add_polygon_attribute_to_points(
            processing,
            context,
            INPUT_Layer=os.path.join(tempfolder, "Obspoint_project2.shp"),
            POLYGONS=Path_products,
            FIELDS="SubId",
            OUTPUT=os.path.join(tempfolder, "Sub_Selected_by_Points.shp"),
        )

        hyshdinfo2 = Dbf_To_Dataframe(
            os.path.join(tempfolder, "Sub_Selected_by_Points.shp")
        )
        SubId_Selected = hyshdinfo2["SubId"].values
        SubId_Selected = SubId_Selected[SubId_Selected > 0]

    Qgs.exit()

    return SubId_Selected


def Select_Routing_product_based_SubId_qgis(
    OutputFolder,
    # Path_Catchment_Polygon="#",
    # Path_River_Polyline="#",
    # Path_Con_Lake_ply="#",
    # Path_NonCon_Lake_ply="#",
    mostdownid=-1,
    mostupstreamid=-1,
    Routing_Product_Folder = '#',
    qgis_prefix_path="#",
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
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)
    tempfolder = os.path.join(
        tempfile.gettempdir(),
        "basinmaker_extsubprod" + str(np.random.randint(1, 10000 + 1)),
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)

    QgsApplication.setPrefixPath(qgis_prefix_path, True)
    Qgs = QgsApplication([], False)
    Qgs.initQgis()
    from processing.core.Processing import Processing
    from qgis import processing

    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())


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
    
    sub_colnm = "SubId"
    down_colnm = "DowSubId"
    ##3
    hyshdinfo2 = Dbf_To_Dataframe(Path_Catchment_Polygon).drop_duplicates(
        sub_colnm, keep="first"
    )
    hyshdinfo = hyshdinfo2[[sub_colnm, down_colnm]].astype("int32").values


    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in hyshdinfo2.columns:
        Gauge_col_Name = "Has_Gauge"
        
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)
        
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
    Selectfeatureattributes(
        processing,
        Input=Path_Catchment_Polygon,
        Output=Outputfilename_cat,
        Attri_NM="SubId",
        Values=HydroBasins_All,
    )
    if Path_River_Polyline != "#":
        Outputfilename_cat_riv = os.path.join(
            OutputFolder, os.path.basename(Path_River_Polyline)
        )
        Selectfeatureattributes(
            processing,
            Input=Path_River_Polyline,
            Output=Outputfilename_cat_riv,
            Attri_NM="SubId",
            Values=HydroBasins_All,
        )

    finalcat_csv = Outputfilename_cat[:-3] + "dbf"
    finalcat_info = Dbf5(finalcat_csv)
    finalcat_info = finalcat_info.to_dataframe().drop_duplicates("SubId", keep="first")

    #### extract lakes

    Connect_Lake_info = finalcat_info.loc[finalcat_info["Lake_Cat"] == 1]
    Connect_Lakeids = np.unique(Connect_Lake_info["HyLakeId"].values)
    Connect_Lakeids = Connect_Lakeids[Connect_Lakeids > 0]

    NConnect_Lake_info = finalcat_info.loc[finalcat_info["Lake_Cat"] == 2]
    NonCL_Lakeids = np.unique(NConnect_Lake_info["HyLakeId"].values)
    NonCL_Lakeids = NonCL_Lakeids[NonCL_Lakeids > 0]


    if len(Connect_Lakeids) > 0 and Path_Con_Lake_ply != "#":
        Selectfeatureattributes(
            processing,
            Input=Path_Con_Lake_ply,
            Output=os.path.join(OutputFolder, os.path.basename(Path_Con_Lake_ply)),
            Attri_NM="Hylak_id",
            Values=Connect_Lakeids,
        )
    if len(NonCL_Lakeids) > 0 and Path_NonCon_Lake_ply != "#":
        Selectfeatureattributes(
            processing,
            Input=Path_NonCon_Lake_ply,
            Output=os.path.join(OutputFolder, os.path.basename(Path_NonCon_Lake_ply)),
            Attri_NM="Hylak_id",
            Values=NonCL_Lakeids,
        )

    Gauge_info = finalcat_info.loc[finalcat_info[Gauge_col_Name] > 0]
    Gauge_NMs = np.unique(Gauge_info["Obs_NM"].values)

    if len(Gauge_NMs) > 0 and Path_obs_gauge_point != "#":
        Selectfeatureattributes(
            processing,
            Input=Path_obs_gauge_point,
            Output=os.path.join(OutputFolder, os.path.basename(Path_obs_gauge_point)),
            Attri_NM="Obs_NM",
            Values=Gauge_NMs,
            Is_str = True,
        )
    
    Qgs.exit()
