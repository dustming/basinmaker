from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import pandas as pd
import numpy as np
import shutil
import tempfile


def simplify_routing_structure_by_filter_lakes_qgis(
    Routing_Product_Folder = '#',
    Thres_Area_Conn_Lakes=-1,
    Thres_Area_Non_Conn_Lakes=-1,
    Selected_Lake_List_in=[],
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):
    """Simplify the routing product by lake area

    Function that used to simplify the routing product by user
    provided lake area thresthold.
    The input catchment polygons is the routing product before
    merging for lakes. It is provided with the routing product.
    The result is the simplified catchment polygons. But
    result from this fuction still not merging catchment
    covering by the same lake. Thus, The result generated
    from this tools need further processed by
    Define_Final_Catchment, or can be further processed by
    Customize_Routing_Topology

    Parameters
    ----------

    Path_final_riv_ply             : string
        Path to the catchment polygon which is the routing product
        before merging lakes catchments and need to be processed before
        used. It is the input for simplify the routing product based
        on lake area or drianage area.
    Path_final_riv                 : string
        Path to the river polyline which is the routing product
        before merging lakes catchments and need to be processed before
        used. It is the input for simplify the routing product based
        on lake area or drianage area.
    Path_Con_Lake_ply              : string
        Path to a connected lake polygon. Connected lakes are lakes that
        are connected by Path_final_riv.
    Path_NonCon_Lake_ply           : string
        Path to a non connected lake polygon. Connected lakes are lakes
        that are not connected by Path_final_riv.
    Thres_Area_Conn_Lakes          : float (optional)
        It is the lake area threshold for connated lakes, in km2
    Thres_Area_Non_Conn_Lakes      : float (optional)
        It is the lake area threshold for non connated lakes, in km2
    Selection_Method               : string
        It is a string indicate lake selection methods
        "ByArea" means lake in the routing product will be selected based
        on two lake area thresthold Thres_Area_Conn_Lakes and
        Thres_Area_Non_Conn_Lakes
        "ByLakelist" means lake in the routing product will be selected
        based on user provided hydrolake id, in Selected_Lake_List_in
    Selected_Lake_List_in          : list
        A list of lake ids that will be keeped in the routing product.
        Lakes not in the list will be removed from routing product.
    OutputFolder                   : string
        Folder name that stores generated extracted routing product

    Notes
    -------
    This function has no return values, instead will generate following
    files. The output tpye will be the same as inputs, but the routing
    network will be simplified by removing lakes.

    os.path.join(OutputFolder,os.path.basename(Path_final_riv_ply))
    os.path.join(OutputFolder,os.path.basename(Path_final_riv))
    os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply))
    os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply))

    Returns:
    -------
    None

    Examples
    -------

    """
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

    
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

    tempfolder = os.path.join(
        tempfile.gettempdir(),
        "basinmaker_sllake_" + str(np.random.randint(1, 10000 + 1)),
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

    Path_final_riv_ply = Path_Catchment_Polygon
    Path_final_riv = Path_River_Polyline
    
    ## copy obs_gauges to output folder 
    for file in os.listdir(Routing_Product_Folder):
        if 'obs_gauges' in file:
            shutil.copy(os.path.join(Routing_Product_Folder, file), os.path.join(OutputFolder, file))

    
    ### read attribute table
    finalcat_info = Dbf_To_Dataframe(Path_final_riv_ply)
    finalcat_info = finalcat_info.fillna(-9999)

    ### Obtain selected lake's attribute info
    (
        Selected_Non_ConnLakes,
        Selected_ConnLakes,
        Un_Selected_ConnLakes_info,
        Un_Selected_Non_ConnL_info,
    ) = Return_Selected_Lakes_Attribute_Table_And_Id(
        finalcat_info,
        Thres_Area_Conn_Lakes*1000*1000,
        Thres_Area_Non_Conn_Lakes*1000*1000,
        Selected_Lake_List_in,
    )

    ### Extract lake polygons
    if len(Selected_Non_ConnLakes) > 0:
        Selectfeatureattributes(
            processing,
            Input=Path_NonCon_Lake_ply,
            Output=os.path.join(OutputFolder, os.path.basename(Path_NonCon_Lake_ply)),
            Attri_NM="Hylak_id",
            Values=Selected_Non_ConnLakes,
        )
    if len(Selected_ConnLakes) > 0:
        Selectfeatureattributes(
            processing,
            Input=Path_Con_Lake_ply,
            Output=os.path.join(OutputFolder, os.path.basename(Path_Con_Lake_ply)),
            Attri_NM="Hylak_id",
            Values=Selected_ConnLakes,
        )

    ### create a copy of shapfiles in temp folder
    Path_Temp_final_rviply = os.path.join(
        tempfolder, "temp_finalriv_ply_selectlake.shp"
    )
    Path_Temp_final_rvi = os.path.join(tempfolder, "temp_finalriv_selectlake.shp")
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_riv,
        FIELD=["SubId"],
        OUTPUT=Path_Temp_final_rvi,
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_riv_ply,
        FIELD=["SubId"],
        OUTPUT=Path_Temp_final_rviply,
    )

    # change lake related attribute for un selected connected lake
    # catchment to -1.2345
    Remove_Unselected_Lake_Attribute_In_Finalcatinfo(
        Path_Temp_final_rviply, Selected_ConnLakes
    )
    # remove lake attributes
    Remove_Unselected_Lake_Attribute_In_Finalcatinfo(
        Path_Temp_final_rvi, Selected_ConnLakes
    )

    # Modify attribute table to merge un selected lake catchment if needed
    finalcat_info_temp = Dbf_To_Dataframe(Path_Temp_final_rviply)
    # change attributes for catchment due to remove of connected lakes
    mapoldnew_info = (
        Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Remove_CL(
            finalcat_info_temp, Un_Selected_ConnLakes_info
        )
    )
    # change attribute for catchment due to remove of non connected lakes
    mapoldnew_info = (
        Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Remove_NCL(
            mapoldnew_info, finalcat_info_temp, Un_Selected_Non_ConnL_info
        )
    )
    
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0,'HyLakeId'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0,'Lake_Cat'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0,'LakeVol'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0,'LakeDepth'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0,'LakeArea'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0,'Laketype'] = 0
    
    # update topology for new attribute table
    UpdateTopology(mapoldnew_info, UpdateStreamorder=-1)
    mapoldnew_info = update_non_connected_catchment_info(mapoldnew_info)


    all_subids = finalcat_info_temp['SubId'].values
    
    copy_data_and_dissolve(all_subids,tempfolder,processing,Path_Temp_final_rviply,Path_Temp_final_rvi,
        mapoldnew_info,COLUMN_NAMES_CONSTANT_CLEAN,OutputFolder,Path_Catchment_Polygon,context,
        Path_final_riv_ply,Path_final_riv)
        
        
    # # copy new attribute table to shpfiles
    # Copy_Pddataframe_to_shpfile(
    #     Path_Temp_final_rviply,
    #     mapoldnew_info,
    #     link_col_nm_shp="SubId",
    #     link_col_nm_df="Old_SubId",
    #     UpdateColNM=["#"],
    # )
    # Copy_Pddataframe_to_shpfile(
    #     Path_Temp_final_rvi,
    #     mapoldnew_info,
    #     link_col_nm_shp="SubId",
    #     link_col_nm_df="Old_SubId",
    #     UpdateColNM=["#"],
    # )
    # 
    # # disslove line and polygon based on new subid
    # qgis_vector_dissolve(
    #     processing,
    #     context,
    #     INPUT=Path_Temp_final_rvi,
    #     FIELD=["SubId"],
    #     OUTPUT=os.path.join(OutputFolder, os.path.basename(Path_final_riv)),
    # )
    # qgis_vector_dissolve(
    #     processing,
    #     context,
    #     INPUT=Path_Temp_final_rviply,
    #     FIELD=["SubId"],
    #     OUTPUT=os.path.join(OutputFolder, os.path.basename(Path_final_riv_ply)),
    # )
    # 
    # # clean attribute table
    # Clean_Attribute_Name(
    #     os.path.join(OutputFolder, os.path.basename(Path_final_riv)),
    #     COLUMN_NAMES_CONSTANT,
    # )
    # Clean_Attribute_Name(
    #     os.path.join(OutputFolder, os.path.basename(Path_final_riv_ply)),
    #     COLUMN_NAMES_CONSTANT,
    # )

    return
