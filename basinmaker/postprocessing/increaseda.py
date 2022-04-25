from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import pandas as pd
import numpy as np
import shutil
import tempfile


def simplify_routing_structure_by_drainage_area_qgis(
    Routing_Product_Folder= '#',
    Area_Min=-1,
    OutputFolder="#",
    qgis_prefix_path="#",
):
    """Simplify the routing product by drainage area

    Function that used to simplify the routing product by
    using user provided minimum subbasin drainage area.
    The input catchment polygons are routing product before
    merging for lakes. It is provided with routing product.
    The result is the simplified catchment polygons. But
    result from this fuction still not merging catchment
    covering by the same lake. Thus, The result generated
    from this tools need further processed by
    Define_Final_Catchment, or can be further processed by
    SelectLakes

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
    Area_Min                       : float
        The minimum drainage area of each catchment in km2
    OutputFolder                   : string
        Folder name that stores generated simplified routing product

    Notes
    -------
    This function has no return values, instead will generate following
    files. The output tpye will be the same as inputs, but the routing
    network will be simplified by increase subbasin size, reduce
    number of subbasins and number of river segments.

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
    #### generate river catchments based on minmum area.
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

    sub_colnm = "SubId"
    down_colnm = "DowSubId"
    DA_colnm = "DrainArea"
    SegID_colnm = "Seg_ID"

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

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



    # overall procedure,
    # 1. first get product attribute table
    # 2. determine which features needs to be merged together to increage
    #    drainage area of the sub basin, for example sub a, b c needs to be merged
    # 3. in the attribute table, change sub a b c 's content to a, assuming sub b and c drainge to a.
    # 4. copy modified attribute table to shafiles
    # 5. dissolve based on subid amd finished
    # overall procedure,

    # read attribute table, and
    Path_final_rviply = Path_final_riv_ply
    Path_final_riv = Path_final_riv
    Path_Conl_ply = Path_Con_Lake_ply
    Path_Non_ConL_ply = Path_NonCon_Lake_ply

    finalriv_info = Dbf_To_Dataframe(Path_final_rviply).drop_duplicates(
        sub_colnm, keep="first"
    )

    if Path_Conl_ply != '#':
        Conn_Lakes_ply = Dbf_To_Dataframe(Path_Conl_ply).drop_duplicates(
            "Hylak_id", keep="first"
        )
    else: 
        Conn_Lakes_ply = pd.DataFrame(np.full((10,1),-9999),columns=["Hylak_id"])
    # change attribute table
    (
        mapoldnew_info,
        Selected_riv_ids,
        Connected_Lake_Mainriv,
        Old_Non_Connect_LakeIds,
        Conn_To_NonConlakeids,
    ) = Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Increase_DA(
        finalriv_info, Conn_Lakes_ply, Area_Min
    )

    # remove river polyline
    Path_Temp_final_rviply = os.path.join(tempfolder, "temp1_finalriv_ply.shp")
    Path_Temp_final_rvi = os.path.join(tempfolder, "temp1_finalriv.shp")
    Selectfeatureattributes(
        processing,
        Input=Path_final_riv,
        Output=Path_Temp_final_rvi,
        Attri_NM="SubId",
        Values=Selected_riv_ids,
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_rviply,
        FIELD=["SubId"],
        OUTPUT=Path_Temp_final_rviply,
    )

    # update topology of new attribute table
    UpdateTopology(mapoldnew_info)
    mapoldnew_info = update_non_connected_catchment_info(mapoldnew_info)

    # create output folder
    outputfolder_subid = OutputFolder
    if not os.path.exists(outputfolder_subid):
        os.makedirs(outputfolder_subid)
        
    all_subids = finalriv_info['SubId'].values
    
    copy_data_and_dissolve(all_subids,tempfolder,processing,Path_Temp_final_rviply,Path_Temp_final_rvi,
        mapoldnew_info,COLUMN_NAMES_CONSTANT_CLEAN,OutputFolder,Path_Catchment_Polygon,context,
        Path_final_rviply,Path_final_riv)
        
        
    # # copy new attribute table to subbasin polyline and polygon
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



    # export lake polygons
    if Path_Conl_ply != '#' and len(Connected_Lake_Mainriv) > 0:
        Selectfeatureattributes(
            processing,
            Input=Path_Conl_ply,
            Output=os.path.join(outputfolder_subid, os.path.basename(Path_Conl_ply)),
            Attri_NM="Hylak_id",
            Values=Connected_Lake_Mainriv,
        )
    if Path_Non_ConL_ply != '#' and len(Old_Non_Connect_LakeIds) >0:
        Selectfeatureattributes(
            processing,
            Input=Path_Non_ConL_ply,
            Output=os.path.join(outputfolder_subid, os.path.basename(Path_Non_ConL_ply)),
            Attri_NM="Hylak_id",
            Values=Old_Non_Connect_LakeIds,
        )

    # Copy connected lakes that are transfered into non-connected
    # lake to non connected lake polygon

    if len(Conn_To_NonConlakeids) > 0 and Path_Conl_ply !='#':
        if Path_Non_ConL_ply != '#':
            Copyfeature_to_another_shp_by_attribute(
                Source_shp=Path_Conl_ply,
                Target_shp=os.path.join(
                    outputfolder_subid, os.path.basename(Path_Non_ConL_ply)
                ),
                Col_NM="Hylak_id",
                Values=Conn_To_NonConlakeids,
                Attributes=Conn_Lakes_ply,
            )
        else:
            if len(os.path.basename(Path_Conl_ply).split('_')) == 4:
                outlake_name = 'sl_non_connected_lake_' + os.path.basename(Path_Conl_ply).split('_')[3]
            else:
                outlake_name = 'sl_non_connected_lake.shp'
            Selectfeatureattributes(
                processing,
                Input=Path_Conl_ply,
                Output=os.path.join(outputfolder_subid, outlake_name),
                Attri_NM="Hylak_id",
                Values=Conn_To_NonConlakeids,
            )
                    
    # # dissolve subbasin polygon based on new subbasin id
    # Path_out_final_rviply = os.path.join(
    #     outputfolder_subid, os.path.basename(Path_final_riv_ply)
    # )
    # Path_out_final_rvi = os.path.join(
    #     outputfolder_subid, os.path.basename(Path_final_riv)
    # )
    # 
    # qgis_vector_dissolve(
    #     processing,
    #     context,
    #     INPUT=Path_Temp_final_rviply,
    #     FIELD=["SubId"],
    #     OUTPUT=Path_out_final_rviply,
    # )
    # qgis_vector_dissolve(
    #     processing,
    #     context,
    #     INPUT=Path_Temp_final_rvi,
    #     FIELD=["SubId"],
    #     OUTPUT=Path_out_final_rvi,
    # )
    # 
    # # clean attribute table and done
    # Clean_Attribute_Name(Path_out_final_rviply, COLUMN_NAMES_CONSTANT_CLEAN)
    # Clean_Attribute_Name(Path_out_final_rvi, COLUMN_NAMES_CONSTANT_CLEAN)

    Qgs.exit()
