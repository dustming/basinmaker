from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import pandas as pd
import numpy as np
import tempfile


def combine_catchments_covered_by_the_same_lake_qgis(
    # OutputFolder, 
    # Path_final_rivply="#", 
    # Path_final_riv="#", 
    Routing_Product_Folder = '#',
    qgis_prefix_path="#"
):
    """Define final lake river routing structure

    Generate the final lake river routing structure by merging subbasin
    polygons that are covered by the same lake.
    The input are the catchment polygons and river segements
    before merging for lakes. The input files can be output of
    any of following functions:
    SelectLakes, Select_Routing_product_based_SubId,
    Customize_Routing_Topology,RoutingNetworkTopologyUpdateToolset_riv
    The result is the final catchment polygon that ready to be used for
    hydrological modeling

    Parameters
    ----------
    OutputFolder                   : string
        Folder name that stores generated extracted routing product
    Path_final_riv_ply             : string
        Path to the catchment polygon which is the routing product
        before merging lakes catchments and need to be processed before
        used. It is the input for simplify the routing product based
        on lake area or drianage area.
        routing product and can be directly used.
    Path_final_riv                 : string
        Path to the river polyline which is the routing product
        before merging lakes catchments and need to be processed before
        used. It is the input for simplify the routing product based
        on lake area or drianage area.

    Notes
    -------
    This function has no return values, instead will generate following
    files. They are catchment polygons and river polylines that can be
    used for hydrological modeling.
    os.path.join(OutputFolder,'finalcat_info.shp')
    os.path.join(OutputFolder,'finalcat_info_riv.shp')

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


    OutputFolder = Routing_Product_Folder

    sub_colnm = "SubId"
    Path_final_rviply = Path_Catchment_Polygon
    Path_final_riv = Path_River_Polyline

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)
    tempfolder = os.path.join(
        tempfile.gettempdir(), "basinmaker_" + str(np.random.randint(1, 10000 + 1))
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)

    ### create a copy of shapfiles in temp folder
    Path_Temp_final_rviply = os.path.join(
        tempfolder,
        "temp_finalriv_ply" + str(np.random.randint(1, 10000 + 1)) + ".shp",
    )
    Path_Temp_final_rvi = os.path.join(
        tempfolder,
        "temp_finalriv" + str(np.random.randint(1, 10000 + 1)) + ".shp",
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_rviply,
        FIELD=["SubId"],
        OUTPUT=Path_Temp_final_rviply,
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_final_riv,
        FIELD=["SubId"],
        OUTPUT=Path_Temp_final_rvi,
    )

    ### read riv ply info
    ### read attribute table
    finalrivply_info = Dbf_To_Dataframe(Path_Temp_final_rviply).drop_duplicates(
        "SubId", keep="first"
    )
    # change attribute table for lake covered catchments,
    mapoldnew_info = change_attribute_values_for_catchments_covered_by_same_lake(
        finalrivply_info
    )
    
    mapoldnew_info.loc[mapoldnew_info['Lake_Cat'] > 0,'RivLength'] = 0
    # update topology for new attribute table
    update_topology(mapoldnew_info, UpdateStreamorder=-1)

    # copy new attribute table to shpfile
    Copy_Pddataframe_to_shpfile(
        Path_Temp_final_rviply,
        mapoldnew_info,
        link_col_nm_shp="SubId",
        link_col_nm_df="Old_SubId",
        UpdateColNM=["#"],
    )
    Copy_Pddataframe_to_shpfile(
        Path_Temp_final_rvi,
        mapoldnew_info,
        link_col_nm_shp="SubId",
        link_col_nm_df="Old_SubId",
        UpdateColNM=["#"],
    )

    # dissolve shpfile based on new subid
    if len(os.path.basename(Path_Catchment_Polygon).split('_')) == 5:
        Path_final_rviply = os.path.join(OutputFolder, "finalcat_info_"+os.path.basename(Path_Catchment_Polygon).split('_')[4])
        Path_final_rvi = os.path.join(OutputFolder, "finalcat_info_"+os.path.basename(Path_Catchment_Polygon).split('_')[4])
    else:
        Path_final_rviply = os.path.join(OutputFolder, "finalcat_info.shp")
        Path_final_rvi = os.path.join(OutputFolder, "finalcat_info.shp")
                
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_Temp_final_rvi,
        FIELD=["SubId"],
        OUTPUT=Path_final_rvi,
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=Path_Temp_final_rviply,
        FIELD=["SubId"],
        OUTPUT=Path_final_rviply,
    )

    # clean attribute table of shpfile
    Clean_Attribute_Name(Path_final_rvi, COLUMN_NAMES_CONSTANT)
    Clean_Attribute_Name(Path_final_rviply, COLUMN_NAMES_CONSTANT)

    # add centroid to new drived polygons
    Add_centroid_to_feature(Path_final_rviply, "centroid_x", "centroid_y")
