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

def simplify_routing_structure_by_filter_lakes_purepy(
    Routing_Product_Folder = '#',
    Thres_Area_Conn_Lakes=-1,
    Thres_Area_Non_Conn_Lakes=-1,
    Selected_Lake_List_in=[],
    OutputFolder="#",
    qgis_prefix_path="#",
    gis_platform="arcgis",
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
    finalcat_info = geopandas.read_file(Path_final_riv_ply)

    ### Obtain selected lake's attribute info
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
        sl_con_lakes = geopandas.read_file(Path_NonCon_Lake_ply)
        sl_con_lakes = sl_con_lakes.loc[sl_con_lakes['Hylak_id'].isin(Selected_Non_ConnLakes)]
        sl_con_lakes.to_file(os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)))

    if len(Selected_ConnLakes) > 0:
        sl_con_lakes = geopandas.read_file(Path_Con_Lake_ply)
        sl_con_lakes = sl_con_lakes.loc[sl_con_lakes['Hylak_id'].isin(Selected_ConnLakes)]
        sl_con_lakes.to_file(os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply)))

    print(" Obtain selected Lake IDs done")
 
    finalcat_ply = geopandas.read_file(Path_final_riv_ply)
    # change lake related attribute for un selected connected lake
    # catchment to -1.2345
    finalcat_ply = Remove_Unselected_Lake_Attribute_In_Finalcatinfo_purepy(
        finalcat_ply, Selected_ConnLakes
    )
    # remove lake attribute
        
    # Modify attribute table to merge un selected lake catchment if needed
    # change attributes for catchment due to remove of connected lakes
    mapoldnew_info = (
        Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Remove_CL(
            finalcat_ply, Un_Selected_ConnLakes_info
        )
    )
    # change attribute for catchment due to remove of non connected lakes
    mapoldnew_info = (
        Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Remove_NCL(
            mapoldnew_info, finalcat_ply, Un_Selected_Non_ConnL_info
        )
    )
    
    # update topology for new attribute table
    mapoldnew_info = UpdateTopology(mapoldnew_info, UpdateStreamorder=-1)
    mapoldnew_info = update_non_connected_catchment_info(mapoldnew_info)
    
    save_modified_attributes_to_outputs(
        mapoldnew_info=mapoldnew_info,
        tempfolder=tempfolder,
        OutputFolder=OutputFolder,
        cat_name=os.path.basename(Path_final_riv_ply),
        riv_name = os.path.basename(Path_final_riv),
        Path_final_riv = Path_final_riv,
    )

    return
