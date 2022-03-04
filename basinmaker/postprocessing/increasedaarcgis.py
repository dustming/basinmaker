import numpy as np
import sys
import os
import csv
import tempfile 
import copy 
import pandas as pd
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import arcpy
from arcpy import env
from arcpy.sa import *
import shutil


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

def simplify_routing_structure_by_drainage_area_arcgis(
    Routing_Product_Folder= '#',
    Area_Min=-1,
    OutputFolder="#",
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

    finalriv_infoply = pd.DataFrame.spatial.from_featureclass(Path_final_rviply)
    finalriv_inforiv = pd.DataFrame.spatial.from_featureclass(Path_final_riv)
    if Path_Conl_ply != '#':
        Conn_Lakes_ply = pd.DataFrame.spatial.from_featureclass(Path_Conl_ply)
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
        finalriv_infoply, Conn_Lakes_ply, Area_Min
    )

    finalriv_inforiv_main = finalriv_inforiv.loc[finalriv_inforiv['SubId'].isin(Selected_riv_ids)]
    finalriv_inforiv_main.spatial.to_featureclass(location=os.path.join(tempfolder,'selected_riv.shp'),overwrite=True,sanitize_columns=False)
    
    UpdateTopology(mapoldnew_info)
    mapoldnew_info = update_non_connected_catchment_info(mapoldnew_info)

    save_modified_attributes_to_outputs(
        mapoldnew_info=mapoldnew_info,
        tempfolder=tempfolder,
        OutputFolder=OutputFolder,
        cat_name=os.path.basename(Path_final_riv_ply),
        riv_name = os.path.basename(Path_final_riv),
        Path_final_riv = os.path.join(tempfolder,'selected_riv.shp'),
    )
   
    # export lakes 
    if Path_Conl_ply == '#':
        Conn_Lakes_ply_select = []
        Conn_Lakes_ply_not_select = []
    else:     
        Conn_Lakes_ply = pd.DataFrame.spatial.from_featureclass(Path_Conl_ply)    
        lake_mask = Conn_Lakes_ply['Hylak_id'].isin(Connected_Lake_Mainriv)
        Conn_Lakes_ply_select = Conn_Lakes_ply.loc[lake_mask].copy()
        Conn_Lakes_ply_not_select = Conn_Lakes_ply.loc[np.logical_not(lake_mask)].copy()

    # export lake polygons
    # export connected lake polygon
    if len(Conn_Lakes_ply_select) > 0:
        Conn_Lakes_ply_select.spatial.to_featureclass(location=os.path.join(OutputFolder,os.path.basename(Path_Conl_ply)),overwrite=True,sanitize_columns=False)
    
    # export non connected polygon 
    if len(Conn_Lakes_ply_not_select) >0 and Path_NonCon_Lake_ply != '#':
        non_conn_Lakes_ply = pd.DataFrame.spatial.from_featureclass(Path_NonCon_Lake_ply)  
        new_non_connected_lake = pd.concat([non_conn_Lakes_ply, Conn_Lakes_ply_not_select], ignore_index=True)
        new_non_connected_lake.spatial.to_featureclass(location=os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)),overwrite=True,sanitize_columns=False)
    if len(Conn_Lakes_ply_not_select) <= 0 and Path_NonCon_Lake_ply != '#':
        non_conn_Lakes_ply = pd.DataFrame.spatial.from_featureclass(Path_NonCon_Lake_ply)
        non_conn_Lakes_ply.spatial.to_featureclass(location=os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply)),overwrite=True,sanitize_columns=False)
    if len(Conn_Lakes_ply_not_select) > 0 and Path_NonCon_Lake_ply == '#':
       
        if len(os.path.basename(Path_Conl_ply).split('_')) == 4:
            outlake_name = 'sl_non_connected_lake_' + os.path.basename(Path_Conl_ply).split('_')[3]
        else:
            outlake_name = 'sl_non_connected_lake.shp'
            
        Conn_Lakes_ply_not_select.spatial.to_featureclass(location=os.path.join(OutputFolder,outlake_name),overwrite=True,sanitize_columns=False)
        
    return 
