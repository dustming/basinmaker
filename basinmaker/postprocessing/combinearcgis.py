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


import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from func.arcgis import *
from func.pdtable import *
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

def combine_catchments_covered_by_the_same_lake_arcgis(
    OutputFolder, Path_final_rivply="#", Path_final_riv="#"
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


    sub_colnm = "SubId"
    Path_final_rviply = Path_final_rivply
    Path_final_riv = Path_final_riv

    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)
        
    tempfolder = os.path.join(
        tempfile.gettempdir(), "basinmaker_" + str(np.random.randint(1, 10000 + 1))
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)
    arcpy.env.workspace = tempfolder

    ### create a copy of shapfiles in temp folder
    Path_Temp_final_rviply = os.path.join(OutputFolder,"temp_finalriv_ply" + str(np.random.randint(1, 10000 + 1)) + ".shp")
    
    Path_Temp_final_rvi = os.path.join(OutputFolder,"temp_finalriv_riv" + str(np.random.randint(1, 10000 + 1)) + ".shp")
    
    ### read riv ply info
    ### read attribute table
    finalrivply_info = pd.DataFrame.spatial.from_featureclass(Path_final_rivply)
    # change attribute table for lake covered catchments,
    finalrivply_info['SubId'] = finalrivply_info['SubId'].astype('int32')
    finalrivply_info['DowSubId'] = finalrivply_info['DowSubId'].astype('int32')
    finalrivply_info['HyLakeId'] = finalrivply_info['HyLakeId'].astype('int32')
#    finalrivply_info['DrainArea'] = finalrivply_info['DrainArea'].astype('float')

    mapoldnew_info = change_attribute_values_for_catchments_covered_by_same_lake(
        finalrivply_info
    )

    # update topology for new attribute table
    mapoldnew_info = update_topology(mapoldnew_info, UpdateStreamorder=-1)    
    
    mapoldnew_info['DowSubId'] = mapoldnew_info['DowSubId'].astype('int32')

    save_modified_attributes_to_outputs(
        mapoldnew_info=mapoldnew_info,
        tempfolder=tempfolder,
        OutputFolder=OutputFolder,
        cat_name='finalcat_info.shp',
        riv_name ='finalcat_info_riv.shp',
        Path_final_riv = Path_final_riv,
    )
    return 
