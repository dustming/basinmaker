import numpy as np
from scipy.optimize import curve_fit
import arcpy
from arcpy import env
from arcpy.sa import *
import copy
import sys
import shutil
import os
import csv
from simpledbf import Dbf5
import pandas as pd
from shutil import copyfile
import tempfile
from basinmaker.func.pdtable import *

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
##### Readed inputs

def Area_Weighted_Mapping_Between_Two_Polygons_Arcgis(Target_Ply_Path,Mapping_Ply_Path,Col_NM="HRU_ID",Output_Folder="#"):
    """Generate Grid polygon from NetCDF file
    Function that used to generate grid polygon from a NetCDF file
    Parameters
    ----------
    Target_Ply_Path                       : string
        It is the path of one inputs HRU polygon file
    Mapping_Ply_Path                      : string
        It is the path of one inputs grid polygon file
    Output_Folder                     : string
        It is the path to a folder to save output polygon shpfiles
    Notes
    -------
    Overlay_Polygons.shp                 : Polygon shpfile (output)
       It is overlay of two input polygons shpfiles
    GriddedForcings2.txt                 : Text file (output)
       It is the polygon area weighted of each polygon in Mapping_Ply_Path
       to each polygon in Target_Ply_Path
    Returns:
    -------
       None
    Examples
    -------
    """
    if not os.path.exists(Output_Folder):
        os.makedirs(Output_Folder)
        
    tempfolder = os.path.join(
        tempfile.gettempdir(), "basinmaker_" + str(np.random.randint(1, 10000 + 1))
    )
    if not os.path.exists(tempfolder):
        os.makedirs(tempfolder)

    Focspre = arcpy.Describe(Mapping_Ply_Path).spatialReference

    arcpy.Project_management(Target_Ply_Path, os.path.join(tempfolder,"temp1.shp"),Focspre)
    arcpy.RepairGeometry_management(os.path.join(tempfolder,"temp1.shp"))
    
    arcpy.Identity_analysis(os.path.join(tempfolder,"temp1.shp"), Mapping_Ply_Path,os.path.join(tempfolder,"finalcat_Frocing.shp"))
    arcpy.Dissolve_management(os.path.join(tempfolder,"finalcat_Frocing.shp"), os.path.join(Output_Folder,"Overlay_Polygons.shp"),[Col_NM, "FGID"
    ,"Row","Col"], "", "", "")
    
    if Focspre.type == "Geographic":
        arcpy.env.CoordinateSystem = arcpy.SpatialReference(3573)####wgs84 - north pore canada
    
    arcpy.AddField_management(os.path.join(Output_Folder,"Overlay_Polygons.shp"),"s_area","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management(os.path.join(Output_Folder,"Overlay_Polygons.shp"),"s_area","!shape.area@squaremeters!","PYTHON_9.3","#")

    dbf1 = Dbf5(Mapping_Ply_Path[:-3] + "dbf")
    Forcinfo = dbf1.to_dataframe()
    Avafgid = Forcinfo["FGID"].values

    dbf2 = Dbf5(os.path.join(Output_Folder,"Overlay_Polygons.shp")[:-3] + "dbf")
    Mapforcing = dbf2.to_dataframe()
    Mapforcing["Map_FGID"] = Mapforcing["FGID"]
    Mapforcing["Map_Row"] = Mapforcing["Row"]
    Mapforcing["Map_Col"] = Mapforcing["Col"]        
    Mapforcing = Mapforcing.loc[Mapforcing[Col_NM] > 0]  ### remove
    Mapforcing = Mapforcing.loc[Mapforcing["Map_FGID"] > 0]  ### remove
    Mapforcing = Mapforcing.loc[Mapforcing["s_area"] > 0.000001]  ### remove 
    
    grid_weight_string = create_grid_weight_main(Mapforcing,Forcinfo)
    
    ####
    grid_weight_file_path = os.path.join(Output_Folder, "GriddedForcings2.txt")
        
    WriteStringToFile(grid_weight_string, grid_weight_file_path, "w")

