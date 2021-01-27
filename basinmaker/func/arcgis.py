from arcgis import GIS
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import arcpy
from arcpy import env
from arcpy.sa import *
#####
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

def select_feature_by_attributes_arcgis(input,Attri_NM,Attri_v,output):
    where_clause = '"%s" IN' % (Attri_NM) 
    where_clause = where_clause + " (" 
    for i in range(0,len(Attri_v)):
        if i == 0:
            where_clause = where_clause + str(Attri_v[i])
        else:
            where_clause = where_clause + "," + str(Attri_v[i])
    where_clause = where_clause + ")"
    
    arcpy.Select_analysis(input, output, where_clause)
    return
##################

