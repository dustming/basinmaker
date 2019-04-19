def dbftocsv(filename,outname):
    if filename.endswith('.dbf'):
#        print "Converting %s to csv" % filename
        csv_fn = outname
        with open(csv_fn,'wb') as csvfile:
            in_db = dbf.Dbf(filename)
            out_csv = csv.writer(csvfile)
            names = []
            for field in in_db.header.fields:
                names.append(field.name)
            out_csv.writerow(names)
            for rec in in_db:
                out_csv.writerow(rec.fieldData)
            in_db.close()
#            print "Done..."
    else:
        print "Filename does not end with .dbf"

###################################################################3
def Defcat(out,outletid):
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((10000000,1),-99999999999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((10000000,1),-999999999999999)
        poshdid = 0
        for i in range(0,len(otsheds)):
            Shedid[psid] = otsheds[i]
            psid = psid + 1
            irow = np.argwhere(rout[:,1]==otsheds[i]).astype(int)
            for j in range(0,len(irow)):
                noutshd[poshdid] = rout[irow[j],0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd>=0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid>=0]
    return Shedid

def writeraster(w_filname,nraster,dataset):
    orvh = open(w_filname,"w")
    ncols = arcpy.GetRasterProperties_management(dataset, "COLUMNCOUNT")
    nrows = arcpy.GetRasterProperties_management(dataset, "ROWCOUNT")
    xllcorner = arcpy.GetRasterProperties_management(dataset, "LEFT")
    yllcorner = arcpy.GetRasterProperties_management(dataset, "BOTTOM")
    orvh.write("ncols      "+str(ncols) + "\n")
    orvh.write("nrows      "+ str(nrows) + "\n")
    orvh.write("xllcorner    "+str(xllcorner) + "\n")
    orvh.write("yllcorner    "+str(yllcorner) + "\n")
    orvh.write("cellsize     "+str(cellSize) + "\n")
    orvh.write("NODATA_value  -9999" + "\n")
    orvh.close()
    f_handle = open(w_filname, 'a')
    np.savetxt(f_handle,nraster,fmt='%i')
    f_handle.close()


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
import csv
from dbfpy import dbf
import pandas as pd
from shutil import copyfile
from simpledbf import Dbf5
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
##### Readed inputs
OutHyID = int(sys.argv[1])
hyshdply = sys.argv[2]
OutputFolder = sys.argv[3] + "/"
OutHyID2 = int(sys.argv[4])
tempinfo = Dbf5(hyshdply[:-3]+'dbf')#np.genfromtxt(hyinfocsv,delimiter=',')
hyshdinfo2 = tempinfo.to_dataframe()
hyshdinfo = hyshdinfo2[['SubId','DowSubId']].values

SptailRef = arcpy.Describe(hyshdply).spatialReference

arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84

if not os.path.exists(OutputFolder):
    os.makedirs(OutputFolder)
HydroBasins1 = Defcat(hyshdinfo,OutHyID)
if OutHyID2 > 0:
    HydroBasins2 = Defcat(hyshdinfo,OutHyID2)
#    arcpy.AddMessage(HydroBasins1)
#    arcpy.AddMessage(HydroBasins2)
    for i in range(len(HydroBasins2)):
        rows =np.argwhere(HydroBasins1 == HydroBasins2[i])
#        arcpy.AddMessage(rows)
        HydroBasins1 = np.delete(HydroBasins1, rows)
#    arcpy.AddMessage(HydroBasins1)
    HydroBasins = HydroBasins1
else:
    HydroBasins = HydroBasins1
#arcpy.AddMessage(HydroBasins)
out_feature_class = OutputFolder +"finalcat_info.shp"
where_clause = '"SubId" IN'+ " ("
for i in range(0,len(HydroBasins)):
    if i == 0:
        where_clause = where_clause + str(HydroBasins[i])
    else:
        where_clause = where_clause + "," + str(HydroBasins[i])
where_clause = where_clause + ")"
arcpy.Select_analysis(hyshdply, out_feature_class, where_clause)
##################
