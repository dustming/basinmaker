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
hyshdply = sys.argv[1]
hyshddem = sys.argv[2]
WidDep = sys.argv[3]
VolThreshold = int(sys.argv[4])
Lakefile = sys.argv[5]
obspoint = sys.argv[6]
OutputFolder = sys.argv[7] + "/"
Landuse = sys.argv[8]
Landuseinfo = sys.argv[9]

arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(4326) ### WGS84

if not os.path.exists(OutputFolder):
    os.makedirs(OutputFolder)

outFill = Fill(hyshddem)
outFill.save(OutputFolder + "dem")

arcpy.Copy_management(hyshdply, OutputFolder + "HyMask.shp")


cellSize = arcpy.Describe(OutputFolder + "dem").children[0].meanCellHeight

outFlowDirection = FlowDirection(outFill, "NORMAL")
outFlowDirection.save(OutputFolder + "dir")

outFlowAccumulation = FlowAccumulation(OutputFolder + "dir")
outFlowAccumulation.save(OutputFolder + "acc")


arcpy.RasterToASCII_conversion(OutputFolder + "dir", OutputFolder + "dir.asc")

######################################################3
arcpy.AddMessage("###################################3"+str(VolThreshold))


if VolThreshold >= 0:
    arcpy.AddMessage("#####################33")
    arcpy.Clip_analysis(Lakefile, OutputFolder +"HyMask.shp", OutputFolder +"HyLake1.shp", "")
    copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"HyLake.prj")
    where_clause = '"Lake_area"> '+ str(VolThreshold)
    arcpy.Select_analysis(OutputFolder +"HyLake1.shp", OutputFolder +"HyLake.shp", where_clause)
    copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"HyLake.prj")
else:
    arcpy.Clip_analysis(Lakefile, OutputFolder +"HyMask.shp", OutputFolder +"HyLake1.shp", "")
    copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"HyLake.prj")
    where_clause = '"Lake_area"> '+ str(1000000000000000)
    arcpy.Select_analysis(OutputFolder +"HyLake1.shp", OutputFolder +"HyLake.shp", where_clause)
    copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"HyLake.prj")
    tdir = np.loadtxt(OutputFolder+ 'dir.asc',dtype = 'i4',skiprows = 6)
    tlake = copy.copy(tdir)
    tlake[:,:] = -9999
    writeraster(OutputFolder + "hylake.asc",tlake,OutputFolder + 'dir')
###################3
####### Set envroment variable to dir

arcpy.env.XYTolerance = cellSize
arcpy.arcpy.env.cellSize = cellSize
arcpy.env.extent = arcpy.Describe(OutputFolder + "dir").extent
arcpy.env.snapRaster = OutputFolder + "dir"
#########################################################

if VolThreshold >= 0:
    arcpy.PolygonToRaster_conversion(OutputFolder +"HyLake.shp", "Hylak_id", OutputFolder + "hylake",
                                 "MAXIMUM_COMBINED_AREA","Hylak_id", cellSize)
    arcpy.RasterToASCII_conversion(OutputFolder + "hylake", OutputFolder + "hylake.asc")
##### dem


arcpy.RasterToASCII_conversion(OutputFolder + "dem", OutputFolder + "dem.asc")
##### acc
arcpy.RasterToASCII_conversion(OutputFolder + "acc", OutputFolder + "acc.asc")
#######land use

arcpy.AddMessage(Landuse)
if Landuse != "#":
    outExtractByMask = ExtractByMask(Landuse, OutputFolder +"dem")
    outExtractByMask.save(OutputFolder + "landuse_1")
    
    arcpy.Resample_management(OutputFolder + "landuse_1", OutputFolder + "landuse", cellSize, "NEAREST")

    arcpy.RasterToASCII_conversion(OutputFolder + "landuse", OutputFolder + "landuse.asc")
    copyfile(Landuseinfo, OutputFolder + "landuseinfo.csv")
else:
    tdir = np.loadtxt(OutputFolder+ 'dir.asc',dtype = 'i4',skiprows = 6)
    tlake = copy.copy(tdir)
    tlake[:,:] = -9999
    writeraster(OutputFolder + "landuse.asc",tlake,OutputFolder + 'dir')
    kk = pd.DataFrame(columns=['RasterV', 'MannV'],index=range(0,1))
    kk.loc[0,'RasterV'] = -9999
    kk.loc[0,'MannV'] = 0.035
    kk.to_csv(OutputFolder + "landuseinfo.csv",sep=",")
######################################################################################
#######hydrobasin
arcpy.PolygonToRaster_conversion(OutputFolder +"HyMask.shp", "FID", OutputFolder + "hybasinfid",
                                 "CELL_CENTER","NONE", cellSize)
arcpy.RasterToASCII_conversion(OutputFolder + "hybasinfid", OutputFolder + "hybasinfid.asc")

##### width and depth
arcpy.Clip_analysis(WidDep, OutputFolder +"HyMask.shp", OutputFolder + "WidDep.shp")
arcpy.PolylineToRaster_conversion(OutputFolder + "WidDep.shp", "WIDTH", OutputFolder + "width",
                                  "MAXIMUM_LENGTH", "NONE", cellSize)
arcpy.PolylineToRaster_conversion(OutputFolder + "WidDep.shp", "DEPTH", OutputFolder + "depth",
                                  "MAXIMUM_LENGTH", "NONE", cellSize)
arcpy.PolylineToRaster_conversion(OutputFolder + "WidDep.shp", "Q_Mean", OutputFolder + "Q_Mean",
                                  "MAXIMUM_LENGTH", "NONE", cellSize)
arcpy.PolylineToRaster_conversion(OutputFolder + "WidDep.shp", "Shape_Leng", OutputFolder + "WD_Len",
                                  "MAXIMUM_LENGTH", "NONE", cellSize)
arcpy.RasterToASCII_conversion( OutputFolder + "depth", OutputFolder + "depth.asc")
arcpy.RasterToASCII_conversion( OutputFolder + "width", OutputFolder + "width.asc")
arcpy.RasterToASCII_conversion( OutputFolder + "Q_Mean", OutputFolder + "Q_Mean.asc")
arcpy.RasterToASCII_conversion( OutputFolder + "WD_Len", OutputFolder + "WD_Len.asc")
########

#########prepare obspoints
if obspoint != "#":
    arcpy.PointToRaster_conversion(obspoint, "FID",
                                OutputFolder + "obs", "MAXIMUM", "", cellSize)
    arcpy.RasterToASCII_conversion( OutputFolder + "obs", OutputFolder + "obs.asc")
else:
    tdir = np.loadtxt(OutputFolder+ 'dir.asc',dtype = 'i4',skiprows = 6)
    tlake = copy.copy(tdir)
    tlake[:,:] = -9999
    writeraster(OutputFolder + "obs.asc",tlake,OutputFolder + 'dir')
#######converting dbf to csv files
dbftocsv( OutputFolder +"HyLake.dbf",OutputFolder +"lakeinfo.csv")
dbftocsv( OutputFolder +"HyMask.dbf",OutputFolder +"hybinfo.csv")
############
