# Copyright 2018-2018 Ming Han - ming.han(at)uwaterloo.ca
#
# License
# This file is part of Ming Han's personal code library.
#
# Ming Han's personal code library is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ming Han's personal code library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public License
# along with Ming Han's personal code library.  If not, see <http://www.gnu.org/licenses/>.
#

def Generatencply(ncfile,thre,clusTol,WorkingFolder):
    dsin2 =  Dataset(ncfile,'r') # sample structure of in nc file converted from fst
    np.savetxt(WorkingFolder + "lon.csv", dsin2.variables['lon'][:,:], delimiter=",")
    np.savetxt(WorkingFolder + "lat.csv", dsin2.variables['lat'][:,:], delimiter=",")
    np.savetxt(WorkingFolder + "rlat.csv", dsin2.variables['rlat'][:], delimiter=",")
    np.savetxt(WorkingFolder + "rlon.csv", dsin2.variables['rlon'][:], delimiter=",")
#    np.savetxt(WorkingFolder + "FI_SFC.csv", dsin2.variables['FI_SFC'][0,:,:], delimiter=",")
    ncols = len(dsin2.variables['lon'][0,:])  ### from 0 to (ncols-1).
    nrows = len(dsin2.variables['lon'][:,0])
    latlonrow = np.full((nrows*ncols,4),-9999)
    for i in range(0,nrows):
        for j in range(0,ncols):
            k = i*ncols + j
            latlonrow[k,0] = dsin2.variables['lon'][i,j]
            latlonrow[k,1] = dsin2.variables['lat'][i,j]
            latlonrow[k,2] = i
            latlonrow[k,3] = j
    pdlatlonrow = pd.DataFrame(latlonrow,columns=['lon','lat','irow','icol'])
    pdlatlonrow.to_csv(WorkingFolder + "Gridcorr.csv",sep = ',',index = False)

    pt = arcpy.Point()
    ptGeoms = []
    for i in range(0,len(pdlatlonrow)):
        pt.X = pdlatlonrow.ix[i]['lon']
        pt.Y = pdlatlonrow.ix[i]['lat']
        ptGeoms.append(arcpy.PointGeometry(pt))
    arcpy.CopyFeatures_management(ptGeoms, WorkingFolder + "Gridncpoint.shp")
    arcpy.DefineProjection_management(WorkingFolder + "Gridncpoint.shp", arcpy.SpatialReference(4326))
    arcpy.AddXY_management(WorkingFolder + "Gridncpoint.shp")

    dbfile = WorkingFolder + "Gridncpoint.shp"
    inFeatures = dbfile
    fieldPrecision = 16
    field_scale = 12
    arcpy.AddField_management(dbfile, "FGID", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Row", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Col", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Gridlon", "Double", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Gridlat", "Double", fieldPrecision,field_scale,"", "", "NULLABLE","","")

    for i in range(0,nrows):
        for j in range(0,ncols):
            lon = dsin2.variables['lon'][i,j]
            lat = dsin2.variables['lat'][i,j]
            addgridinfo2point(dbfile,lon,lat,i,j,thre,ncols,nrows)

    features = []
    feature_info1 = np.full((nrows +1,ncols,2),-1)
    k = 0
    for irow in range(0,nrows):
        if irow == nrows -1:
            lon1 = dsin2.variables['lon'][irow,:]
            lon2 = dsin2.variables['lon'][irow,:]
            lat1 = dsin2.variables['lat'][irow,:]
            lat2 = dsin2.variables['lat'][irow,:]
        else:
            lon1 = dsin2.variables['lon'][irow,:]
            lon2 = dsin2.variables['lon'][irow + 1,:]
            lat1 = dsin2.variables['lat'][irow,:]
            lat2 = dsin2.variables['lat'][irow + 1,:]
        for i in range(0,ncols):
            feature_info1[k,i,0] = (lon1[i] + lon2[i])/2
            feature_info1[k,i,1] = (lat1[i] + lat2[i])/2
        k = k + 1
        if irow == 0:
            lon1 = dsin2.variables['lon'][irow,:]
            lon2 = dsin2.variables['lon'][irow ,:]
            lat1 = dsin2.variables['lat'][irow,:]
            lat2 = dsin2.variables['lat'][irow ,:]
            for i in range(0,ncols):
                feature_info1[k,i,0] = (lon1[i] + lon2[i])/2
                feature_info1[k,i,1] = (lat1[i] + lat2[i])/2
            k = k + 1

    for feature in feature_info1:
        features.append(
        arcpy.Polyline(arcpy.Array([arcpy.Point(*coords) for coords in feature])))

    feature_info2 = np.full((ncols+1,nrows,2),-1)
    k=0
    for icol in range(0,ncols):
        if icol == ncols - 1:
            lon1 = dsin2.variables['lon'][:,icol]
            lon2 = dsin2.variables['lon'][:,icol]
            lat1 = dsin2.variables['lat'][:,icol]
            lat2 = dsin2.variables['lat'][:,icol]
        else:
            lon1 = dsin2.variables['lon'][:,icol]
            lon2 = dsin2.variables['lon'][:,icol+1]
            lat1 = dsin2.variables['lat'][:,icol]
            lat2 = dsin2.variables['lat'][:,icol + 1]
        for i in range(0,nrows):
            feature_info2[k,i,0] = (lon1[i] + lon2[i])/2
            feature_info2[k,i,1] = (lat1[i] + lat2[i])/2
        k = k + 1
        if icol == 0:
            lon1 = dsin2.variables['lon'][:,icol]
            lon2 = dsin2.variables['lon'][:,icol]
            lat1 = dsin2.variables['lat'][:,icol]
            lat2 = dsin2.variables['lat'][:,icol]
            for i in range(0,nrows):
                feature_info2[k,i,0] = (lon1[i] + lon2[i])/2
                feature_info2[k,i,1] = (lat1[i] + lat2[i])/2
            k = k + 1
    for feature in feature_info2:
        features.append(
        arcpy.Polyline(arcpy.Array([arcpy.Point(*coords) for coords in feature])))
    dsin2.close()
    arcpy.CopyFeatures_management(features, WorkingFolder + "Gridline.shp")
    arcpy.DefineProjection_management(WorkingFolder + "Gridline.shp", arcpy.SpatialReference(4326))
    arcpy.FeatureToPolygon_management(WorkingFolder + "Gridline.shp", WorkingFolder + "Gridply.shp", clusTol, "NO_ATTRIBUTES", "")
    arcpy.SpatialJoin_analysis(WorkingFolder + "Gridply.shp", WorkingFolder + "Gridncpoint.shp", WorkingFolder + "Gridncply.shp", "JOIN_ONE_TO_ONE", "#", "#","INTERSECT","0.01")


def addgridinfo2point(dbfile,lon,lat,i,j,thre,ncols,nrows):
	rows = arcpy.UpdateCursor(dbfile)
	for row in rows:
		if abs (row.POINT_X - lon) < thre and abs(row.POINT_Y - lat) < thre:
			row.FGID = i*(ncols)  + j   # ncols = max(colnumers) + 1
			row.Row = i
			row.Col = j
			row.Gridlon = round(lon, 5)
			row.Gridlat = round(lat, 5)
			rows.updateRow(row)
	del row
	del rows

import datetime as dt  # Python standard library datetime  module
import numpy as np
import pandas as pd
from netCDF4 import Dataset, num2date, date2num
from datetime import datetime, timedelta
import os
import arcpy
from arcpy import env
from arcpy.sa import *
from simpledbf import Dbf5
from dbfpy import dbf
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

############### part1 create empty nc sample file with all possible variable
WorkingFolder = sys.argv[1]
ncfilename = sys.argv[2]
catply = sys.argv[3]
ascinput  = -1
############### Change the input in this section
#WorkingFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/RDRS'
#ncfilename = 'RDRS_CaPA24hr_forcings_final.nc'
#catply = 'finalcat.shp'
###############################################
clusTol = "0.01"
thre = 0.001
WorkingFolder = WorkingFolder + '/'
arcpy.env.workspace = WorkingFolder
os.chdir(WorkingFolder)
outfolder = WorkingFolder

ncfile = ncfilename
Generatencply(ncfile,thre,clusTol,WorkingFolder)
