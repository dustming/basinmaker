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
def Maphru2forceply(forcingply,outfolder,catply,outFolderraven,Boundaryply,missrow,misscol):
    dbf1 = Dbf5(forcingply[:-3]+'dbf')
    Forcinfo = dbf1.to_dataframe()
    Focspre = arcpy.Describe(forcingply).spatialReference
    # run the tool
    if Boundaryply != "#":
        arcpy.Project_management(Boundaryply, outfolder+ "Boundary_freferen.shp",Focspre)
        arcpy.Identity_analysis(outfolder+ "Boundary_freferen.shp", forcingply,outfolder+ "Boundary_Frocing.shp")
        dbf3 = Dbf5(outfolder+ "Boundary_Frocing.dbf")
        BounForc = dbf3.to_dataframe()
        Avafgid = BounForc['FGID'].values
        Avafgid = np.unique(Avafgid)
    else:
        Avafgid = Forcinfo['FGID'].values
    arcpy.Project_management(catply, outfolder+ "finalcat_freferen.shp",Focspre)
    arcpy.Identity_analysis(outfolder+ "finalcat_freferen.shp", forcingply,outfolder+ "finalcat_Frocing.shp")
    arcpy.env.CoordinateSystem = arcpy.SpatialReference(3573)####wgs84 - north pore canada
    arcpy.AddField_management(outfolder +"finalcat_Frocing.shp","s_area","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management(outfolder +"finalcat_Frocing.shp","s_area","!shape.area@squaremeters!","PYTHON_9.3","#")
    dbf2 = Dbf5(outfolder+ "finalcat_Frocing.dbf")
    Mapforcing = dbf2.to_dataframe()
    catids = Mapforcing['SubId'].values
    catids = np.unique(catids)
    Lakeids = Mapforcing['HyLakeId'].values
    Lakeids = np.unique(Lakeids)
    Lakeids = Lakeids[Lakeids>0]
    ogridforc = open(outFolderraven+"GriddedForcings2.txt","w")
    ogridforc.write(":GridWeights" +"\n")
    ogridforc.write("   #      " +"\n")
    ogridforc.write("   # [# HRUs]"+"\n")
    sNhru = len(catids) + len(Lakeids)
    ogridforc.write("   :NumberHRUs       "+ str(sNhru) + "\n")
    sNcell = (max(Forcinfo['Row'].values)+1+missrow) * (max(Forcinfo['Col'].values)+1+misscol)
    ogridforc.write("   :NumberGridCells  "+str(sNcell)+"\n")
    ogridforc.write("   #            "+"\n")
    ogridforc.write("   # [HRU ID] [Cell #] [w_kl]"+"\n")
    maparray = np.full((1000000,5),-9)
    imap = 0
    maxcatid = max(catids)
#    arcpy.AddMessage(catids)
    for i in range(len(catids)):
#        print i
        catid = catids[i]
        cats = Mapforcing.loc[Mapforcing['SubId'] == catid]
        cats = cats[cats['FGID'].isin(Avafgid)]
#        print cats
        if len(cats) <= 0:
            cats = Mapforcing.loc[Mapforcing['SubId'] == catid]
            arcpy.AddMessage("Following Grid has to be inluded:.......")
            arcpy.AddMessage(cats['FGID'])
        tarea = sum(cats['s_area'].values)
        fids = cats['FGID'].values
        fids = np.unique(fids)
        sumwt = 0.0
#        arcpy.AddMessage(cats)
        for j in range(len(fids)):
            scat = cats[cats['FGID'] == fids[j]]
#            arcpy.AddMessage(str(j) + str(fids[j]))
#            arcpy.AddMessage(scat)
            if j < len(fids) - 1:
                sarea = sum(scat['s_area'].values)
                wt = float(sarea)/float(tarea)
                sumwt = sumwt + wt
            else:
                wt = 1- sumwt
#            arcpy.AddMessage(scat)
            if(len(scat['Row'].values) > 1):
                arcpy.AddMessage(str(catid)+"need check.......")
            Strcellid = str(int(scat['Row'].values[0] * (max(Forcinfo['Col'].values) + 1 +misscol) + scat['Col'].values[0])) + "      "
                ### str((ncrowcol[0,0] * ncncols + ncrowcol[0,1]))
            ogridforc.write("    "+str(int(catid)) + "     "+Strcellid+str(wt) +"\n")
            maparray[imap,0] = int(catid)
            maparray[imap,1] = int(Strcellid)
            maparray[imap,2] = wt
            maparray[imap,3] = scat['Row'].values[0]
            maparray[imap,4] = scat['Col'].values[0]
            imap = imap + 1
#        arcpy.AddMessage(cats)
        if cats['IsLake'].values[0] > 0:
#        arcpy.AddMessage(Mapforcing.loc[Mapforcing['SubId'] == catid])
            tarea = sum(cats['s_area'].values)
            fids = cats['FGID'].values
            fids = np.unique(fids)
            sumwt = 0.0
            for j in range(len(fids)):
                scat = cats[cats['FGID'] == fids[j]]
                if j < len(fids) - 1:
                    sarea = sum(scat['s_area'].values)
                    wt = float(sarea)/float(tarea)
                    sumwt = sumwt + wt
                else:
                    wt = 1- sumwt
#            arcpy.AddMessage(scat)
                if(len(scat['Row'].values) > 1):
                    arcpy.AddMessage(str(catid)+"need check......")
                Strcellid = str(int(scat['Row'].values * (max(Forcinfo['Col'].values)+1 + misscol) + scat['Col'].values)) + "      "
                ### str((ncrowcol[0,0] * ncncols + ncrowcol[0,1]))
                ogridforc.write("    "+str(int(catid) + int(maxcatid)) + "     "+Strcellid+str(wt) +"\n")
    ogridforc.write(":EndGridWeights")
    ogridforc.close()
    return maparray

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
forcingply =  WorkingFolder + "Gridncply.shp"
forcinggrid = "#"
Boundaryply = "#"
missrow = 0
misscol= 0

maparray = Maphru2forceply(forcingply,outfolder,catply,WorkingFolder,Boundaryply,missrow,misscol)

#########Generate HYPE DAILY INPUTS
varnames = ['PR0_SFC','FB_SFC','Tave','Tmin','Tmax']
varnames = ['PR0_SFC','FB_SFC','Tave','Tmin','Tmax']
outfilenm = ['Pobs.txt','Rnobs.txt','Tobs.txt','TMINobs.txt','TMAXobs.txt']
unitscale = [1000,1,1,1,1]
dbf1 = Dbf5(catply[:-3]+'dbf')
catinfo = dbf1.to_dataframe()
print len(catinfo['SubId'].values),len(np.unique(catinfo['SubId'].values))

ncin =  Dataset(ncfile,'r')
datevar = []
t_unit  = ncin.variables['time'].units
nctime = ncin.variables['time'][:]
t_cal =ncin.variables['time'].calendar
datevar = [num2date( x ,units = t_unit,calendar = t_cal) for x in nctime]

for ivar in range(0,len(varnames)):
    varname = varnames[ivar]
    outdatahr_var = pd.DataFrame(np.nan,index = datevar,columns=["date"])
    for icat in range(0,len(catinfo['SubId'].values)):
        arcpy.AddMessage(str(varname) + "     " + str(icat))
        catid = int(catinfo['SubId'].values[icat])
        outdatahr_var[str(catid)] = np.nan
        imaparray = maparray[maparray[:,0] == catid,]
        sumvalues = np.full(len(nctime),0.000000)
        for m in range(0,len(imaparray)):
            mrow = int(imaparray[m,3])
            mcol = int(imaparray[m,4])
            mwt = imaparray[m,2]
            if varname == 'Tave' or varname == 'Tmax' or varname == 'Tmin':
                 data = ncin.variables['TT_40m'][:,mrow,mcol]
                 sumvalues[:] = sumvalues[:] + data*mwt*unitscale[ivar]
            elif varname == 'FB_SFC':
                data1 = ncin.variables[varname][:,mrow,mcol]
                data2 = ncin.variables['FI_SFC'][:,mrow,mcol]
                data1[data1<0] = 0.0
                data2[data2<0] = 0.0
                data = data1 + data2
                sumvalues[:] = sumvalues[:] + data*mwt*unitscale[ivar]
            else:
                data = ncin.variables[varname][:,mrow,mcol]
                data[data<0] = 0.0
                sumvalues[:] = sumvalues[:] + data*mwt*unitscale[ivar]
        outdatahr_var.loc[:,str(catid)] =sumvalues
    if varname == 'Tave':
        outdataday_var = outdatahr_var.resample('D').mean()
        outdatahr_var.to_csv(WorkingFolder+"checkT.csv",sep='\t',index = None)
    elif varname == 'Tmin':
        outdataday_var = outdatahr_var.resample('D').min()
    elif varname == 'Tmax':
        outdataday_var = outdatahr_var.resample('D').max()
    else:
        outdataday_var = outdatahr_var.resample('D').sum()
    outdataday_var.loc[:,'date'] = outdataday_var.index[:].strftime('%Y-%m-%d')
    outdataday_var.to_csv(WorkingFolder+outfilenm[ivar],sep='\t',index = None)
ncin.close()
