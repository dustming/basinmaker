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
    

import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
from simpledbf import Dbf5
import pandas as pd
import copy

arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

Outputfolder1 = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Lakeofwoods_merit"
OutputFolder = Outputfolder1 + "/" + "Result8470/"





arcpy.CopyFeatures_management(OutputFolder + "/" + "finalcat_info.shp", OutputFolder + "/" + "finalcat_info2")

tempinfo = Dbf5(OutputFolder + "/" + "finalcat_info.dbf")#np.genfromtxt(hyinfocsv,delimiter=',')
hyshdinfo2 = tempinfo.to_dataframe()
hyshdinfo = hyshdinfo2[['SubId','DowSubId']].values

dbfile = OutputFolder + 'finalcat_info2.shp'
inFeatures = dbfile
fieldPrecision = 10
field_scale = 3
arcpy.AddField_management(dbfile, "DA", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
arcpy.AddField_management(dbfile, "R_DAl_DAw", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
arcpy.AddField_management(dbfile, "R_Al_DAl", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
arcpy.AddField_management(dbfile, "NlRDAlDAw", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
arcpy.AddField_management(dbfile, "NlRAlDAl", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
arcpy.AddField_management(dbfile, "Lake_zone", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")

rows = arcpy.UpdateCursor(dbfile)
Watarea = sum(hyshdinfo2['Area2'].values)/1000.00/1000.00
k = 0
for row in rows:
    k = k +1
    print(k)
    ### obtain upstream catchment IDS
#    HydroBasins1 = Defcat(hyshdinfo,row.SubId)
#    hytemp = hyshdinfo2.loc[hyshdinfo2['SubId'].isin(HydroBasins1)]
    
    ### calculate upstream drainage area 
#    row.DA = sum(hytemp['Area2'].values)
    if row.LakeArea > 0:   ##3 check if it is a lake catchment
    
        HydroBasins1 = Defcat(hyshdinfo,row.SubId)
        hytemp = hyshdinfo2.loc[hyshdinfo2['SubId'].isin(HydroBasins1)]
    
    ### calculate upstream drainage area 
        row.DA = sum(hytemp['Area2'].values)/1000.00/1000.00
 
        row.R_DAl_DAw = row.DA / Watarea
        row.R_Al_DAl = min(row.LakeArea / row.DA,1)
        row.NlRDAlDAw = -np.log10(row.DA / Watarea)
        row.NlRAlDAl = -np.log10(min(row.LakeArea / row.DA,1))
    
        x = row.NlRDAlDAw
        y = row.NlRAlDAl
#        print(x,y,row.DA,row.LakeArea,row.R_DAl_DAw,row.R_Al_DAl)
        for j in range (0,7):
            print(j,-x +j, y,-x + j + 1)
            if -x +j <= y and -x + j + 1 >= y:
                row.Lake_zone = j + 1
#                print(row.Lake_zone,j)
                break; 
    rows.updateRow(row)
del row
del rows
    
