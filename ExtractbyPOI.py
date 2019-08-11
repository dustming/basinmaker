
import numpy as np
import arcpy
from arcpy import env
from arcpy.sa import *
import sys
import os
import csv
import pandas as pd
from simpledbf import Dbf5
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

##### Readed inputs
POI_FIDS_file = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/poi.csv"
OutputF = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Outputs/"
Pathoftoolbox = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Code/Toolbox/An automated ArcGIS toolbox for watershed delineation with lakes.tbx"

arcpy.ImportToolbox(Pathoftoolbox)


pois =  pd.read_csv(POI_FIDS_file,sep=',')### read fid of pois

#### check if need process with several delineation results
if os.path.exists(OutputF + 'finalcat_info.shp'):
    Allfolders = ['']
else:
    Allfolders = os.listdir(OutputF)
    
for i in range(0,len(Allfolders)):
    hyshdply = OutputF +Allfolders[i] + '/'+'finalcat_info.shp'
    hyshdply.replace('//','/')
    if not os.path.exists(OutputF + Allfolders[i]+'/'+ 'finalcat_info.shp'):
        print (OutputF + Allfolders[i]+'/'+ 'finalcat_info.shp')
        continue
    ####read dbf file of the catchment of ith acc
    tempinfo = Dbf5(hyshdply[:-3]+'dbf')#np.genfromtxt(hyinfocsv,delimiter=',')
    hyshdinfo2 = tempinfo.to_dataframe()

#### loop for each poi
    for j in range(0,len(pois)):
        ipoi = pois['POI'][j]

        tsubid = hyshdinfo2.loc[hyshdinfo2['IsObs'] == ipoi,]['SubId']

        if len(tsubid) > 0:
            ### define output folder for catchment of jth point of interest
            OutputfoldercaseRes =  OutputF +Allfolders[i]+'/'+'POI_'+str(ipoi)+'/'
            if not os.path.exists(OutputfoldercaseRes):
                os.makedirs(OutputfoldercaseRes)
            print(tsubid,ipoi,hyshdply)
            arcpy.AutomatedLocalRoutingNetworkExtractionToolset(str(int(tsubid.values[0])),hyshdply,OutputfoldercaseRes,'-1')
        else:
            print('Point of interest is not inculded in subbasins:     ' + str(ipoi))



