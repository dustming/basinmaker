
import numpy as np
import sys
import os
import csv
import pandas as pd
from simpledbf import Dbf5


##### Readed inputs
POI_FIDS_file = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/poi.csv"
OutputF = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Outputs/"
RefFNam = 'RivLloss_ACC_100'

pois =  pd.read_csv(POI_FIDS_file,sep=',')### read fid of pois

pois['Total_DA'] = np.nan

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
    
    pois[Allfolders[i]] = 0.00
#### loop for each poi
    for j in range(0,len(pois)):
        ipoi = pois['POI'][j]

        tsubid = hyshdinfo2.loc[hyshdinfo2['IsObs'] == ipoi,]['SubId']

        if len(tsubid) > 0:
            ### define output folder for catchment of jth point of interest
            ij_shp =  OutputF +Allfolders[i]+'/'+'POI_'+str(ipoi)+'/'+'finalcat_info.shp'
            ij_tempinfo = Dbf5(hyshdply[:-3]+'dbf')
            ij_info = ij_tempinfo.to_dataframe()
            
            for k in range(0,len(ij_info)):
                pois.loc[j,'Total_DA'] = np.sum(ij_info['Area2'])
                pois.loc[j,Allfolders[i]] = pois.loc[j,Allfolders[i]] + ij_info['Rivlen'].values[k]*ij_info['Area2'].values[k]
            pois.loc[j,Allfolders[i]] = pois.loc[j,Allfolders[i]].values[0]/pois.loc[j,'Total_DA'].values[0]
        else:
            print('Point of interest is not inculded in subbasins:     ' + str(ipoi))



