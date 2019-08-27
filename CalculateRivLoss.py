
import numpy as np
import sys
import os
import csv
import pandas as pd
from simpledbf import Dbf5
import matplotlib
import matplotlib.pyplot as plt

##### Readed inputs
#POI_FIDS_file = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/poi.csv"
#OutputF = "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/Outputs1/"
#RefFNam = 'RivLloss_ACC_1000'

POI_FIDS_file = sys.argv[1]
RefFNam = sys.argv[2]
OutputF = sys.argv[3] + '/'


pois =  pd.read_csv(POI_FIDS_file,sep=',')### read fid of pois
pois_new =  pd.read_csv(POI_FIDS_file,sep=',')
pois_new2 =  pd.read_csv(POI_FIDS_file,sep=',')

pois['Total_DA'] = np.nan
pois_new['Total_DA'] = np.nan
pois_new2['Total_DA'] = np.nan

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
    pois_new[Allfolders[i]] = 0.00
    pois_new2[Allfolders[i]] = 0.00

#### loop for each poi
    for j in range(0,len(pois)):
        ipoi = pois['POI'][j]

        tsubid = hyshdinfo2.loc[hyshdinfo2['IsObs'] == ipoi,]['SubId']

        if len(tsubid) > 0:
            ### define output folder for catchment of jth point of interest
            ij_shp =  OutputF +Allfolders[i]+'/'+'POI_'+str(ipoi)+'/'+'finalcat_info.shp'
            if not os.path.exists(ij_shp):
                continue
            ij_tempinfo = Dbf5(ij_shp[:-3]+'dbf')
            ij_info = ij_tempinfo.to_dataframe()
            if np.isnan(pois.loc[j,'Total_DA']):
                pois.loc[j,'Total_DA'] = np.sum(ij_info['Area2'].values) ### calcuate DA for each point of interest
                pois_new.loc[j,'Total_DA'] = np.sum(ij_info['Area2'].values)
                pois_new2.loc[j,'Total_DA'] = np.sum(ij_info['Area2'].values)
            for k in range(0,len(ij_info)):
                catid = ij_info['SubId'].values[k]
                upcatid = ij_info.loc[ij_info['DowSubId'] == catid,]['SubId'].values
                if len(upcatid) > 0:
                    pois.loc[j,Allfolders[i]] = pois.loc[j,Allfolders[i]] + max(ij_info['Rivlen'].values[k],0)*ij_info['Area2'].values[k]
                    pois_new2.loc[j,Allfolders[i]] = pois_new2.loc[j,Allfolders[i]] + max(ij_info['Rivlen'].values[k],0)
                else:
                    pois.loc[j,Allfolders[i]] = pois.loc[j,Allfolders[i]] + 0.0 ##upstream cat
                    pois_new2.loc[j,Allfolders[i]] = pois_new2.loc[j,Allfolders[i]] + 0.0 ##upstream cat
############calculate new methods
                ccatid = catid
                downcatid = ij_info.loc[ij_info['SubId'] == ccatid,]['DowSubId'].values
                while len(downcatid) > 0 and downcatid[0] > 0 and len(ij_info.loc[ij_info['SubId'] == downcatid[0],]['Rivlen'].values) > 0:
#                    arcpy.AddMessage(max(ij_info.loc[ij_info['SubId'] == downcatid[0],]['Rivlen'].values,0))
                    pois_new.loc[j,Allfolders[i]] = pois_new.loc[j,Allfolders[i]] + max(ij_info.loc[ij_info['SubId'] == downcatid[0],]['Rivlen'].values,0)*ij_info['Area2'].values[k]
                    ccatid = downcatid[0]
                    downcatid = ij_info.loc[ij_info['SubId'] == ccatid,]['DowSubId'].values
############
            pois.loc[j,Allfolders[i]] = pois.loc[j,Allfolders[i]]/pois.loc[j,'Total_DA']
            pois_new.loc[j,Allfolders[i]] = pois_new.loc[j,Allfolders[i]]/pois_new.loc[j,'Total_DA']
        else:
            print('Point of interest is not inculded in subbasins:     ' + str(ipoi))
pois.to_csv(OutputF+'rivweightlength_hongli.csv')
pois_new.to_csv(OutputF+'rivweightlength_new.csv')
pois_new2.to_csv(OutputF+'rivwlength_alone.csv')
for i in range(0,len(Allfolders)):
    if RefFNam == Allfolders[i] or Allfolders[i] not in pois.columns:
        continue
    else:
        pois[Allfolders[i]] = pois[RefFNam].values - pois[Allfolders[i]].values
        pois_new[Allfolders[i]] = pois_new[RefFNam].values - pois_new[Allfolders[i]].values
        pois_new2[Allfolders[i]] = pois_new2[RefFNam].values - pois[Allfolders[i]].values
pois.to_csv(OutputF+'rivlenloss_hongli.csv')
pois_new.to_csv(OutputF+'rivlenloss_new.csv')
pois_new2.to_csv(OutputF+'rivlenloss_rivalone.csv')

