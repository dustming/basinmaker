# -*- coding: utf-8 -*-
"""
Created on Fri May  3 15:35:42 2019

@author: m43han
"""

import os
import pandas as pd
from simpledbf import Dbf5
import numpy as np
####### Required parameters
countthreshold = 50
#WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Grand River Basin/'
WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/Dataprocess/'
InputsFolder =  WorkFolder + 'NithRiver/'

hrucat = pd.read_csv(InputsFolder+ 'hrusub.csv')
subwt = pd.read_csv(InputsFolder+ 'subwt.csv')
orvh2 = open(InputsFolder+"GriddedForcings2.txt","w")
orvh2.write(":GridWeights"+ "\n")
orvh2.write("   :NumberHRUs       1731"+ "\n")
orvh2.write("   :NumberGridCells  144"+ "\n")
orvh2.write("   # [HRU ID] [Cell #] [w_kl]"+ "\n")
tab = "   "
for i in range(0,len(hrucat)):
    hruid = tab + str(hrucat['hruid'].values[i]) + tab
    catid = hrucat['subid'].values[i]
    isubwt = subwt[subwt['hruid'] == catid]
    for j in range(0,len(isubwt)):
        gridid = str(isubwt['gridid'].values[j])+tab
        wt = str(isubwt['wt'].values[j])+tab
        orvh2.write(hruid + gridid + wt + "\n")
orvh2.write(":EndGridWeights"+ "\n")
orvh2.close()