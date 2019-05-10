import os
import pandas as pd
from simpledbf import Dbf5
import numpy as np
####### Required parameters
countthreshold = 50
#WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Grand River Basin/'
WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/Dataprocess/NithRiver/'
modelname = 'NithCrop'
soilaynm = ['_Top','_Fast','_Slow']
soiltype = ['Soil1','Soil3','Soil4','Soil5','Soil8','Soil9']
soilwt = pd.read_csv(WorkFolder+'SoilParweight.csv')
soilpars = pd.read_csv(WorkFolder+'SoilParm.csv')
orvh2 = open(WorkFolder + 'RavenInputs/'+modelname+".rvp","w")
f1 = open(WorkFolder+'NithCrop.rvp_sample',"r")
tab = "         "
soilpsec = 0
for x in f1:
    chanpar = 0.0
    if ":RainSnowTransition" in x:
        orvh2.write(":RainSnowTransition" + tab + str(soilpars['RainSnowTransition'].values[0]) + tab + "1" + "\n")
        chanpar = 1
    if ":SoilParameterList" in x:
        soilpsec = 1
        orvh2.write(x)
        continue
    if ":EndSoilParameterList" in x:
        soilpsec = 0
        orvh2.write(x)
        continue

    if "Soil1_Top" in x and soilpsec == 1:
        orvh2.write(":Parameters         POROSITY         FIELD_CAPACITY         SAT_WILT         HBV_BETA         MAX_PERC_RATE         BASEFLOW_COEFF         BASEFLOW_N         MAX_CAP_RISE_RATE"+"\n")
        orvh2.write(":Units          none         none         none         none         mm/d         1/d         none         mm/d"+"\n")
        for itp in range(0,len(soiltype)):
            for ily in range(0,len(soilaynm)):
                Soilname = soiltype[itp] + soilaynm[ily]
                wt = soilwt[soilwt['Soil Name'] == Soilname]['Sand'].values[0]
                if ily == 0:
                    POROSITY = str(round(soilpars['POROSITYa'].values[0]*wt + soilpars['POROSITYb'].values[0],5)) + tab ### for top layer
                    FIELD_CAPACITY = str(round(soilpars['FIELD_CAPACITYa'].values[0]*wt + soilpars['FIELD_CAPACITYb'].values[0],5)) + tab ### for top layer
                    SAT_WILT  = str(round(soilpars['SAT_WILTa'].values[0]*wt + soilpars['SAT_WILTb'].values[0],5))+ tab  ### for top layer
                    HBV_BETA  = str(round(soilpars['HBV_BETAa'].values[0]*wt + soilpars['HBV_BETAb'].values[0],5)) + tab ### for top layer
                    MAX_PERC_RATE = str(round(0,5)) + tab
                    BASEFLOW_COEFF = str(round(0,5)) + tab
                    BASEFLOW_N = str(round(0,5)) + tab
                    MAX_CAP_RISE_RATE = str(round(0,5)) + tab
                else:
                    POROSITY = str(round(soilpars['POROSITYa'].values[1]*wt + soilpars['POROSITYb'].values[1],5)) + tab ### for top layer
                    FIELD_CAPACITY = str(round(soilpars['FIELD_CAPACITYa'].values[1]*wt + soilpars['FIELD_CAPACITYb'].values[1],5)) + tab ## for rest two layer
                    SAT_WILT  = str(round(soilpars['SAT_WILTa'].values[1]*wt + soilpars['SAT_WILTb'].values[1],5))+ tab  ### for top layer
                    HBV_BETA  = str(round(soilpars['HBV_BETAa'].values[1]*wt + soilpars['HBV_BETAb'].values[1],5)) + tab ### for top layer
                    if ily == 1:
                        MAX_PERC_RATE = str(round(soilpars['MAX_PERC_RATEa'].values[0]*wt + soilpars['MAX_PERC_RATEb'].values[0],5)) + tab
                        BASEFLOW_COEFF = str(round(soilpars['BASEFLOW_COEFFa'].values[0]*wt + soilpars['BASEFLOW_COEFFb'].values[0],5)) + tab
                        BASEFLOW_N = str(round(soilpars['BASEFLOW_Na'].values[0]*wt + soilpars['BASEFLOW_Nb'].values[0],5)) + tab
                        MAX_CAP_RISE_RATE = str(round(soilpars['MAX_CAP_RISE_RATEa'].values[0]*wt + soilpars['MAX_CAP_RISE_RATEb'].values[0],5)) + tab
                    if ily == 2:
                        MAX_PERC_RATE = str(round(soilpars['MAX_PERC_RATEa'].values[1]*wt + soilpars['MAX_PERC_RATEb'].values[1],5)) + tab
                        BASEFLOW_COEFF = str(round(soilpars['BASEFLOW_COEFFa'].values[1]*wt + soilpars['BASEFLOW_COEFFb'].values[1],5)) + tab
                        BASEFLOW_N = str(round(soilpars['BASEFLOW_Na'].values[1]*wt + soilpars['BASEFLOW_Nb'].values[1],5)) + tab
                        MAX_CAP_RISE_RATE = str(round(soilpars['MAX_CAP_RISE_RATEa'].values[1]*wt + soilpars['MAX_CAP_RISE_RATEb'].values[1],5)) + tab                        
                orvh2.write(Soilname+tab+POROSITY+FIELD_CAPACITY+SAT_WILT+HBV_BETA+
                            MAX_PERC_RATE+BASEFLOW_COEFF+BASEFLOW_N+MAX_CAP_RISE_RATE+"\n")
        chanpar = 1
        continue
    if chanpar == 0.0 and soilpsec == 0:
        orvh2.write(x)
f1.close()
orvh2.close
#runmodel = WorkFolder + 'RavenInputs/BaseModel/'+'Raven.exe    NithCrop'
#os.system(runmodel)
