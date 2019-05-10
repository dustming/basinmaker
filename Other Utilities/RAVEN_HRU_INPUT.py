def combinervhrvp(foldername):
    orvh1 = open(Ourfolder+foldername+ modelname + ".rvh","w")
    f1 = open(Ourfolder+"test.rvh","r")
    for x in f1:
        orvh1.write(x)
        if ":EndSubBasins" in x:
            break
    f1.close()
    f1 = open(Ourfolder+"SubBasinProp.rvh","r")
    for x in f1:
        orvh1.write(x)
    f1.close()
    f1 = open(Ourfolder+"hru.rvh","r")
    for x in f1:
        orvh1.write(x)
    f1.close()
    orvh1.close()  

    orvh1 = open(Ourfolder+foldername+ modelname + ".rvp","w")
    f1 = open(Ourfolder+"soil.rvp","r")
    for x in f1:
        orvh1.write(x)
    f1.close()
    f1 = open(Ourfolder+"landuse.rvp","r")
    for x in f1:
        orvh1.write(x)
    f1.close()
    orvh1.close()
    
def Generatesoilrvpfile(Ourfolder,catinfo,hrus,landuseinfo,soilinfo):
    modelname = 'soil'
    tab2 = ',    '
    catids = np.unique(catinfo['SubId'].values)
    hruin = hrus[hrus['CATCHMENTS'].isin(catids)]
    orvh2 = open(Ourfolder+modelname+".rvp","w")
    orvh2.write("# --------------------------------------------  "+"\n")
    orvh2.write("# Raven Classed Parameter File"+"\n")
    orvh2.write("# For soil, land use, vegetation, and aquifer classes   " + "\n")
    orvh2.write("# HBV-EC Grand River emulation:" + "\n")
    orvh2.write("# Author: Ming Han Modified from M. Shafii    " + "\n")  
    orvh2.write("# --------------------------------------------    " + "\n")
    orvh2.write(":RainSnowTransition 1 1.7395       " + "\n")
    orvh2.write(":AdiabaticLapseRate 6.5     " + "\n")
    orvh2.write(":IrreducibleSnowSaturation 0.05    " + "\n")
    orvh2.write(":AvgAnnualRunoff 100    " + "\n")
    orvh2.write("# ----Soil Classes----------------------------     " + "\n")            
    orvh2.write(":SoilClasses          " + "\n") 
    orvh2.write(":Attributes		%SAND	%CLAY	%SILT	%ORGANIC" + "\n") 
    orvh2.write("	:Units			none	none	none	none" + "\n")  
    for i in range(0,len(soilinfo)):
        for j in range(0,len(soilaynm)):
            soilarea = sum(hruin[hruin['SOILCA'] == soilinfo['Soilid'].values[i]]['COUNT'].values)*30*30/1000.0/1000.0
#            print(soilarea,soilinfo['Soilid'].values[i])
            if soilarea <= 0.0:
                continue
            soillyname = soilinfo['SoilName'].values[i] + soilaynm[j]
            sand = soilinfo['VFSAND'].values[i]/100.0 + soilinfo['TSAND'].values[i]/100.0
            clay = soilinfo['TCLAY'].values[i]/100.0
            silt = soilinfo['TSILT'].values[i]/100.0
            organ = soilinfo['ORGCARB'].values[i]/100.0
            sum1 = sand+clay+ silt
            orvh2.write(tab+soillyname+tab2 +str(sand/sum1)+tab2+str(clay/sum1)+tab2+str(silt/sum1)+tab2
                        +str(organ)+tab2+ "\n")
    orvh2.write(":EndSoilClasses" + "\n")
    orvh2.write("# --- Soil Profiles ----------------------------"+ "\n")
    orvh2.write(":SoilProfiles"+ "\n")
    for i in range(0,len(soilinfo)):
        soilarea = sum(hruin[hruin['SOILCA'] == soilinfo['Soilid'].values[i]]['COUNT'].values)*30*30/1000.0/1000.0
#            print(soilarea,soilinfo['Soilid'].values[i])
        if soilarea <= 0.0:
            continue
        orvh2.write(tab+ soilinfo['SoilName'].values[i]+tab2+str(3)
        +tab2+soilinfo['SoilName'].values[i] + soilaynm[0]+tab2+str(0.15)
        +tab2+soilinfo['SoilName'].values[i] + soilaynm[1]+tab2+str(0.4)
        +tab2+soilinfo['SoilName'].values[i] + soilaynm[2]+tab2+str(5)+ "\n")
    orvh2.write("       LAKE,    0"+ "\n")
    orvh2.write(":EndSoilProfiles"+"\n")

    orvh2.write(":SoilParameterList"+ "\n")
    orvh2.write(":Parameters,	POROSITY,		FIELD_CAPACITY,	SAT_WILT,	HBV_BETA,	MAX_PERC_RATE,	"+
                "BASEFLOW_COEFF,	BASEFLOW_N,		"+
                "PET_CORRECTION, MAX_CAP_RISE_RATE"+ "\n")

    orvh2.write(":Units     ,	none,			none,			none,		none,		mm/d,"+
                "			1/d,				none,			"+
                "		    	fract,"+"   mm/d"+ "\n")

    orvh2.write("[DEFAULT],		1.0,			1.0,			0.0,		0.57569,		1.446,"+
                "		0.044002,	1.0,"+
                "		0.58314,"+"   0.1  "+ "\n")

# =============================================================================
    for i in range(0,len(soilinfo)):
        for j in range(0,len(soilaynm)):
            soilarea = sum(hruin[hruin['SOILCA'] == soilinfo['Soilid'].values[i]]['COUNT'].values)*30*30/1000.0/1000.0
#            print(soilarea,soilinfo['Soilid'].values[i])
            if soilarea <= 0.0:
                continue
            soillyname = soilinfo['SoilName'].values[i] + soilaynm[j] + tab2
            POROSITY = str(soilinfo['KP0'].values[i]/100.00) + tab2
            tt = soilinfo['KP0'].values[i]/100.00
            FIELD_CAPACITY = str(soilinfo['KP33'].values[i]/100.00/tt) + tab2
            SAT_WILT = str(soilinfo['KP1500'].values[i]/100.00/tt) + tab2
            HBV_BETA = '_DEFAULT' + tab2
            MAX_PERC_RATE = '_DEFAULT' + tab2
            BASEFLOW_COEFF = '_DEFAULT' + tab2
            MAX_BASEFLOW_RATE = '_DEFAULT' + tab2
            BASEFLOW_N = '_DEFAULT' + tab2
            BASEFLOW_THRESH = '_DEFAULT' + tab2
            PET_CORRECTION = '_DEFAULT' + tab2
            MAX_CAP_RISE_RATE = '_DEFAULT' + tab2
            orvh2.write(soillyname+POROSITY+FIELD_CAPACITY+SAT_WILT+HBV_BETA+MAX_PERC_RATE+
                        BASEFLOW_COEFF+ MAX_BASEFLOW_RATE+BASEFLOW_N+
                        BASEFLOW_THRESH+PET_CORRECTION+MAX_CAP_RISE_RATE+"\n")
# =============================================================================
        
    orvh2.write(":EndSoilParameterList"+"\n")
    orvh2.close()
    
def Generatervhfiles(Ourfolder,catinfo,hrus,landuseinfo,soilinfo):
    orvh1 = open(Ourfolder+"SubBasinProp.rvh","w")
    orvh1.write(":SubBasinProperties"+"\n")
    orvh1.write(":Parameters	TIME_CONC	TIME_TO_PEAK	TIME_LAG"+"\n")
    orvh1.write(":Units	d	d	d" + "\n")
    for i in range(0,len(catinfo)):
        subid = str(catinfo['SubId'].values[i]) + tab
        rivlen = catinfo['Rivlen'].values[i]
        timconc = str(rivlen/1000.0*0.0102) + tab
        timpeak = str(rivlen/1000.0*0.0058 + 0.0061) + tab
        if(rivlen/1000.0*0.0102 < rivlen/1000.0*0.0058 + 0.0061):
            timconc = str(rivlen/1000.0*0.0058 + 0.0061 + 0.01) + tab
        timelag = str(rivlen/1000.0*0.0058 + 0.0061) + tab
        orvh1.write(subid + timconc + timpeak + timelag + "\n")
    orvh1.write(":EndSubBasinProperties"+"\n")
    orvh1.close()
    orvh = open(Ourfolder+"hru.croprvh","w")
    orvh.write(":HRUs"+"\n")
    orvh.write("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT"+"\n")
    orvh.write("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg"+"\n")
    cats = np.unique(catinfo['SubId'].values)
    hruidk = 0
    for i in range(0,len(cats)):
        catid = cats[i]
        cathrus = hrus[hrus['CATCHMENTS'] == catid]
        icat = catinfo[catinfo['SubId'] == catid]
        catslope = str(max(icat['BasinSlope'].values[0],0.0001)) + tab
        lat = str(icat['INSIDE_Y'].values[0])+tab
        lon = str(icat['INSIDE_X'].values[0])+tab
        AQUIFER_PROFILE ='[NONE]'+tab
        TERRAIN_CLASS ='[NONE]'+tab
        subid = str(int(catid)) + tab
        ASPECT = '200'+tab
        StrGidelev = str(icat['MeanElev'].values[0]) + tab
        waterhruarea = 0.0
        otherarea = 0.0
        for j in range(len(cathrus.index)):
            idx = cathrus.index[j]
            hruid = str(hruidk) + tab
            soilid = cathrus.loc[idx]['SOILCA']
            landuseid = cathrus.loc[idx]['LANDUSE']
            isoil = soilinfo[soilinfo['Soilid'] == soilid]
            ilanduse = landuseinfo[landuseinfo['LandID'] == landuseid]
            
            hruarea = str(cathrus.loc[idx]['COUNT']*cellsize*cellsize/1000/1000) + tab
            LAND_USE_CLASS = str(ilanduse['Landuse'].values[0]) + tab
            VEG_CLASS = str(ilanduse['Veg'].values[0]) +tab
            
            if LAND_USE_CLASS == "Water" + tab:
                SOIL_PROFILE = "LAKE" + tab
                waterhruarea = waterhruarea + cathrus.loc[idx]['COUNT']*cellsize*cellsize/1000/1000
            else:
                otherarea = otherarea + cathrus.loc[idx]['COUNT']*cellsize*cellsize/1000/1000
                SOIL_PROFILE =str(isoil['SoilName'].values[0]) + tab 
                hruidk = hruidk + 1
                orvh.write("  "+hruid+hruarea+StrGidelev+lat+lon+subid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+catslope+ASPECT+"\n")
#### write Water hru,
        if waterhruarea > 0:
            waterhruarea = max(icat['Area2'].values[0]/1000/1000 - otherarea,0.1)
            hruid = str(hruidk) + tab
            hruarea = str(waterhruarea) + tab
            LAND_USE_CLASS = "Water" + tab
            VEG_CLASS = "Water" + tab
            SOIL_PROFILE = "LAKE" + tab
            hruidk = hruidk + 1
            orvh.write("  "+hruid+hruarea+StrGidelev+lat+lon+subid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+catslope+ASPECT+"\n")
    print("total hru :",  hruidk - 1 )     
    orvh.write(":EndHRUs"+"\n")
    orvh.close()
#########################start the program

def Generatecrophru(Ourfolder,catinfo,hrus,landuseinfo,soilinfo):
    orvh1 = open(Ourfolder+"crop.hrucrop","w")
    orvh1.write("subid	hruid	rotationname	temp_ave"+"\n")
    cats = np.unique(catinfo['SubId'].values)
    hruidk = 0
    for i in range(0,len(cats)):
        catid = cats[i]
        subid = str(catid) + tab
        cathrus = hrus[hrus['CATCHMENTS'] == catid]
        subid = str(int(catid)) + tab
        for j in range(len(cathrus.index)):
            idx = cathrus.index[j]
            hruid = str(hruidk + 1) + tab
            
            landuseid = cathrus.loc[idx]['LANDUSE']
            ilanduse = landuseinfo[landuseinfo['LandID'] == landuseid]  
            LAND_USE_CLASS = str(ilanduse['Landuse'].values[0]) + tab
            VEG_CLASS = str(ilanduse['Veg'].values[0]) +tab
            
            if LAND_USE_CLASS == "Water" + tab:
                continue
            else:
                orvh1.write("  "+subid+hruid+VEG_CLASS + "\n")   
                hruidk = hruidk + 1
    print("total hru :",  hruidk - 1 )     
    orvh1.close()
    
################################3
    
def Generatesoiltplfile(Ourfolder,catinfo,hrus,landuseinfo,soilinfo):
    modelname = 'soil'
    tab2 = ',    '
    hruin = hrus[hrus['CATCHMENTS'].isin(catids)]
    orvh2 = open(Ourfolder+modelname+".rvp","w")
    orvh2.write("# --------------------------------------------  "+"\n")
    orvh2.write("# Raven Classed Parameter File"+"\n")
    orvh2.write("# For soil, land use, vegetation, and aquifer classes   " + "\n")
    orvh2.write("# HBV-EC Grand River emulation:" + "\n")
    orvh2.write("# Author: Ming Han Modified from M. Shafii    " + "\n")  
    orvh2.write("# --------------------------------------------    " + "\n")
    orvh2.write(":RainSnowTransition 1 1.7395       " + "\n")
    orvh2.write(":AdiabaticLapseRate 6.5     " + "\n")
    orvh2.write(":IrreducibleSnowSaturation 0.05    " + "\n")
    orvh2.write(":AvgAnnualRunoff 100    " + "\n")
    orvh2.write("# ----Soil Classes----------------------------     " + "\n")            
    orvh2.write(":SoilClasses          " + "\n") 
    orvh2.write(":Attributes		%SAND	%CLAY	%SILT	%ORGANIC" + "\n") 
    orvh2.write("	:Units			none	none	none	none" + "\n")  
    for i in range(0,len(soilinfo)):
        for j in range(0,len(soilaynm)):
            soilarea = sum(hruin[hruin['SOILCA'] == soilinfo['Soilid'].values[i]]['COUNT'].values)*30*30/1000.0/1000.0
#            print(soilarea,soilinfo['Soilid'].values[i])
            if soilarea <= 0.0:
                continue
            soillyname = soilinfo['SoilName'].values[i] + soilaynm[j]
            sand = soilinfo['VFSAND'].values[i]/100.0 + soilinfo['TSAND'].values[i]/100.0
            clay = soilinfo['TCLAY'].values[i]/100.0
            silt = soilinfo['TSILT'].values[i]/100.0
            organ = soilinfo['ORGCARB'].values[i]/100.0
            sum1 = sand+clay+ silt
            orvh2.write(tab+soillyname+tab2 +str(sand/sum1)+tab2+str(clay/sum1)+tab2+str(silt/sum1)+tab2
                        +str(organ)+tab2+ "\n")
    orvh2.write(":EndSoilClasses" + "\n")
    orvh2.write("# --- Soil Profiles ----------------------------"+ "\n")
    orvh2.write(":SoilProfiles"+ "\n")
    for i in range(0,len(soilinfo)):
        soilarea = sum(hruin[hruin['SOILCA'] == soilinfo['Soilid'].values[i]]['COUNT'].values)*30*30/1000.0/1000.0
#            print(soilarea,soilinfo['Soilid'].values[i])
        if soilarea <= 0.0:
            continue
        orvh2.write(tab+ soilinfo['SoilName'].values[i]+tab2+str(3)
        +tab2+soilinfo['SoilName'].values[i] + soilaynm[0]+tab2+str(0.15)
        +tab2+soilinfo['SoilName'].values[i] + soilaynm[1]+tab2+str(0.4)
        +tab2+soilinfo['SoilName'].values[i] + soilaynm[2]+tab2+str(5)+ "\n")
    orvh2.write("       LAKE,    0"+ "\n")
    orvh2.write(":EndSoilProfiles"+"\n")

    orvh2.write(":SoilParameterList"+ "\n")
    orvh2.write(":Parameters,	POROSITY,		FIELD_CAPACITY,	SAT_WILT,	HBV_BETA,	MAX_PERC_RATE,	"+
                "BASEFLOW_COEFF,	BASEFLOW_N,		"+
                "PET_CORRECTION, MAX_CAP_RISE_RATE"+ "\n")

    orvh2.write(":Units     ,	none,			none,			none,		none,		mm/d,"+
                "			1/d,				none,			"+
                "		    	fract,"+"   mm/d"+ "\n")

    orvh2.write("[DEFAULT],		1.0,			1.0,			0.0,		0.57569,		1.446,"+
                "		0.044002,	1.0,"+
                "		0.58314,"+"   0.1  "+ "\n")

# =============================================================================
    for i in range(0,len(soilinfo)):
        for j in range(0,len(soilaynm)):
            soilarea = sum(hruin[hruin['SOILCA'] == soilinfo['Soilid'].values[i]]['COUNT'].values)*30*30/1000.0/1000.0
#            print(soilarea,soilinfo['Soilid'].values[i])
            if soilarea <= 0.0:
                continue
            soillyname = soilinfo['SoilName'].values[i] + soilaynm[j] + tab2
            POROSITY = str(soilinfo['KP0'].values[i]/100.00) + tab2
            tt = soilinfo['KP0'].values[i]/100.00
            FIELD_CAPACITY = str(soilinfo['KP33'].values[i]/100.00/tt) + tab2
            SAT_WILT = str(soilinfo['KP1500'].values[i]/100.00/tt) + tab2
            HBV_BETA = '_DEFAULT' + tab2
            MAX_PERC_RATE = '_DEFAULT' + tab2
            BASEFLOW_COEFF = '_DEFAULT' + tab2
            MAX_BASEFLOW_RATE = '_DEFAULT' + tab2
            BASEFLOW_N = '_DEFAULT' + tab2
            BASEFLOW_THRESH = '_DEFAULT' + tab2
            PET_CORRECTION = '_DEFAULT' + tab2
            MAX_CAP_RISE_RATE = '_DEFAULT' + tab2
            orvh2.write(soillyname+POROSITY+FIELD_CAPACITY+SAT_WILT+HBV_BETA+MAX_PERC_RATE+
                        BASEFLOW_COEFF+ MAX_BASEFLOW_RATE+BASEFLOW_N+
                        BASEFLOW_THRESH+PET_CORRECTION+MAX_CAP_RISE_RATE+"\n")
# =============================================================================
        
    orvh2.write(":EndSoilParameterList"+"\n")
    orvh2.close()

##### tempk

    
import os
import pandas as pd
from simpledbf import Dbf5
import numpy as np
####### Required parameters
countthreshold = 50
#WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Grand River Basin/'
WorkFolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/Dataprocess/'
InputsFolder =  WorkFolder + 'Project/'
obsshp = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/Data/Obspoint/Obsfinal_infof_inmodel_f.shp'
GRCAfolder = 'C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/GrandRiverProject/Data/GRIN/'
cellsize = 30
maxsubnum = 500
modelname = 'NithCrop'
soilaynm = ['_Top','_Fast','_Slow']
###################Read inputs
Ourfolder =WorkFolder + "NithRiver/" + "RavenInputs/"
if not os.path.exists(Ourfolder):
    os.makedirs(Ourfolder)
tab = '\t'
soilinfo = pd.read_csv(InputsFolder+"Soilinfo.csv",sep=",",low_memory=False)
landuseinfo = pd.read_csv(InputsFolder+"Landuseinfo.csv",sep=",",low_memory=False)
dbf2 = Dbf5(WorkFolder+ "NithRiver/"+"finalcat_info.dbf")
catinfo = dbf2.to_dataframe()
dbf1 = Dbf5(InputsFolder+ "HRU_COMBINE.dbf")
hrus = dbf1.to_dataframe()
routinfo = catinfo[['SubId','DowSubId']].values
dbf2 = Dbf5(obsshp[:-3]+'dbf')
obsinfo = dbf2.to_dataframe()
modelstart ='2010-01-01'
enddate = '2015-01-01'
index = pd.date_range(modelstart, enddate,freq='H')
lakeobsinfo = pd.read_csv(InputsFolder+"Lakeobs.csv",sep=",",low_memory=False)
########   define Raven Inputs
Generatervhfiles(Ourfolder,catinfo,hrus,landuseinfo,soilinfo)
#Generatesoilrvpfile(Ourfolder,catinfo,hrus,landuseinfo,soilinfo)
#combinervhrvp('BaseModel/')
#Generatecrophru(Ourfolder,catinfo,hrus,landuseinfo,soilinfo)