
from qgis.core import *
import qgis
from qgis.analysis import QgsNativeAlgorithms
from simpledbf import Dbf5
import os
import sys
import numpy as np
import shutil
from shutil import copyfile
import tempfile
import copy
import pandas as pd


def writechanel(chname,chwd,chdep,chslope,orchnl,elev,floodn,channeln,iscalmanningn):
    ### Following SWAT instructions, assume a trapezoidal shape channel, with channel sides has depth and width ratio of 2. zch = 2
    zch = 2
    sidwd = zch * chdep ###river side width
    tab = "          "
    botwd = chwd - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*chwd
        sidwd = 0.5*0.5*chwd
        zch = (chwd - botwd)/2/chdep
    if iscalmanningn >= 0:
        mann = str(channeln)
    else:
        mann = str(0.035)
    zfld = 4 + elev
    zbot = elev - chdep
    sidwdfp = 4/0.25
    Channame = ":ChannelProfile"+tab+chname+tab
    orchnl.write(Channame+"\n")
    Chanslop = "  :Bedslope"+tab+str(chslope)
    orchnl.write(Chanslop+"\n")
    orchnl.write("  :SurveyPoints"+"\n")
    orchnl.write("    0"+tab+str(zfld)+"\n")
    orchnl.write("    "+str(sidwdfp)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd)+tab+str(zbot)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd + botwd)+tab+str(zbot)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(2*sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(zfld)+"\n")
    orchnl.write("  :EndSurveyPoints"+"\n")
    orchnl.write("  :RoughnessZones"+"\n")
    orchnl.write("    0" + tab + str(floodn) +"\n")
    orchnl.write("    " + str(sidwdfp + 2*chwd)+ tab + mann +"\n")
    orchnl.write("    " + str(sidwdfp + 2*chwd + 2*sidwd + botwd)+ tab + str(floodn) +"\n")
    orchnl.write("  :EndRoughnessZones"+"\n")
    orchnl.write(":EndChannelProfile"+"\n")
    orchnl.write("\n")
    orchnl.write("##############new channel ##############################\n")
#########################################################################################################33



def writelake(catinfo,outFolderraven):
    f2 = open(os.path.join(outFolderraven,"TestLake.rvh"),"w")
    tab = '       '
    maxcatid = max(catinfo['SubId'].values)
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HyLakeId'] > 0:
            lakeid = int(catinfo.iloc[i]['HyLakeId'])
            catid = catinfo.iloc[i]['SubId']
            if float(catinfo.iloc[i]['Area'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LakeArea']):
                A = float(catinfo.iloc[i]['Area'])*0.95
            else:
                A = float(catinfo.iloc[i]['LakeArea'])*1000*1000
#            A = catinfo.iloc[i]['LakeArea']*1000*1000
            h0 = catinfo.iloc[i]['LakeDepth']
            WeirCoe = 0.6
            hruid = int(catinfo.iloc[i]['SubId']) + int(maxcatid)
            Crewd = catinfo.iloc[i]['BkfWidth']
#            if slakeinfo.iloc[0]['Wshd_area'] < 6000 and slakeinfo.iloc[0]['Wshd_area'] > 0:
        ######write lake information to file
            f2.write(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## " +"\n")
            f2.write("  :SubBasinID  "+str(int(catid))+ "\n")
            f2.write("  :HRUID   "+str(int(hruid))+ "\n")
            f2.write("  :Type RESROUTE_STANDARD   "+"\n")
            f2.write("  :WeirCoefficient  "+str(WeirCoe)+ "\n")
            f2.write("  :CrestWidth "+str(Crewd)+ "\n")
            f2.write("  :MaxDepth "+str(h0)+ "\n")
            f2.write("  :LakeArea    "+str(A)+ "\n")
            f2.write(":EndReservoir   "+"\n")
            f2.write("#############################################"+"\n")
            f2.write("###New Lake starts"+"\n")
    f2.close()
    #### write lake input files for different lake zone
#    arcpy.AddMessage(catinfo.columns) 
    if 'LAKE_ZONE' in catinfo.columns:  ### write output file for each lake zone 
        maxzoneid = int(max(catinfo['LAKE_ZONE'].values))
        minzoneid = int(min(catinfo['LAKE_ZONE'].values[np.nonzero(catinfo['LAKE_ZONE'].values)]))
        arcpy.AddMessage('minium zone id is    ' + str(minzoneid)+'     maximum zone id is: ' + str(maxzoneid))
        for izone in range(minzoneid,maxzoneid+1):
            filename = "TestLake_"+str(izone)+".rvh"
            f3 = open(os.path.join(outFolderraven,filename),"w") ## open a file to save lake output
            tab = '       '
            maxcatid = max(catinfo['SubId'].values)
            catinfozone = catinfo.loc[catinfo['LAKE_ZONE'] != izone]
            arcpy.AddMessage("# of lakes not in zone  " + str(izone)+"  :   " +str(len(catinfozone)))
            for i in range(0,len(catinfozone.index)):
                if catinfozone.iloc[i]['HyLakeId'] > 0:
                    lakeid = int(catinfozone.iloc[i]['HyLakeId'])
                    catid = catinfozone.iloc[i]['SubId']
                    if float(catinfozone.iloc[i]['Area'])/1000.00/1000.00 <= float(catinfozone.iloc[i]['LakeArea']):
                        A = float(catinfozone.iloc[i]['Area'])*0.95
                    else:
                        A = float(catinfozone.iloc[i]['LakeArea'])*1000*1000
#                    A = catinfozone.iloc[i]['LakeArea']*1000*1000
                    h0 = catinfozone.iloc[i]['LakeDepth']
                    WeirCoe = 0.6
                    hruid = int(catinfozone.iloc[i]['SubId']) + int(maxcatid)
                    Crewd = catinfozone.iloc[i]['BkfWidth']
#            if slakeinfo.iloc[0]['Wshd_area'] < 6000 and slakeinfo.iloc[0]['Wshd_area'] > 0:
        ######write lake information to file
                    f3.write(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## " +"\n")
                    f3.write("  :SubBasinID  "+str(int(catid))+ "\n")
                    f3.write("  :HRUID   "+str(int(hruid))+ "\n")
                    f3.write("  :Type RESROUTE_STANDARD   "+"\n")
                    f3.write("  :WeirCoefficient  "+str(WeirCoe)+ "\n")
                    f3.write("  :CrestWidth "+str(Crewd)+ "\n")
                    f3.write("  :MaxDepth "+str(h0)+ "\n")
                    f3.write("  :LakeArea    "+str(A)+ "\n")
                    f3.write(":EndReservoir   "+"\n")
                    f3.write("#############################################"+"\n")
                    f3.write("###New Lake starts"+"\n")
            f3.close()

def Writervhchanl(ocatinfo,outFolder,lenThres,iscalmanningn):
    catinfo = copy.copy(ocatinfo)
#    print int(catinfo.iloc[0]['SUBID']),len(catinfo.index)
    ochn = open(os.path.join(outFolder,"modelchannel.rvp"),"w")
##################3
    orvh = open(os.path.join(outFolder,"test.rvh"),"w")
    orvh.write("# --------------------------------------------"+"\n")
    orvh.write("# Raven HRU Input file"+"\n")
    orvh.write("#  lake catchment emulation"+"\n")
    orvh.write("# --------------------------------------------"+"\n")
    orvh.write(":SubBasins"+"\n")
    orvh.write("  :Attributes   NAME  DOWNSTREAM_ID       PROFILE REACH_LENGTH  GAUGED"+"\n")
    orvh.write("  :Units        none           none          none           km    none"+"\n")
    tab = "     "
    for i in range(0,len(catinfo.index)):
        ### Get catchment width and dpeth
        catid = int(catinfo.iloc[i]['SubId'])
        temp = catinfo.iloc[i]['Rivlen']
        if (float(temp) >= lenThres):
            catlen = float(temp)/1000 #### in km
            strRlen = str(catlen)
        else:
            catlen = -9999
            strRlen = 'ZERO-'
        if catinfo.iloc[i]['IsLake'] >= 0 :
            strRlen = 'ZERO-'
        #####################################################3
        Strcat = str(catid)
        if catid == catinfo.iloc[i]['DowSubId']:
            StrDid = str(-1)
        else:
            StrDid = str(int(catinfo.iloc[i]['DowSubId']))
        pronam = 'Chn_'+ Strcat
        chslope = max(catinfo.iloc[i]['RivSlope'],0.00001)
        if chslope < 0:
            chslope = catinfo.iloc[i]['BasinSlope']
        writechanel(pronam,max(catinfo.iloc[i]['BkfWidth'],1),max(catinfo.iloc[i]['BkfDepth'],1),
                    chslope,ochn,catinfo.iloc[i]['MeanElev'],catinfo.iloc[i]['FloodP_n'],catinfo.iloc[i]['Ch_n'],iscalmanningn)
        if catinfo.iloc[i]['IsObs'] >= 0 :
            Guage = '1'
        else:
            Guage = '0'
        orvh.write("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage+"\n")
    orvh.write(":EndSubBasins"+"\n")
    orvh.write("\n")
##########################################
    orvh.write(":HRUs"+"\n")
    orvh.write("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT"+"\n")
    orvh.write("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg"+"\n")
    maxcatid = max(catinfo['SubId'].values)
    for i in range(0,len(catinfo.index)):
        hruid = int(catinfo.iloc[i]['SubId'])
        catslope = catinfo.iloc[i]['BasinSlope']
        if catinfo.iloc[i]['IsLake'] > 0:
            if float(catinfo.iloc[i]['Area'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LakeArea']):
                catarea2 = float(catinfo.iloc[i]['Area'])*0.05/1000.00/1000.00
            else:
                catarea2 = float(catinfo.iloc[i]['Area'])/1000.00/1000.00 - float(catinfo.iloc[i]['LakeArea'])
        else:
            catarea2 = float(catinfo.iloc[i]['Area'])/1000.00/1000.00
        StrGid =  str(hruid)+tab
        catid = str(int(catinfo.iloc[i]['SubId']))+tab
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(catinfo.iloc[i]['MeanElev'])+tab
        lat = str(catinfo.iloc[i]['centroid_x'])+tab
        lon = str(catinfo.iloc[i]['centroid_y'])+tab
        LAND_USE_CLASS = 'FOREST'+tab
        VEG_CLASS = 'FOREST'+tab
        SOIL_PROFILE ='SOILPROF'+tab
        AQUIFER_PROFILE ='[NONE]'+tab
        TERRAIN_CLASS ='[NONE]'+tab
        SLOPE = str(catslope)+tab
        ASPECT = '200'+tab
        orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
        if catinfo.iloc[i]['IsLake'] > 0:
            hruid = int(catinfo.iloc[i]['SubId']) + int(maxcatid)
            catslope = catinfo.iloc[i]['BasinSlope']
            if float(catinfo.iloc[i]['Area'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LakeArea']):
                catarea2 = float(catinfo.iloc[i]['Area'])*0.95/1000/1000
            else:
                catarea2 = float(catinfo.iloc[i]['LakeArea'])
            StrGid =  str(hruid)+tab
            catid = str(int(catinfo.iloc[i]['SubId']))+tab
            StrGidarea = str(catarea2)+tab
            StrGidelev = str(catinfo.iloc[i]['MeanElev'])+tab
            lat = str(catinfo.iloc[i]['centroid_x'])+tab
            lon = str(catinfo.iloc[i]['centroid_y'])+tab
            LAND_USE_CLASS = 'WATER'+tab
            VEG_CLASS = 'WATER'+tab
            SOIL_PROFILE ='SOILPROF'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab
            SLOPE = str(catslope)+tab
            ASPECT = '200'+tab
            orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
    orvh.write(":EndHRUs"+"\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return catinfo




def calculateChannaln(width,depth,Q,slope):
    zch = 2
    sidwd = zch * depth ###river side width
    tab = "          "
    botwd = width - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*width
        sidwd = 0.5*0.5*width
        zch = (width - botwd)/2/depth
    Ach = botwd*depth + 2*zch*depth*depth/2
#    arcpy.AddMessage(depth)
#    arcpy.AddMessage(zch)
#    arcpy.AddMessage(botwd)
#    arcpy.AddMessage(width)
#    arcpy.AddMessage(slope)

    Pch = botwd + 2*depth*(1+zch**2)**0.5
    Rch = float(Ach)/float(Pch)  ### in meter
    V = float(Q)/float(Ach)
    n = (Rch**(2.0/3.0))*(slope**(1.0/2.0))/V
    return n




def Writecatinfotodbf(catinfo):

    for i in range(0,len(catinfo)):
        if catinfo['SubId'].values[i] == catinfo['DowSubId'].values[i]:
            catinfo.loc[i,'DowSubId'] = -1
        
        if catinfo['BkfWidth'].values[i] < 0:           #### if no bankfulll width data avaiable for this catchment
            twidth = catinfo['BkfWidth'].values[i]      
            ccurid = catinfo['SubId'].values[i]      ### ccurid  is the current catchment id 
            isdown = -1
            while(twidth < 0 and ccurid > 0):        
                downid = catinfo[catinfo['SubId'] == ccurid]['DowSubId'].values[0]   ### get the downstream if of current catchment 
                if downid <= 0:                                   ### if no donwstream catchment exist 
                    twidth = 1.2345                               #### define a default value if not downstream exist and upstream do not have bankfull width data
                    ccurid = -1
                    isdown = -1
                else:                                              ### if down stream exist
                    if ccurid == downid:                           ### if downstream id = current catchment id;; downstream catchemnt is not exist
                        twidth = 1.2345
                        ccurid = -1
                        isdown = -1
                    else:
                        isdown = 1                                    
                        dowcatinfo = catinfo[catinfo['SubId'] ==downid]  ### get downstream id  
                        twidth = dowcatinfo['BkfWidth'].values[0]                      ### update twidth catchment with  the downstream id 
                        ccurid = dowcatinfo['SubId'].values[0]                       ### set currid with
#                arcpy.AddMessage(str(sinfo[0,0]) +"       "+ str(twidth))
            if twidth == 1.2345 and ccurid == -1:
                catinfo.loc[i,'BkfWidth'] = 1.2345
                catinfo.loc[i,'BkfDepth'] = 1.2345
                catinfo.loc[i,'Q_Mean'] = 1.2345
            elif isdown == 1:
                catinfo.loc[i,'BkfWidth'] =dowcatinfo['BkfWidth'].values[0]
                catinfo.loc[i,'BkfDepth'] =dowcatinfo['BkfDepth'].values[0]
                catinfo.loc[i,'Q_Mean'] = dowcatinfo['Q_Mean'].values[0]
                            
        catinfo.loc[i,'Ch_n'] = calculateChannaln(catinfo['BkfWidth'].values[i],catinfo['BkfDepth'].values[i],
                                           catinfo['Q_Mean'].values[i],catinfo['RivSlope'].values[i])

    return catinfo



def Checkcat2(prow,pcol,nrows,ncols,lid,lake):
    noout=0
    nearcat = np.full(8,-9)
    ### if the point  (prow,pcol) is at the boundary of catchment
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if not lake[prow-1,pcol+1] == lid:
            noout=1
            nearcat[0] = lake[prow-1,pcol+1]
        if not lake[prow-1,pcol-1] == lid:
            noout=1
            nearcat[1] = lake[prow-1,pcol-1]
        if not lake[prow-1,pcol] == lid:
            noout=1
            nearcat[2] = lake[prow-1,pcol]
        if not lake[prow,pcol+1] == lid:
            noout=1
            nearcat[3] = lake[prow,pcol+1]
        if not lake[prow,pcol-1] == lid:
            noout=1
            nearcat[4] = lake[prow,pcol-1]
        if not lake[prow+1,pcol-1] == lid:
            noout=1
            nearcat[5] = lake[prow+1,pcol-1]
        if not lake[prow+1,pcol+1] == lid:
            noout=1
            nearcat[6] = lake[prow+1,pcol+1]
        if not lake[prow+1,pcol] == lid:
            noout=1
            nearcat[7] = lake[prow+1,pcol]
    nearcat = nearcat[nearcat > 0]
    return noout,nearcat


def Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo):
    catidx = finalcat == catid
    rivids = rivlen > 0
    rivincat = np.logical_and(catidx, rivids)
    Landtypes = landuse[rivincat]
    Landtypeid = np.unique(Landtypes)
    Landtypeid1 = Landtypeid[Landtypeid >= 0]
    Landtypeid2 = Landtypeid1[Landtypeid1 > 0]
    Landtypes = Landtypes[Landtypes > 0]
    if len(Landtypes) > 0 and float(len(Landtypeid2))/float(len(Landtypeid1)) >= 0.1:
        sum = 0.0
        for i in range(0,len(Landtypeid2)):
            iid = Landtypeid2[i]
            sum = sum + landuseinfo[landuseinfo['RasterV'] == iid]['MannV'].values*len(np.argwhere(Landtypes == iid))
        floodn = sum/(len(Landtypes))
    else:
        Landtypes2 = landuse[catidx]
        Landtypes2 = Landtypes2[Landtypes2>0]
        Landtypeiicat = np.unique(Landtypes2)
        Landtypeiicat = Landtypeiicat[Landtypeiicat > 0]
        if len(Landtypes2) > 0:
            sum = 0.0
            for i in range(0,len(Landtypeiicat)):
                iid = Landtypeiicat[i]
                sum = sum + landuseinfo[landuseinfo['RasterV'] == iid]['MannV'].values*len(np.argwhere(Landtypes2 == iid))
            floodn = sum/(len(Landtypes2))
        else:
            floodn = 0.035
#    arcpy.AddMessage(floodn)
    return float(floodn)

def Getcatwd(catid,finalcat,width,depth,Q_Mean,DA,rivpath):
    catregs = finalcat == catid
    riverp = rivpath > 0
    rivincat = np.logical_and(catregs, riverp)
    wd = width[rivincat]
    dp = depth[rivincat]
    Q = Q_Mean[rivincat]
    wd = wd[wd > 0]
    dp = dp[dp > 0]
    Q =  Q[Q > 0]
    if len(wd) > 0:
        unique, counts = np.unique(wd, return_counts=True)
        catwd = np.average(unique, weights=counts)
        unique, counts = np.unique(dp, return_counts=True)
        catdps = np.average(unique, weights=counts)
        unique, counts =  np.unique(Q, return_counts=True)
        catQ = np.average(unique, weights=counts)
    else:
        catwd = -9
        catdps = -9
        catQ = -9
    return catwd,catdps,catQ
############################################################



def Getcatrivlenslope(catrow,catcol,rivlen,dem,fac,hydir,finalcat,trow,tcol,nrows,ncols,slope,rivpath):
    rivs = rivlen[catrow,catcol]
    rivs = np.unique(rivs)  ### get all river length within the catchment 
    rivs = rivs[rivs > 0]
    rivcinfo = np.full((len(catrow),4),-999999999999.99999)
    rivcinfo[:,0] = rivlen[catrow,catcol]   #### store river length 
    rivcinfo[:,1] = fac[catrow,catcol]   ### flow accumulation 
    rivcinfo[:,2] = catrow    #### row 
    rivcinfo[:,3] = catcol    ### col 
    rivout =  np.full((len(rivs),4),-999999999999.000000) ### store riv info for each river path started at boundary of catchment 
    rivout2 =  np.full((len(rivs),4),-999999999999.000000) ### store riv info for each river path started at inside of catchment 
    for i in range(0,len(rivs)):  ## loop for each river segments within the catchment  each river segment was labled with their river legth 
        rivsid = rivs[i] 
        rivcinfo2 = rivcinfo[rivcinfo[:,0]==rivsid,]   ### get information of the river segment i
        rivcinfo2 = rivcinfo2[rivcinfo2[:,1].argsort()]
        prow = rivcinfo2[0,2].astype(int)   ## find the first grids of the river segment 
        pcol = rivcinfo2[0,3].astype(int)   ## find the first grids of the river segment  
        lid = finalcat[prow,pcol]     #### catid of start point of stream cell 
        nout, nearcat = Checkcat2(prow,pcol,nrows,ncols,lid,finalcat)  #### check if the point (prow,pcol) is close to the catchment boundary, most upstream cell 
        rivtemp = np.full((len(catrow),4),-9999999999.999999)
        icell = 0
        iscase1 = 0
        if len(nearcat) > 0:  ### check if one of the near cat is the upstream catchment
            for incat in range(0,len(nearcat)):
                inearcat = nearcat[incat]
                incat_trow,incat_tcol = Getbasinoutlet(inearcat,finalcat,fac,hydir,nrows,ncols)
                incat_nrow,incat_ncol = Nextcell(hydir,incat_trow,incat_tcol)### get the downstream catchment id
                if incat_nrow >= nrows or incat_nrow < 0 or incat_ncol >= ncols or incat_ncol < 0:
                    continue
                if finalcat[incat_nrow,incat_ncol] == lid:
                    iscase1 = 1
        if nout > 0 and iscase1 == 1:  #### this means the strem connect at other catchments. or  this catchment only have one stream some head stream was included
#            arcpy.AddMessage("in riv  case 1   " + str(rivsid) + "     " + str(iscase1))
            nrow = prow
            ncol = pcol
            rivpath[nrow,ncol] = 1
            while finalcat[nrow,ncol] == finalcat[trow,tcol]:
                flen_orow,flen_ocol = nrow,ncol
                if flen_orow < 0 or flen_ocol<0 or icell >= len(rivtemp):
                    break
                rivpath[nrow,ncol] = 1
                rivtemp[icell,0] = rivlen[nrow,ncol]
                rivtemp[icell,1] = dem[nrow,ncol]
                rivtemp[icell,3] = slope[nrow,ncol]
                if icell > 0:
                    if rivtemp[icell,0] != rivtemp[icell - 1,0]:
                        rivtemp[icell,2] = rivtemp[icell,0] + rivtemp[icell - 1,2]
                    else:
                        rivtemp[icell,2] = rivtemp[icell-1,2]
                else:
                    rivtemp[icell,2] = rivtemp[icell,0]
                icell = icell + 1
                nrow,ncol = Nextcell(hydir,int(flen_orow),int(flen_ocol))
                if nrow < 0 or ncol < 0:
                    nrow,ncol = Nextcell(hydir,int(trow),int(tcol))
                    print("warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol)
                    nrow,ncol = 0,0
                if nrow >= nrows or ncol >= ncols:
#                    arcpy.AddMessage("out of the boundary")
                    break                

            rivtemp = rivtemp[rivtemp[:,0]>0,]
            if icell > 0:
                icell = min(icell,len(rivtemp))
                rivout[i,0] = rivtemp[icell-1,2]
                rivout[i,2] = float(max(abs(rivtemp[0,1] - rivtemp[icell-1,1]),0.5))/float(rivtemp[icell-1,2])
#                arcpy.AddMessage("rivsid     case1  " + str(rivtemp[0,1]) + "       end dem " + str(rivtemp[icell-1,1]) + '     river len    ' + str(rivtemp[icell-1,2]))
#                rivout[i,2] = (max(rivtemp[:,1]) - min(rivtemp[:,1]))/rivtemp[icell-1,2]
                rivtemp = rivtemp[rivtemp[:,3]>=0,]
                if len(rivtemp) > 0:
                    rivout[i,1] = np.mean(rivtemp[:,3])
                else:
                    rivout[i,1] = -9999
        rivtemp2 = np.full((len(catrow),4),-9999999999.999999)
        icell = 0
        if len(rivs) > 0 and iscase1 == 0: # for river segment start inside the catchemnt and do not rieve water from upstream catchment 
#            arcpy.AddMessage("in riv  case 2   " + str(rivsid))
            nrow = prow
            ncol = pcol
            rivpath[nrow,ncol] = 1
            while finalcat[nrow,ncol] == finalcat[trow,tcol]:
                flen_orow,flen_ocol = nrow,ncol
                if flen_orow < 0 or flen_ocol<0 or icell >= len(rivtemp2):
                    break
                rivpath[nrow,ncol] = 1
                rivtemp2[icell,0] = rivlen[nrow,ncol]    #### store riv length of each channel cell 
                rivtemp2[icell,1] = dem[nrow,ncol]       #### store riv dem of each channel cell 
                rivtemp2[icell,3] = slope[nrow,ncol]     #### store slope  of each channel cell 
                if icell > 0: ### start from the second cell 
                    if rivtemp2[icell,0] != rivtemp2[icell - 1,0]:   ### come to a new river segment 
                        rivtemp2[icell,2] = rivtemp2[icell,0] + rivtemp2[icell - 1,2]    ### store cumulated river length of different river segment
                    else:
                        rivtemp2[icell,2] = rivtemp2[icell-1,2]   ### still in old river segment, do not update cumulated river length 
                else:
                    rivtemp2[icell,2] = rivtemp2[icell,0]  ### stroe the river length of the first river segment 
                icell = icell + 1   ## move to next cell 
                nrow,ncol = Nextcell(hydir,int(flen_orow),int(flen_ocol))
                if nrow < 0 or ncol < 0:  ### if the next cell move out of the domain 
                    nrow,ncol = Nextcell(hydir,int(trow),int(tcol))
                    print("warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol)
                    nrow,ncol = 0,0
                if nrow >= nrows or ncol >= ncols: ### if the next cell move out of the domain 
#                    arcpy.AddMessage("out of the boundary")
                    break 
            rivtemp2 = rivtemp2[rivtemp2[:,0]>0,]    ### remove the emmpty positions in the array 
            if icell > 0:   ## means there is river grids  icell should equal with number of river grids 
                icell = min(icell,len(rivtemp2))   
                rivout2[i,0] = rivtemp2[icell-1,2]   ####  get the acuumulated river length for river i
                rivout2[i,2] = float(max(abs(rivtemp2[0,1] - rivtemp2[icell-1,1]),0.5))/float(rivtemp2[icell-1,2])  ### (dem_begin - dem_end)/river length, the difference minimum value is 0.5 m
#                arcpy.AddMessage("rivsid   case 2   " + str(rivtemp2[0,1]) + "       end dem " + str(rivtemp2[icell-1,1]) + '     river len    ' + str(rivtemp2[icell-1,2]))
#                rivout2[i,2] = (max(rivtemp2[:,1]) - min(rivtemp2[:,1]))/rivtemp2[icell-1,2] ### ####  get the river slope  for river i 
                rivtemp2 = rivtemp2[rivtemp2[:,3]>=0,]    
                if len(rivtemp2) > 0:
                    rivout2[i,1] = np.mean(rivtemp2[:,3])  ####  get the river slope  for river average 
                else:
                    rivout2[i,1] = -9999
    rivout = rivout[rivout[:,0]>0,]
    rivout2 = rivout2[rivout2[:,0]>0,]
    if len(rivout) > 0:
        rivout = rivout[rivout[:,0].argsort()]   ### sort with river length of each river segment, the river slopes from longest river segment was used
        outrivlen = rivout[len(rivout)-1,0]
        outrivslp =  rivout[len(rivout)-1,2] ### slope np.mean(slope along the river channel )
        outrivslp2 = rivout[len(rivout)-1,1]   ### slope max (den_b - dem end , 0.5)/ river length 
#        arcpy.AddMessage('final slope ' + str(outrivslp) + "     " +  str(outrivslp2) + "    " + str(outrivlen))
    elif len(rivout2) > 0:
        rivout2 = rivout2[rivout2[:,0].argsort()]
        outrivlen = rivout2[len(rivout2)-1,0]
        outrivslp =  rivout2[len(rivout2)-1,2]
        outrivslp2 = rivout2[len(rivout2)-1,1]
#        arcpy.AddMessage('final slope ' + str(outrivslp) + "     " +  str(outrivslp2))
    else:
        outrivlen = -9999.00
        outrivslp = -9999.00
        outrivslp2 = -9999.00
    return outrivlen, outrivslp,outrivslp2,rivpath
######################################################


def Generatecatinfo(Watseds,fac,fdir,lake,dem,area,catinfo,allcatid,lakeinfo,width,depth,
                    rivlen,obs,nrows,ncols,slope,landuse,landuseinfo,Q_Mean):
    finalcat = copy.copy(Watseds)
    rivpath = copy.copy(Watseds)
    rivpath[:,:] = -9999
    for i in range(0,len(allcatid)):
#        print("subid      " + str(allcatid[i].astype(int)))
        catid = allcatid[i].astype(int)
        catinfo.loc[i,'SubId'] = catid
        rowcol = np.argwhere(finalcat==catid).astype(int)
        trow,tcol = Getbasinoutlet(catid,finalcat,fac,fdir,nrows,ncols)
        nrow,ncol = Nextcell(fdir,trow,tcol)### get the downstream catchment id
        if nrow < 0 or ncol < 0:
            catinfo.loc[i,'DowSubId'] = -1
        elif nrow >= nrows or ncol >= ncols:
            catinfo.loc[i,'DowSubId'] = -1
        elif finalcat[nrow,ncol] <= 0:
            catinfo.loc[i,'DowSubId'] = -1
        else:
            catinfo.loc[i,'DowSubId'] = finalcat[nrow,ncol]
#        catinfo[i,2] = trow
#        catinfo[i,3] = tcol
################################## Get lake information
        lakeid = lake[trow,tcol]
        if lakeid > 0:
            slakeinfo = lakeinfo.loc[lakeinfo['Hylak_id'] == lakeid]
            catinfo.loc[i,'IsLake'] = 1
            catinfo.loc[i,'HyLakeId'] = lakeid
            catinfo.loc[i,'LakeVol'] = slakeinfo.iloc[0]['Vol_total']
            catinfo.loc[i,'LakeArea']= slakeinfo.iloc[0]['Lake_area']
            catinfo.loc[i,'LakeDepth']= slakeinfo.iloc[0]['Depth_avg']
            catinfo.loc[i,'Laketype'] = slakeinfo.iloc[0]['Lake_type']
#            catinfo[i,29] = min(float(len(lake[lake == lakeid]))/float(len(finalcat[finalcat == catid])),1.0)
########Check if it is observation points
        if obs[trow,tcol]  >= 0:
#            arcpy.AddMessage(str(catid)+"      "+str(obs[trow,tcol]))
            catinfo.loc[i,'IsObs'] =  obs[trow,tcol]
########Got basin width and depth
        catrivlen,catrivslp,catrivslp2,rivpath = Getcatrivlenslope(rowcol[:,0],rowcol[:,1],rivlen,dem,fac,fdir,finalcat,
                                                trow,tcol,nrows,ncols,slope,rivpath)
        catwidth,catdepth,catQ = Getcatwd(catid,finalcat,width,depth,Q_Mean,-1,rivlen) ### width depth in m
#        arcpy.AddMessage("catid is    " + str(catid) + "    " + str(catwidth))
        catinfo.loc[i,'MeanElev'] = float(sum(dem[rowcol[:,0],rowcol[:,1]])/float(len(rowcol))) ### average elevation
#        catinfo[i,13] = float(sum(area[rowcol[:,0],rowcol[:,1]]))/1000/1000  #### maximum area in km^2
#        catinfo[i,14] = max(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
#        catinfo[i,15] = min(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
#        catinfo[i,16] = dem[trow,tcol] #### outlet elevation
        catinfo.loc[i,'BkfWidth'] = catwidth
        catinfo.loc[i,'BkfDepth'] = catdepth
#        catinfo[i,19] = 0.035
#######Got basin area and rivlen
#        catinfo[i,11] = np.mean(area[rowcol[:,0],rowcol[:,1]])
        catinfo.loc[i,'Rivlen'] = catrivlen
        slopet = slope[rowcol[:,0],rowcol[:,1]]
        slopet = slopet[slopet>0,]
        catinfo.loc[i,'BasinSlope'] = np.mean(slopet)
        if catrivslp < 0:
            catinfo.loc[i,'RivSlope'] = np.mean(slopet)
        else:
            catinfo.loc[i,'RivSlope'] = catrivslp
#        catinfo[i,27] = catrivslp2
        if len(slopet) < 1:
            catinfo.loc[i,'RivSlope'] = 0.001
            catinfo.loc[i,'BasinSlope'] = 0.001
        catinfo.loc[i,'FloodP_n'] = Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo)
        catinfo.loc[i,'Q_Mean'] = catQ
#    writeraster(outputFolder + "/" + "rivpath.asc",rivpath,OutputFolder + "/" + "dir")
    return catinfo,rivpath




def checklakeboundary(lake1,lid,p_row,p_col):
    numnonlake = 0
    nonlakedir = np.full(10,-999)
    if lake1[p_row + 0,p_col + 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 0
    if lake1[p_row + 1,p_col + 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 1
    if lake1[p_row + 1,p_col + 0] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 2
    if lake1[p_row + 1,p_col - 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 3
    if lake1[p_row + 0,p_col - 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 4
    if lake1[p_row - 1,p_col - 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 5
    if lake1[p_row - 1,p_col + 0] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 6
    if lake1[p_row - 1,p_col + 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 7
    return numnonlake,nonlakedir
    


########## Return row and column number of a specific catchment's outlet  
# hylake       : a dataframe from HydroLAKE database containing all infromation about lakes 
# noncnlake    : a two dimension array of lakes that are not connected by river network 
# NonConLThres : the flow accumulation 2D array  
# dir    : the flow direction 2D array
# nrows, ncols  : maximum number of rows and columns in the basin 2D array   
# return crow,ccol  : row and column id of the catchment outlet.

def Dirpoints2(N_dir,p_row,p_col,lake1,lid,goodpoint,k,ncols,nrows):
    ndir = copy.copy(N_dir)
    ip = copy.copy(k) + 1
#dir 1
    if lake1[p_row + 0,p_col + 1] == lid:  #### it is a lake cell
        if p_row + 0 == nrows-1 or p_col + 1 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 0,p_col + 1)
        if numnonlake != 0: # it is a boundary lake cell
            tt = goodpoint[goodpoint[:,0] == p_row + 0,]
            if len(tt[tt[:,1] == p_col + 1,]) < 1:#### the point not exist in good points which store all boundary points
                ndir[p_row + 0,p_col + 1] = 16   ### change the flow direction of new point to old points
                goodpoint[ip,0] = p_row + 0
                goodpoint[ip,1] = p_col + 1
                ip = ip + 1
#dir 2
    if lake1[p_row + 1,p_col + 1] == lid:
        if p_row + 1 == nrows-1 or p_col + 1 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 1,p_col + 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 1,]
            if len(tt[tt[:,1] == p_col + 1,]) < 1:
                ndir[p_row + 1,p_col + 1] = 32
                goodpoint[ip,0] = p_row + 1
                goodpoint[ip,1] = p_col + 1
                ip = ip + 1
#dir 3
    if lake1[p_row + 1,p_col + 0] == lid:
        if p_row + 1 == nrows-1 or p_col + 0 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 1,p_col + 0)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 1,]
            if len(tt[tt[:,1] == p_col + 0,]) < 1:
                ndir[p_row + 1,p_col + 0] = 64
                goodpoint[ip,0] = p_row + 1
                goodpoint[ip,1] = p_col + 0
                ip = ip + 1
#dir 4
    if lake1[p_row + 1,p_col - 1] == lid:
        if p_row + 1 == nrows-1 or p_col - 1 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 1,p_col - 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 1,]
            if len(tt[tt[:,1] == p_col - 1,]) < 1:
                ndir[p_row + 1,p_col - 1] = 128
                goodpoint[ip,0] = p_row + 1
                goodpoint[ip,1] = p_col - 1
                ip = ip + 1
#dir 5
    if lake1[p_row + 0,p_col - 1] == lid:
        if p_row + 0 == nrows-1 or p_col - 1 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 0,p_col - 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 0,]
            if len(tt[tt[:,1] == p_col - 1,]) < 1:
                ndir[p_row + 0,p_col - 1] = 1
                goodpoint[ip,0] = p_row + 0
                goodpoint[ip,1] = p_col - 1
                ip = ip + 1
#dir 6
    if lake1[p_row - 1,p_col - 1] == lid:
        if p_row - 1 == nrows-1 or p_col - 1 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row - 1,p_col - 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row - 1,]
            if len(tt[tt[:,1] == p_col - 1,]) < 1:
                ndir[p_row - 1,p_col - 1] = 2
                goodpoint[ip,0] = p_row - 1
                goodpoint[ip,1] = p_col - 1
                ip = ip + 1
#dir 7
    if lake1[p_row - 1,p_col + 0] == lid:
        if p_row - 1 == nrows-1 or p_col + 0 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row - 1,p_col + 0)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row - 1,]
            if len(tt[tt[:,1] == p_col + 0,]) < 1:
                ndir[p_row - 1,p_col + 0] = 4
                goodpoint[ip,0] = p_row - 1
                goodpoint[ip,1] = p_col + 0
                ip = ip + 1
#dir 8
    if lake1[p_row - 1,p_col + 1] == lid:
        if p_row - 1 == nrows-1 or p_col + 1 == ncols - 1:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row - 1,p_col + 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row - 1,]
            if len(tt[tt[:,1] == p_col + 1,]) < 1:
                ndir[p_row - 1,p_col + 1] = 8
                goodpoint[ip,0] = p_row - 1
                goodpoint[ip,1] = p_col + 1
                ip = ip + 1
#    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir,goodpoint,ip

def ChangeDIR(dir,lake1,acc,ncols,nrows,outlakeids,nlakegrids):
    ndir = copy.copy(dir)
    for i in range(0,len(outlakeids)):
        lid = outlakeids[i]
        goodpoint = np.full((20000,2),-99999)
        lrowcol = np.argwhere(lake1==lid).astype(int)
        if len(lrowcol) > nlakegrids:
            continue
#        arcpy.AddMessage(str(lid) + "    " +  str(len(lrowcol)) + "     " + str(i))
        print("start modify lake flow direction, the lake id is    " + str(int(lid)))
        prow,pcol = Getbasinoutlet(lid,lake1,acc,dir,nrows,ncols)
        goodpoint[0,0] = prow
        goodpoint[0,1] = pcol
        if prow >= nrows - 1 or pcol == ncols  - 1:
            continue
        ip = 0
        k = 0
        ipo = -1
        while ip > ipo:
            for i in range(0,len(goodpoint[goodpoint[:,0]>0,])):
                if i > ipo:
#                    arcpy.AddMessage("start of checking:    " + str(ip) + "     "+ str(ipo)+ "   ")
#                    arcpy.AddMessage(goodpoint[0:ip,0])
                    trow = goodpoint[i,0]
                    tcol = goodpoint[i,1]
                    if trow >= nrows - 1 or tcol == ncols  - 1:
                        continue
                    ndir,goodpoint,k1= Dirpoints2(ndir,trow,tcol,lake1,lid,goodpoint,k,ncols,nrows)
#                    arcpy.AddMessage("start of checking:    " + str(k1) + "     "+ str(len(goodpoint[goodpoint[:,0]>0,]))+ "   ")
                    k = k1 - 1
            ipo = ip
            ip = len(goodpoint[goodpoint[:,0]>0,]) - 1
            k = ip
#            arcpy.AddMessage("start of checking:    " + str(ip) + "     "+ str(ipo)+ "   " + str(len(lrowcol)))
#        goodpoint[goodpoint[:,0]>0,]
#        ndir[goodpoint[goodpoint[:,0]>0,0].astype(int),goodpoint[goodpoint[:,0]>0,1].astype(int)] = 9
    return ndir

  
def Checklake(prow,pcol,nrows,ncols,lid,lake):
    noout=0
    ### double check if the head stream cell is nearby the lake, if near by the lake the stream was ignored
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if lake[prow-1,pcol+1] == lid:
            noout=1
        if lake[prow-1,pcol-1] == lid:
            noout=1
        if lake[prow-1,pcol] == lid:
            noout=1
        if lake[prow,pcol+1] == lid:
            noout=1
        if lake[prow,pcol-1] == lid:
            noout=1
        if lake[prow+1,pcol-1] == lid:
            noout=1
        if lake[prow+1,pcol+1] == lid:
            noout=1
        if lake[prow+1,pcol] == lid:
            noout=1
    return noout
############################################################################33

####################################################################3
def Checkcat(prow,pcol,nrows,ncols,lid,lake):
    noout=0
    ### double check if the head stream cell is nearby the lake, if near by the lake the stream was ignored
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if not lake[prow-1,pcol+1] == lid:
            noout=1 + noout
        if not lake[prow-1,pcol-1] == lid:
            noout=1 + noout
        if not lake[prow-1,pcol] == lid:
            noout=1 + noout
        if not lake[prow,pcol+1] == lid:
            noout=1 + noout
        if not lake[prow,pcol-1] == lid:
            noout=1 + noout
        if not lake[prow+1,pcol-1] == lid:
            noout=1 + noout
        if not lake[prow+1,pcol+1] == lid:
            noout=1 + noout
        if not lake[prow+1,pcol] == lid:
            noout=1 + noout
    return noout

###################################################################3

def selectlake(hylake,noncnlake,NonConLThres,hylakeinfo):
#    arcpy.AddMessage("123asdfasdfasdfasd")
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(noncnlake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(noncnlake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['Hylak_id'] == sl_lid]
        if len(slakeinfo) <=0:
            continue
        if slakeinfo.iloc[0]['Lake_area'] >= NonConLThres:
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = sl_lid
    return sl_lake

def selectlake2(hylake,Lakehres,hylakeinfo):
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(sl_lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(sl_lake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['Hylak_id'] == sl_lid]
        if len(slakeinfo)<=0:
#            print("Lake excluded     asdfasd " + str(sl_lid))
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999
            continue
        if slakeinfo.iloc[0]['Lake_area'] < Lakehres:
#            print("Lake excluded     due to area " + str(sl_lid))
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999
    return sl_lake

##################################################################3

def Addobspoints(obs,pourpoints,boid,cat):
    obsids = np.unique(obs)
    obsids = obsids[obsids>=0]
    for i in range(0,len(obsids)):
        rowcol = np.argwhere(obs==obsids[i]).astype(int)
        if cat[rowcol[0,0],rowcol[0,1]] > 0:
            if pourpoints[rowcol[0,0],rowcol[0,1]] < 0:
                pourpoints[rowcol[0,0],rowcol[0,1]] = boid + obsids[i]
    return pourpoints

####################################################################
def CE_mcat4lake(cat1,lake,fac,fdir,bsid,nrows,ncols,Pourpoints):
    #####adjust for some lakes was divided into two catchment beacuse of flow direction and in stream. double check each lake
    ##### and decide if need to merge these catchment into one lake catchment.
    cat = copy.copy(cat1)
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    outlakeids = np.full(1000000,-99999)
    outi = 0
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        lakacc = np.full((len(lrowcol),3),-9999)
        lakacc[:,0] = lrowcol[:,0]
        lakacc[:,1] = lrowcol[:,1]
        lakacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakacc = lakacc[lakacc[:,2].argsort()]
        lorow = lakacc[len(lakacc)-1,0]
        locol = lakacc[len(lakacc)-1,1]  ###### lake outlet row and col
        arclakeid = cat[lorow,locol]  ####### lake catchment id 
        lakecatrowcol = np.argwhere(cat==arclakeid).astype(int)    #### lake catchment cells  
        lakenrow,lakencol = Nextcell(fdir,lorow,locol)
        if lakenrow < 0 or lakencol < 0:
            lakedowcatid = -999
        else:
            if lakenrow >= nrows or lakencol >= ncols:
                continue
            lakedowcatid = cat[lakenrow,lakencol]
        # if not arclakeid < bsid and arclakeid > blid:
        #     continue
        arcatid,catcounts = np.unique(cat[lrowcol[:,0],lrowcol[:,1]],return_counts=True) ###### all catchment id containing this lake
        tarid = 0
        ### if there are more than 1 catchment in cat1, determine if they need to be combined
        ### check if these catchment flow into the lake if it is true, change catchment id into lake catchment id
        if len(arcatid)>1:  #
#            if float(len(lakecatrowcol))/float(len(lrowcol)) < 0.9: #and len(lrowcol) < 10000: #: float(max(catcounts))/float(len(lrowcol)) < 0.8 and 
            outlakeids[outi] = lakeid
            outi = outi + 1
            for j in range(0,len(arcatid)):
                crowcol = np.argwhere(cat==arcatid[j]).astype(int)
                catacc = np.full((len(crowcol),3),-9999)
                catacc[:,0] = crowcol[:,0]
                catacc[:,1] = crowcol[:,1]
                catacc[:,2] = fac[crowcol[:,0],crowcol[:,1]]
                catacc = catacc[catacc[:,2].argsort()]
                catorow = catacc[len(catacc)-1,0]
                catocol = catacc[len(catacc)-1,1] ### catchment outlet
                Lakeincat = lake[crowcol[:,0],crowcol[:,1]]
                nlake = np.argwhere(Lakeincat==lakeid).astype(int)
                nrow,ncol = Nextcell(fdir,catorow,catocol) #####Get the next row and col of downstream catchment
                if nrow < 0 or ncol < 0:
                    continue
                if nrow < nrows and ncol < ncols:
               ### if downstream catchment is target lake,and this catchment is an lakeinflow catchment combine them
                    if cat[nrow,ncol] == arclakeid and float(len(nlake))/float(len(crowcol)) > 0.1 and cat[catorow,catocol] > bsid:
                        cat[crowcol[:,0],crowcol[:,1]] = arclakeid
                    if float(len(nlake))/float(len(lrowcol)) > 0.1 and cat[catorow,catocol] > bsid and cat[catorow,catocol] != lakedowcatid:
#                        arcpy.AddMessage("2")
                        cat[crowcol[:,0],crowcol[:,1]] = arclakeid
#                        if cat[catorow,catocol] != arclakeid and cat[nrow,ncol] != arclakeid:
#                            print lakeid
                    if cat[nrow,ncol] > bsid and arcatid[j] > bsid:  #### lake input cat route to another lake input catch
                        cat[crowcol[:,0],crowcol[:,1]] = cat[nrow,ncol]
        pp = Pourpoints[lrowcol[:,0],lrowcol[:,1]]
        pp = np.unique(pp)
        pp = pp[pp > 0]
        if len(pp) == 1:
            cat[lrowcol[:,0],lrowcol[:,1]] = arclakeid
    outlakeids= outlakeids[outlakeids > 0]
    return cat,outlakeids
###################################################33
def CE_mcat4lake2(cat1,lake,fac,fdir,bsid,nrows,ncols,Pourpoints):
    cat = copy.copy(cat1)
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        lakacc = np.full((len(lrowcol),3),-9999)
        lakacc[:,0] = lrowcol[:,0]
        lakacc[:,1] = lrowcol[:,1]
        lakacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakacc = lakacc[lakacc[:,2].argsort()]
        lorow = lakacc[len(lakacc)-1,0]
        locol = lakacc[len(lakacc)-1,1]  ###### lake outlet row and col
        arclakeid = cat1[lorow,locol]
        pp = Pourpoints[lorow,locol]
        pp = np.unique(pp)
        pp = pp[pp > 0]
        if len(pp) == 1:
            if arclakeid < 0:
                cat[lrowcol[:,0],lrowcol[:,1]] = pp
            else:
                cat[lrowcol[:,0],lrowcol[:,1]] = arclakeid
#### For some reason pour point was missing in non-contribute catchment 
    Pours = np.unique(Pourpoints)
    Pours = Pours[Pours>0]
    for i in range(0,len(Pours)):
        pourid = Pours[i]
        rowcol = Pourpoints == pourid
        if cat[rowcol] < 0:
            nout = Checkcat(rowcol[0,0],rowcol[0,1],nrows,ncols,pourid,cat)
            if len(cat[cat == pourid]) > 0 and nout < 8:
                cat[rowcol] = pourid
    rowcol1 = fac > 0
    rowcol2 = cat < 0
    noncontribuite = np.logical_and(rowcol1, rowcol2)
    cat[noncontribuite] = 2*max(np.unique(cat)) + 1
    return cat
######################################################
def CE_Lakeerror(fac,fdir,lake,cat2,bsid,blid,boid,nrows,ncols,cat):
    Watseds = copy.copy(cat2)
    Poups = np.unique(Watseds)
    Poups = Poups[Poups>=0]
    ##### Part 2, remove some small catchment which is not lake catchment
    out = np.full((len(Poups),4),-9999)
    for i in range(0,len(Poups)):
        catid = Poups[i]
        if catid > boid:
            continue #### do nothing for observation catchments
        rowcol = np.argwhere(Watseds==catid).astype(int)
        catacc = np.full((len(rowcol),3),-9999)
        catacc[:,0] = rowcol[:,0]
        catacc[:,1] = rowcol[:,1]
        catacc[:,2] = fac[rowcol[:,0],rowcol[:,1]]
        catacc = catacc[catacc[:,2].argsort()].astype(int)
        rowcol[0,0] = catacc[len(catacc)-1,0]
        rowcol[0,1] = catacc[len(catacc)-1,1]
        nrow,ncol = Nextcell(fdir,rowcol[0,0],rowcol[0,1])### get the downstream catchment id
        if nrow < 0 or ncol < 0:
            continue
        if nrow < nrows and ncol < ncols:
            if len(rowcol) < 10 and Watseds[rowcol[0,0],rowcol[0,1]] > bsid:
                Watseds[catacc[:,0],catacc[:,1]] = Watseds[nrow,ncol]
            if len(rowcol) < 10 and Watseds[rowcol[0,0],rowcol[0,1]] < blid:
                Watseds[catacc[:,0],catacc[:,1]] = Watseds[nrow,ncol]
    return Watseds

#########################################33
def GenerateFinalPourpoints(fac,fdir,lake,cat3,bsid,blid,boid,nrows,ncols,cat,obs):
    Poups = copy.copy(cat3)
    Poups[:,:]=-9999
    GWat = copy.copy(cat3)
    GWatids = np.unique(cat3)
    GWatids = GWatids[GWatids>=0]
    ncatid = 1
    for i in range(0,len(GWatids)):
        trow,tcol = Getbasinoutlet(GWatids[i],GWat,fac,fdir,nrows,ncols)
#        arcpy.AddMessage("result     " +str(trow) +"    " +str(tcol))
        Poups[trow,tcol] = ncatid
        ncatid = ncatid + 1
    OWat = copy.copy(cat)
    OWatids = np.unique(cat)
    OWatids = OWatids[OWatids>=0]
    for i in range(0,len(OWatids)):
        trow,tcol = Getbasinoutlet(OWatids[i],OWat,fac,fdir,nrows,ncols)
        if not GWat[trow,tcol] >= blid:
            if Poups[trow,tcol] < 0:
                Poups[trow,tcol] = ncatid
                ncatid = ncatid + 1
    obsids = np.unique(obs)
    obsids = obsids[obsids>=0]
    for i in range(0,len(obsids)):
        rowcol = np.argwhere(obs==obsids[i]).astype(int)
        if Poups[rowcol[0,0],rowcol[0,1]] < 0 and lake[rowcol[0,0],rowcol[0,1]] < 0:
            Poups[rowcol[0,0],rowcol[0,1]] = ncatid
            ncatid = ncatid + 1
    return Poups
#######
####
# ####################################################33


def GenerPourpoint(cat,lake,Str,nrows,ncols,blid,bsid,bcid,fac,hydir):
    GP_cat = copy.copy(cat)
    sblid = copy.copy(blid)
    ############### Part 1 Get all pourpoints of hydroshed catchment
    arcatid = np.unique(cat)#### cat all catchment idd
    arcatid = arcatid[arcatid>=0]
    catoutloc = np.full((len(arcatid),3),-9999)
    for i in range(0,len(arcatid)):
        catid = arcatid[i]
        catrowcol = np.argwhere(cat==catid).astype(int)
        trow,tcol = Getbasinoutlet(catid,cat,fac,hydir,nrows,ncols)
        GP_cat[catrowcol[:,0],catrowcol[:,1]]=-9999   ### set the catment cells into null
        GP_cat[trow,tcol]=bcid #### change the outlet of catchment into wid
        bcid = bcid + 1
        catoutloc[i,0] = catid  ## catchment id
        catoutloc[i,1] = trow  #### catchment pourpont row
        catoutloc[i,2] = tcol  #### catchment pourpont col
#    writeraster(outFolder+subid+"_Pourpoints_1.asc",GP_cat)
    ##################Part 2 Get pourpoints of Lake inflow streams
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)): #### loop for each lake
        lid = arlakeid[i]
        rowcol = np.argwhere(lake==lid.astype(int))
        nrow = rowcol.shape[0]
        Stridinlake = np.full(nrow,-9999)
        Stridinlake[:] = Str[rowcol[:,0],rowcol[:,1]]
        Strid_L  = np.unique(Stridinlake[np.argwhere(Stridinlake > 0).astype(int)]) ### Get all stream that intercept with lake
        ##### find the intercept point of stream and lake
        for j in range(0,len(Strid_L)):  #### loop for each stream intercept with lake
            strid = Strid_L[j]
            strrowcol = np.argwhere(Str == strid).astype(int)
            nstrrow = strrowcol.shape[0]
            Strchek = np.full((nstrrow,4),-9999)##### 0 row, 1 col, 2 fac,
            Strchek[:,0] = strrowcol[:,0]
            Strchek[:,1] = strrowcol[:,1]
            Strchek[:,2] = fac[strrowcol[:,0],strrowcol[:,1]]
            Strchek = Strchek[Strchek[:,2].argsort()].astype(int)
            ibg=-99
            for irowst in range(nstrrow):#### search from smallest acc stream cell
                if lake[Strchek[irowst,0],Strchek[irowst,1]] == lid and ibg==-99: ### if the begining of stream in lake the stream is ingnored
                    ibg = 1
                    if irowst != 0:
                        if lake[Strchek[irowst-1,0],Strchek[irowst-1,1]] == -9999:
                        ##### this means the stream connect two lakes, so must assign an pourpoints
                            if len(np.unique(lake[Strchek[:,0],Strchek[:,1]])) >= 3:
                                GP_cat[Strchek[irowst-1,0],Strchek[irowst-1,1]] = bsid
                                bsid = bsid + 1
                        ### double check if the head stream cell is nearby the lake, if near by the lake the stream was ignored
                        if Strchek[0,0] != 0 and Strchek[0,0] != nrows -1 and Strchek[0,1] != 0 and Strchek[0,1] != ncols-1:
                            noout = Checklake(Strchek[0,0],Strchek[0,1],nrows,ncols,lid,lake)
                        ##### the head stream celll is not nearby the lake
                        if noout == 0:
                            GP_cat[Strchek[irowst-1,0],Strchek[irowst-1,1]] = bsid
                            bsid = bsid + 1
                        #### it is possible that two steam combine together near the lake, double check if the stream conncet to
                        # anotehr stream and this steam is not witin the lake
                    if irowst == 0 or noout == 1:
                        nostr = Str[Strchek[0,0],Strchek[0,1]]
                        a = 0
                        orowcol = np.full((8,3),-9999)
                        if Strchek[0,0] != 0 and Strchek[0,0] != nrows -1 and Strchek[0,1] != 0 and Strchek[0,1] != ncols-1:
                            if Str[Strchek[0,0]-1,Strchek[0,1]+1] != -9999 and lake[Strchek[0,0]-1,Strchek[0,1]+1] == -9999:
                                orowcol[a,0] = Strchek[0,0]-1
                                orowcol[a,1] = Strchek[0,1]+1
                                orowcol[a,2] = Str[Strchek[0,0]-1,Strchek[0,1]+1]
                                a = a+1
                            if Str[Strchek[0,0]-1,Strchek[0,1]-1] != -9999 and lake[Strchek[0,0]-1,Strchek[0,1]-1] == -9999:
                                orowcol[a,0] = Strchek[0,0]-1
                                orowcol[a,1] = Strchek[0,1]-1
                                orowcol[a,2] = Str[Strchek[0,0]-1,Strchek[0,1]-1]
                                a = a+1
                            if Str[Strchek[0,0]-1,Strchek[0,1]] != -9999 and lake[Strchek[0,0]-1,Strchek[0,1]] == -9999:
                                orowcol[a,0] = Strchek[0,0]-1
                                orowcol[a,1] = Strchek[0,1]-0
                                orowcol[a,2] = Str[Strchek[0,0]-1,Strchek[0,1]-0]
                                a = a+1
                            if Str[Strchek[0,0],Strchek[0,1]+1] != -9999 and lake[Strchek[0,0],Strchek[0,1]+1] == -9999:
                                orowcol[a,0] = Strchek[0,0]-0
                                orowcol[a,1] = Strchek[0,1]+1
                                orowcol[a,2] = Str[Strchek[0,0]-0,Strchek[0,1]+1]
                                a = a+1
                            if Str[Strchek[0,0],Strchek[0,1]-1] != -9999 and lake[Strchek[0,0],Strchek[0,1]-1] == -9999:
                                orowcol[a,0] = Strchek[0,0]-0
                                orowcol[a,1] = Strchek[0,1]-1
                                orowcol[a,2] = Str[Strchek[0,0]-0,Strchek[0,1]-1]
                                a = a+1
                            if Str[Strchek[0,0]+1,Strchek[0,1]-1] != -9999 and lake[Strchek[0,0]+1,Strchek[0,1]-1] == -9999:
                                orowcol[a,0] = Strchek[0,0]+1
                                orowcol[a,1] = Strchek[0,1]-1
                                orowcol[a,2] = Str[Strchek[0,0]+1,Strchek[0,1]-1]
                                a = a+1
                            if Str[Strchek[0,0]+1,Strchek[0,1]+1] != -9999 and lake[Strchek[0,0]+1,Strchek[0,1]+1] == -9999:
                                orowcol[a,0] = Strchek[0,0]+1
                                orowcol[a,1] = Strchek[0,1]+1
                                orowcol[a,2] = Str[Strchek[0,0]+1,Strchek[0,1]+1]
                                a = a+1
                            if Str[Strchek[0,0]+1,Strchek[0,1]] != -9999 and lake[Strchek[0,0]+1,Strchek[0,1]] == -9999:
                                orowcol[a,0] = Strchek[0,0]+1
                                orowcol[a,1] = Strchek[0,1]-0
                                orowcol[a,2] = Str[Strchek[0,0]+1,Strchek[0,1]-0]
                                a = a+1
                            if a > 0:
                                for ka in range(0,a):
                                    nostr =orowcol[ka,2]
                                    srowcol = np.argwhere(Str==nostr).astype(int)
                                    snrow = srowcol.shape[0]
                                    iStrchek = np.full((snrow,4),-9999)##### 0 row, 1 col, 2 fac,
                                    iStrchek[:,0] = srowcol[:,0]
                                    iStrchek[:,1] = srowcol[:,1]
                                    iStrchek[:,2] = fac[srowcol[:,0],srowcol[:,1]]
                                    iStrchek = iStrchek[iStrchek[:,2].argsort()]
                                    noout = Checklake(iStrchek[0,0],iStrchek[0,1],nrows,ncols,lid,lake)
                                    Lakinstr = np.full(snrow,-9999)
                                    Lakinstr[:] = lake[srowcol[:,0],srowcol[:,1]]
                                    d = np.argwhere(Lakinstr==lid).astype(int)  #### the connected stream should not within the lake
                                    if len(d) < 1 and noout == 0:
                                        GP_cat[orowcol[ka,0],orowcol[ka,1]] = bsid
                                        bsid = bsid + 1
################################################################################
################## Part 3Get Lake pourpoint id and remove cat pourpoint that contribute to lake
#    writeraster(outFolder+"Pourpoints_2.asc",GP_cat)
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        arcatid = np.unique(cat[lrowcol[:,0],lrowcol[:,1]]) ## Get all catchment that intercept with lake
        lakeacc = np.full((len(lrowcol),3),-9999)
        lakeacc[:,0] = lrowcol[:,0]
        lakeacc[:,1] = lrowcol[:,1]
        lakeacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakeacc = lakeacc[lakeacc[:,2].argsort()]
        maxcatid = cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]] ### Get the catchment id that will keeped
        if lakeacc[len(lakeacc)-1,0] != 0 and lakeacc[len(lakeacc)-1,0] != nrows -1 and lakeacc[len(lakeacc)-1,1] != 0 and lakeacc[len(lakeacc)-1,1] != ncols-1:
            GP_cat[lakeacc[len(lakeacc)-1,0]-1,lakeacc[len(lakeacc)-1,1]-1]=-9999 ### remove the all pourpoints close to lake pourpoints
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]-1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]-1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]+1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]+1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]-1,lakeacc[len(lakeacc)-1,1]]=-9999
        for j in range(0,len(arcatid)):
            if arcatid[j] != maxcatid:
                crowcol = np.argwhere(cat==arcatid[j]).astype(int)
                checkcat = np.full((len(crowcol),4),-9999)
                checkcat[:,0] = crowcol[:,0]
                checkcat[:,1] = crowcol[:,1]
                checkcat[:,2] = GP_cat[crowcol[:,0],crowcol[:,1]]
                checkcat[:,3] = fac[crowcol[:,0],crowcol[:,1]]
                checkcat = checkcat[checkcat[:,3].argsort()]
                dele = checkcat[checkcat[:,2]<blid].astype(int)   #### do not delete the lake and strem ids
                GP_cat[dele[:,0],dele[:,1]]=-9999
                nrow,ncol = Nextcell(hydir,checkcat[len(checkcat)-1,0],checkcat[len(checkcat)-1,1])
                if nrow > 0 or ncol >0:
                    if nrow >= nrows or ncols >= ncols:
                        continue
                    if cat[nrow,ncol] < 0:
                        GP_cat[checkcat[len(checkcat)-1,0],checkcat[len(checkcat)-1,1]] = bcid
                        bcid = bcid + 1
        GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]]= sblid
        sblid = sblid + 1
    return GP_cat
###################################################################3


def Nextcell(N_dir,N_row,N_col):
    if N_dir[N_row,N_col] == 1:
        N_nrow = N_row + 0
        N_ncol = N_col + 1
    elif N_dir[N_row,N_col] == 2:
        N_nrow = N_row + 1
        N_ncol = N_col + 1
    elif N_dir[N_row,N_col] == 4:
        N_nrow = N_row + 1
        N_ncol = N_col + 0
    elif N_dir[N_row,N_col] == 8:
        N_nrow = N_row + 1
        N_ncol = N_col - 1
    elif N_dir[N_row,N_col] == 16:
        N_nrow = N_row + 0
        N_ncol = N_col - 1
    elif N_dir[N_row,N_col] == 32:
        N_nrow = N_row - 1
        N_ncol = N_col - 1
    elif N_dir[N_row,N_col] == 64:
        N_nrow = N_row - 1
        N_ncol = N_col + 0
    elif N_dir[N_row,N_col] == 128:
        N_nrow = N_row - 1
        N_ncol = N_col + 1
    else:
        N_nrow = -9999
        N_ncol = -9999
    return N_nrow,N_ncol

##################################################################3
def Getbasinoutlet(ID,basin,fac,dir,nrows,ncols):
    import numpy as np
    catrowcol = np.argwhere(basin==ID).astype(int)
    catacc = np.full((len(catrowcol),3),-9999)
    catacc[:,0] = catrowcol[:,0]
    catacc[:,1] = catrowcol[:,1]
    catacc[:,2] = fac[catrowcol[:,0],catrowcol[:,1]]
    catacc = catacc[catacc[:,2].argsort()]
    ### check if it is a real basin outlet 
    crow = catacc[len(catrowcol)-1,0]
    ccol = catacc[len(catrowcol)-1,1]
          
    nrow,ncol =  Nextcell(dir,crow,ccol)
    
    if nrow < 0 or ncol < 0:
        return crow, ccol
    elif nrow >= nrows or ncol >= ncols:
        return crow, ccol
    elif basin[nrow,ncol] < 0:
        return crow, ccol
    elif basin[nrow,ncol] != ID:   #  all above means the outlet is the real loutlet 
        return crow, ccol
    else:
        crow = nrow 
        ccol = ncol 
        ifound = 0
        for i in range(0,1000): #### find next 1000 grids, to find the basin outlet 
            nrow,ncol =  Nextcell(dir,crow,ccol)
            if nrow < 0 or ncol < 0:
                ifound = 1
                break
            elif nrow >= nrows or ncol >= ncols:
                ifound = 1
                break
            elif basin[nrow,ncol] < 0:
                ifound = 1
                break
            elif basin[nrow,ncol] != ID:
                ifound =  1 #     all above means the outlet is the real loutlet 
                break
            else:
                crow = nrow
                ccol = ncol
                continue
        if ifound == 0: 
            print(" true basin outlet not found for ID...."+ str(ID))
        return crow,ccol        


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
###########


    
class LRRT:
    def __init__(self, dem_in = '#',dir_in = '#',hyshdply = '#',WidDep = '#',Lakefile = '#'
                                     ,Landuse = '#',Landuseinfo = '#',obspoint = '#',OutHyID = -1 ,OutHyID2 = -1, 
                                     OutputFolder = '#',ProjectNM = '#'):
        self.Path_dem_in = dem_in
        self.Path_dir_in = dir_in
        self.Path_hyshdply_in = hyshdply
        self.Path_WiDep_in = WidDep
        self.Path_Lakefile_in = Lakefile
        self.Path_Landuse_in = Landuse
        self.Path_Landuseinfo_in = Landuseinfo
        self.Path_obspoint_in = obspoint
        self.Path_OutputFolder = OutputFolder
        

        
        self.OutHyID = OutHyID
        self.OutHyID2 = OutHyID2
        
        self.ProjectNM = ProjectNM
        
        self.OutputFolder = os.path.join(self.Path_OutputFolder,self.ProjectNM)
        
        if not os.path.exists(self.OutputFolder):
	           os.makedirs(self.OutputFolder)         
        
        self.Raveinputsfolder = self.OutputFolder + '/'+'RavenInput/'
                
        self.qgisPP = os.environ['QGISPrefixPath']
        self.RoutingToolPath = os.environ['RoutingToolFolder']
                
        self.grassdb =os.path.join(tempfile.gettempdir(), 'grassdata_toolbox',self.ProjectNM)
        if not os.path.exists(self.grassdb):
	           os.makedirs(self.grassdb) 
        os.environ['GISDBASE'] = self.grassdb 
		
        self.grass_location_geo = 'Geographic'
        self.grass_location_pro = 'Projected'
        
        self.tempfolder = os.path.join(tempfile.gettempdir(), 'grassdata_toolbox_temp',self.ProjectNM)
        
        if not os.path.exists(self.tempfolder):
	           os.makedirs(self.tempfolder)  

        self.Path_Maskply = '#'
        self.Path_dem = os.path.join(self.tempfolder,'dem.tif')
        self.Path_demproj = os.path.join(self.tempfolder,'dem_proj.tif')
        self.Path_allLakeply = os.path.join(self.tempfolder,'Hylake.shp')
        self.Path_WidDepLine = os.path.join(self.tempfolder,'WidDep.shp')
        self.Path_ObsPoint = os.path.join(self.tempfolder,'obspoint.shp')
        self.Path_Landuseinfo = os.path.join(self.tempfolder,'landuseinfo.csv')
        self.Path_allLakeRas = os.path.join(self.tempfolder,'hylakegdal.tif')
        self.Path_finalcatinfo = os.path.join(self.tempfolder,'catinfo.csv')
        self.Path_alllakeinfoinfo = os.path.join(self.tempfolder,'hylake.csv')
		
        self.cellSize = -9.9999
        self.SpRef_in = '#'
        self.ncols = -9999
        self.nrows = -9999
########################################################################################
### Remove tempfolders

#########################################################################################
##### function to enerate mask based on the most donwstream polygon id
##### Output: self.Path_Maskply    
    def Generatmaskregion(self):
        #### g
        ### Set up QGIS enviroment 
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        
        hyinfocsv = self.Path_hyshdply_in[:-3] + "dbf"
        tempinfo = Dbf5(hyinfocsv)
        hyshdinfo = tempinfo.to_dataframe().values 
           
        HydroBasins1 = Defcat(hyshdinfo,self.OutHyID) ### return fid of polygons that needs to be select 
        if self.OutHyID > 0:
            HydroBasins2 = Defcat(hyshdinfo,self.OutHyID2)            
    ###  exculde the Ids in HydroBasins2 from HydroBasins1
            for i in range(len(HydroBasins2)):
                rows =np.argwhere(HydroBasins1 == HydroBasins2[i])
                HydroBasins1 = np.delete(HydroBasins1, rows)
            HydroBasins = HydroBasins1            
        else:
            HydroBasins = HydroBasins1
        
    ### Load HydroSHED Layers 
        hyshedl12 = QgsVectorLayer(self.Path_hyshdply_in, "")
        
    ### Build qgis selection expression
        where_clause = '"HYBAS_ID" IN'+ " ("
        for i in range(0,len(HydroBasins)):
            if i == 0:
                where_clause = where_clause + str(HydroBasins[i])
            else:
                where_clause = where_clause + "," + str(HydroBasins[i])
        where_clause = where_clause + ")"
    
        req = QgsFeatureRequest().setFlags( QgsFeatureRequest.NoGeometry )
        req.setFilterExpression(where_clause)
        it = hyshedl12.getFeatures( req )
        
        ### obtain all feature id of selected polygons
        selectedFeatureID = []
    
        for feature in it:
    #        print(feature.id())
            selectedFeatureID.append(feature.id())
            
        hyshedl12.select(selectedFeatureID)   ### select with polygon id
        
        # Save selected polygons to output 
        _writer = QgsVectorFileWriter.writeAsVectorFormat(hyshedl12, os.path.join(self.tempfolder, 'HyMask.shp'), "UTF-8", hyshedl12.crs(), "ESRI Shapefile", onlySelected=True)
        processing.run('gdal:dissolve', {'INPUT':os.path.join(self.tempfolder, 'HyMask.shp'),'FIELD':'MAIN_BAS','OUTPUT':os.path.join(self.tempfolder, 'HyMask2.shp')})
        self.Path_Maskply = os.path.join(self.tempfolder, 'HyMask2.shp')
        del tempinfo
        del hyshedl12
        Qgs.exitQgis()
        
        return 

##################################################################################################  
#### functions to preprocess data, Output:
##  self.Path_dem,self.cellSize,self.SpRef_in
##  
    def Generateinputdata(self):

        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
    
        r_dem_layer = QgsRasterLayer(self.Path_dem_in, "") ### load DEM raster as a  QGIS raster object to obtain attribute        
        self.cellSize = float(r_dem_layer.rasterUnitsPerPixelX())  ### Get Raster cell size
        self.SpRef_in = r_dem_layer.crs().authid()   ### get Raster spatialReference id
        print("Working with a  sptail reference  :   " , r_dem_layer.crs().description(), "      ", self.SpRef_in )
        print("The cell cize is   ",self.cellSize)
    
        if  self.Path_hyshdply_in != '#':   #### input is using hydroshed DEM and hydroshed polygons 
            params = {'INPUT': self.Path_dem_in,'MASK': self.Path_Maskply,'NODATA': -9999,'ALPHA_BAND': False,'CROP_TO_CUTLINE': True,
                                                                    'KEEP_RESOLUTION': True,
                                                                    'OPTIONS': 'COMPRESS=LZW',
                                                                    'DATA_TYPE': 0,  # Byte
                                                                    'OUTPUT': self.Path_dem}
            dem = processing.run('gdal:cliprasterbymasklayer',params)  #### extract dem
    
        else:
            params = {'INPUT': self.Path_dem_in, 'format': 'GTiff', 'OUTPUT': self.Path_dem}
            processing.run('gdal:translate',params)
            
#### clip vector fiels  
        processing.run("native:clip", {'INPUT':self.Path_Lakefile_in,'OVERLAY':self.Path_Maskply,'OUTPUT':self.Path_allLakeply})
        processing.run("native:clip", {'INPUT':self.Path_WiDep_in,'OVERLAY':self.Path_Maskply,'OUTPUT':self.Path_WidDepLine})
        processing.run("native:clip", {'INPUT':self.Path_obspoint_in,'OVERLAY':self.Path_Maskply,'OUTPUT':self.Path_ObsPoint})
    
###### set up GRASS environment for translate vector to rasters and clip rasters
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session
    
       
        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
    
        grass.run_command("r.import", input = self.Path_dem, output = 'dem', overwrite = True)
        grass.run_command('g.region', raster='dem')
        grass.run_command('r.mask'  , raster='dem', maskcats = '*',overwrite = True)

    
        if self.Path_dir_in != '#':
            grass.run_command("r.external", input = self.Path_dir_in, output = 'dir_in', overwrite = True)
            grass.run_command("r.clip", input = 'dir_in', output = 'dir_Arcgis', overwrite = True, flags = 'r')
            grass.run_command('r.reclass', input='dir_Arcgis',output = 'dir_Grass',rules =os.path.join(self.RoutingToolPath,'Arcgis2GrassDIR.txt'), overwrite = True)
        else:
            grass.run_command('r.watershed',elevation = 'dem', drainage = 'dir_Grass', overwrite = True)
            grass.run_command('r.reclass', input='dir_Grass',output = 'dir_Arcgis',rules = os.path.join(self.RoutingToolPath,'Grass2ArcgisDIR.txt'), overwrite = True)
    
        copyfile(self.Path_Landuseinfo_in, self.Path_Landuseinfo) 
        grass.run_command("r.external", input = self.Path_Landuse_in, output = 'landuse_in', overwrite = True)
        grass.run_command("r.clip", input = 'landuse_in', output = 'landuse', overwrite = True)        
    
        grass.run_command("v.import", input = self.Path_WidDepLine, output = 'WidDep', overwrite = True)
        grass.run_command("v.import", input = self.Path_ObsPoint, output = 'obspoint', overwrite = True)
        grass.run_command("v.import", input = self.Path_allLakeply, output = 'Hylake', overwrite = True)
    
    
    
        os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(self.cellSize) + "  " +str(self.cellSize) +"   " + "\"" +  self.Path_allLakeply +"\""+ "    "+ "\""+self.Path_allLakeRas+"\"")
    
        grass.run_command("r.in.gdal", input = self.Path_allLakeRas, output = 'alllakeraster_in', overwrite = True)
#    grass.run_command("r.null", input = 'alllakeraster_in', setnull = -9999)

        grass.run_command('v.to.rast',input = 'WidDep',output = 'width',use = 'attr',attribute_column = 'WIDTH',overwrite = True)
        grass.run_command('v.to.rast',input = 'WidDep',output = 'depth',use = 'attr',attribute_column = 'DEPTH',overwrite = True)
        grass.run_command('v.to.rast',input = 'WidDep',output = 'qmean',use = 'attr',attribute_column = 'Q_Mean2',overwrite = True)
        grass.run_command('v.to.rast',input = 'obspoint',output = 'obs',use = 'attr',attribute_column = 'Obs_ID',overwrite = True)
#    grass.run_command('v.to.rast',input = 'Hylake',output = 'alllake',use = 'attr',attribute_column = 'Hylak_id',overwrite = True)

        grass.run_command('r.mapcalc',expression = 'alllake = int(alllakeraster_in)',overwrite = True)
        grass.run_command("r.null", map = 'alllake', setnull = -9999)
    
        PERMANENT.close()
        del r_dem_layer
        Qgs.exitQgis()
#####################################################################################################

####################################################################################################3
                
    def WatershedDiscretizationToolset(self,accthresold):
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        grass.run_command('g.region', raster='dem')


        grass.run_command('r.accumulate', direction='dir_Grass',format = '45degree',accumulation ='acc_grass',
                  stream = 'str_grass_v',threshold = accthresold, overwrite = True)

        grass.run_command('v.to.rast',input = 'str_grass_v',output = 'str_grass_r2',use = 'cat',overwrite = True)
    
        strtemp_array = garray.array(mapname="str_grass_r2")
        acc_array = garray.array(mapname="acc_grass")
        dirarc_array = garray.array(mapname="dir_Arcgis")
    
    #####Correct stream network
        strids = np.unique(strtemp_array)
        strids = strids[strids >= 0] 
        self.ncols = int(strtemp_array.shape[1])
        self.nrows = int(strtemp_array.shape[0])
 
        for i in range(0,len(strids)):
            strid = strids[i]

            trow,tcol = Getbasinoutlet(strid,strtemp_array,acc_array,dirarc_array,self.ncols,self.nrows)
            nrow,ncol = Nextcell(dirarc_array,trow,tcol)### get the downstream catchment id
            nstrid = strtemp_array[nrow,ncol]
        
            rowcol = np.argwhere(strtemp_array==nstrid).astype(int)
            catacc = np.full((len(rowcol),3),-9999)
            catacc[:,0] = rowcol[:,0]
            catacc[:,1] = rowcol[:,1]
            catacc[:,2] = acc_array[rowcol[:,0],rowcol[:,1]]
            catacc = catacc[catacc[:,2].argsort()]
    
            if nrow != catacc[0,0] or ncol != catacc[0,1]:  ### If the stream 1 did not connect to the begining of the downstream stream (2) 
                nnrow,nncol = Nextcell(dirarc_array,nrow,ncol) ### The last cell of 2 will be changed. 
                if nnrow <= 0 or nncol <=0 or nnrow >=self.nrows or nncol >= self.ncols:
                    continue
                nnstrid = strtemp_array[nnrow,nncol]
                strtemp_array[nrow,ncol] = nnstrid
        ##### end modify stream grid and store new grids
        temparray = garray.array()
        temparray[:,:] = 0
        temparray[:,:] = strtemp_array[:,:]
        temparray.write(mapname="str_grass_rf", overwrite=True)
        grass.run_command('r.null', map='str_grass_rf',setnull=0)
        grass.run_command('r.mapcalc',expression = 'str_grass_rfn = int(str_grass_rf)',overwrite = True)
        grass.run_command('r.thin',input = 'str_grass_rfn', output = 'str_grass_r',overwrite = True)
        grass.run_command('r.stream.basins',direction = 'dir_Grass', stream = 'str_grass_r', basins = 'cat1',overwrite = True)
        grass.run_command('r.to.vect',  input = 'str_grass_r',output = 'str', type ='line' ,overwrite = True)
        
################ check connected lakes  and non connected lakes 
        grass.run_command('v.select',ainput = 'Hylake',binput = 'str_grass_v',output = 'lake_str',overwrite = True)
        grass.run_command('v.out.ogr', input = 'lake_str',output = os.path.join(self.tempfolder, "Connect_lake.shp"),format= 'ESRI_Shapefile',overwrite = True)
        os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(self.cellSize) + "  " +str(self.cellSize) +"   " + "\"" +  os.path.join(self.tempfolder, "Connect_lake.shp") +"\""+ "    "+ "\""+os.path.join(self.tempfolder, "cnhylakegdal.tif")+"\"")
        grass.run_command("r.in.gdal", input = os.path.join(self.tempfolder, "cnhylakegdal.tif"), output = 'cnlakeraster_in', overwrite = True)
        grass.run_command('r.mapcalc',expression = 'Connect_Lake = int(cnlakeraster_in)',overwrite = True)    
        grass.run_command('r.mapcalc',expression = 'Nonconnect_Lake = if(isnull(Connect_Lake),alllake,-9)',overwrite = True)
    
        PERMANENT.close()
###########################################################################################3

############################################################################################
    def AutomatedWatershedsandLakesFilterToolset(self,Thre_Lake_Area_Connect = 0,Thre_Lake_Area_nonConnect = -1,MaximumLakegrids = 10000):

        tempinfo = Dbf5(self.Path_allLakeply[:-3] + "dbf")
        allLakinfo = tempinfo.to_dataframe()
        allLakinfo.to_csv(self.Path_alllakeinfoinfo,index = None, header=True)
        
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        grass.run_command('g.region', raster='dem')
    
        VolThreshold =Thre_Lake_Area_Connect ### lake area thresthold for connected lakes 
        NonConLThres = Thre_Lake_Area_nonConnect ### lake area thresthold for non connected lakes 
        nlakegrids = MaximumLakegrids 
    
        temparray = garray.array()
        temparray[:,:] = -9999
        temparray.write(mapname="tempraster", overwrite=True)
##### begin processing 

        blid = 1000000    #### begining of new lake id
        bcid = 1          ## begining of new cat id of hydrosheds
        bsid = 2000000    ## begining of new cat id of in inflow of lakes
        blid2 = 3000000
        boid = 4000000
    
        noncnlake_arr = garray.array(mapname="Nonconnect_Lake")
        conlake_arr = garray.array(mapname="Connect_Lake")
        cat1_arr = garray.array(mapname="cat1")
        str_array = garray.array(mapname="str_grass_r")
        acc_array = garray.array(mapname="acc_grass")
        dir_array = garray.array(mapname="dir_Arcgis")
        obs_array = garray.array(mapname="obs")

###### generate selected lakes 
        hylake1 = selectlake2(conlake_arr,VolThreshold,allLakinfo) ### remove lakes with lake area smaller than the VolThreshold from connected lake raster 
        if NonConLThres >= 0:
            Lake1 = selectlake(hylake1,noncnlake_arr,NonConLThres,allLakinfo) ### remove lakes with lake area smaller than the NonConLThres from non-connected lake raster 
        else:
            Lake1 = hylake1
        temparray[:,:] = Lake1[:,:]
        temparray.write(mapname="SelectedLakes", overwrite=True)
        grass.run_command('r.null', map='SelectedLakes',setnull=-9999)

####    
        Pourpoints = GenerPourpoint(cat1_arr,Lake1,str_array,self.nrows,self.ncols,blid,bsid,bcid,acc_array,dir_array)
        temparray[:,:] = Pourpoints[:,:]
        temparray.write(mapname="Pourpoints_1", overwrite=True)
        grass.run_command('r.null', map='Pourpoints_1',setnull=-9999)
        grass.run_command('r.to.vect', input='Pourpoints_1',output='Pourpoints_1_F',type='point', overwrite = True)
        grass.run_command('r.stream.basins',direction = 'dir_Grass', points = 'Pourpoints_1_F', basins = 'cat2',overwrite = True)
    
        cat2_array =  garray.array(mapname="cat2")
        temcat,outlakeids =CE_mcat4lake(cat2_array,Lake1,acc_array,dir_array,bsid,self.nrows,self.ncols,Pourpoints) 
        temcat2 = CE_Lakeerror(acc_array,dir_array,Lake1,temcat,bsid,blid,boid,self.nrows,self.ncols,cat1_arr)
        temparray[:,:] = temcat2[:,:]
        temparray.write(mapname="cat3", overwrite=True)
        grass.run_command('r.null', map='cat3',setnull=-9999)
    
        nPourpoints = GenerateFinalPourpoints(acc_array,dir_array,Lake1,temcat2,bsid,blid,boid,self.nrows,self.ncols,cat1_arr,obs_array)
        temparray[:,:] = Pourpoints[:,:]
        temparray.write(mapname="Pourpoints_2", overwrite=True)
        grass.run_command('r.null', map='Pourpoints_2',setnull=-9999)
        grass.run_command('r.to.vect', input='Pourpoints_2',output='Pourpoints_2_F',type='point', overwrite = True)    
    
        ### Modify lake flow directions
        ndir = ChangeDIR(dir_array,Lake1,acc_array,self.nrows,self.ncols,outlakeids,nlakegrids)
        temparray[:,:] = ndir[:,:]
        temparray.write(mapname="ndir_Arcgis", overwrite=True)
        grass.run_command('r.null', map='ndir_Arcgis',setnull=-9999)    
        grass.run_command('r.reclass', input='ndir_Arcgis',output = 'ndir_Grass',rules =os.path.join(self.RoutingToolPath,'Arcgis2GrassDIR.txt'),overwrite = True)
     
    
        grass.run_command('r.stream.basins',direction = 'ndir_Grass', points = 'Pourpoints_2_F', basins = 'cat4',overwrite = True)
        cat4_array =  garray.array(mapname="cat4")
        rowcols = np.argwhere(cat4_array == 0)
        cat4_array[rowcols[:,0],rowcols[:,1]] = -9999
        finalcat = CE_mcat4lake2(cat4_array,Lake1,acc_array,dir_array,bsid,self.nrows,self.ncols,nPourpoints)
        temparray[:,:] = finalcat[:,:]
        temparray.write(mapname="finalcat", overwrite=True)
        grass.run_command('r.null', map='finalcat',setnull=-9999)    
        grass.run_command('r.to.vect', input='finalcat',output='finalcat_F1',type='area', overwrite = True)    

####   dissolve final catchment polygons    
        grass.run_command('v.db.addcolumn', map= 'finalcat_F1', columns = "Gridcode VARCHAR(40)")
        grass.run_command('v.db.update', map= 'finalcat_F1', column = "Gridcode",qcol = 'value')
#    grass.run_command('v.reclass', input= 'finalcat_F1', column = "Gridcode",output = 'finalcat_F',overwrite = True)
        grass.run_command('v.dissolve', input= 'finalcat_F1', column = "Gridcode",output = 'finalcat_F',overwrite = True)    
#    grass.run_command('v.select',ainput = 'str',binput = 'finalcat_F',output = 'str_finalcat',overwrite = True)
        grass.run_command('v.overlay',ainput = 'str',binput = 'finalcat_F1',operator = 'and',output = 'str_finalcat',overwrite = True)  
    
        grass.run_command('r.to.vect', input='SelectedLakes',output='SelectedLakes_F',type='area', overwrite = True) 
#        grass.run_command('v.out.ogr', input = 'SelectedLakes_F',output = OutputFolder  + 'SelectedLakes.shp',format= 'ESRI_Shapefile',overwrite = True)
#    grass.run_command('v.db.join', map= 'SelectedLakes_F',column = 'value', other_table = 'result',other_column ='SubId', overwrite = True)

        PERMANENT.close()
                
############################################################################3

    def RoutingNetworkTopologyUpdateToolset(self,projection = 'default'):
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
    
    
        if projection != 'default':

            os.system('gdalwarp ' + "\"" + self.Path_dem +"\""+ "    "+ "\""+self.Path_demproj+"\""+ ' -t_srs  ' + "\""+projection+"\"")

            project = Session()
            project.open(gisdb=self.grassdb, location=self.grass_location_pro,create_opts=projection)
            grass.run_command("r.import", input = self.Path_demproj, output = 'dem_proj', overwrite = True)
            grass.run_command('g.region', raster='dem_proj')  
        
            grass.run_command('v.proj', location=self.grass_location_geo,mapset = 'PERMANENT', input = 'str_finalcat',overwrite = True)
            grass.run_command('v.proj', location=self.grass_location_geo,mapset = 'PERMANENT', input = 'finalcat_F',overwrite = True)
        
            grass.run_command('v.db.addcolumn', map= 'finalcat_F', columns = "Area double precision") 
            grass.run_command('v.db.addcolumn', map= 'str_finalcat', columns = "Length double precision")
              
            grass.run_command('v.to.db', map= 'finalcat_F',option = 'area',columns = "Area", units = 'meters') 
            grass.run_command('v.to.db', map= 'str_finalcat',option = 'length', columns = "Length",units = 'meters')
        
            grass.run_command('r.slope.aspect', elevation= 'dem_proj',slope = 'slope',aspect = 'aspect',precision = 'DCELL',overwrite = True)
            grass.run_command('r.mapcalc',expression = 'tanslopedegree = tan(slope) ',overwrite = True) 
            project.close 
        
            PERMANENT = Session()
            PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
            grass.run_command('g.region', raster='dem')
            grass.run_command('v.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'str_finalcat',overwrite = True)
            grass.run_command('v.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'finalcat_F',overwrite = True) 
            grass.run_command('r.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'tanslopedegree',overwrite = True) 
            grass.run_command('r.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'aspect',overwrite = True) 
             
        else:
            PERMANENT = Session()
            PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
            grass.run_command('g.region', raster='dem')
        
            grass.run_command('v.db.addcolumn', map= 'finalcat_F', columns = "Area double precision") 
            grass.run_command('v.db.addcolumn', map= 'str_finalcat', columns = "Length double precision")
              
            grass.run_command('v.to.db', map= 'finalcat_F',option = 'area',columns = "Area", units = 'meters') 
            grass.run_command('v.to.db', map= 'str_finalcat',option = 'length', columns = "Length",units = 'meters')
        
            grass.run_command('r.slope.aspect', elevation= 'dem_proj',slope = 'slope',aspect = 'aspect',precision = 'DCELL',overwrite = True)
            grass.run_command('r.mapcalc',expression = 'tanslopedegree = tan(slope) ',overwrite = True)         

        grass.run_command('v.to.rast',input = 'finalcat_F',output = 'Area',use = 'attr',attribute_column = 'Area',overwrite = True)  
        grass.run_command('v.to.rast',input = 'str_finalcat',output = 'Length',use = 'attr',attribute_column = 'Length',overwrite = True)
    
#        grass.run_command('r.out.gdal', input = 'Length',output = outputFolder  + 'rivlength.tif',format= 'GTiff',overwrite = True)
#        grass.run_command('r.out.gdal', input = 'finalcat',output = outputFolder  + 'finalcat.tif',format= 'GTiff',overwrite = True)
#        grass.run_command('r.out.gdal', input = 'tanslopedegree',output = outputFolder  + 'slope.tif',format= 'GTiff',overwrite = True)
#        grass.run_command('r.out.gdal', input = 'aspect',output = outputFolder  + 'aspect.tif',format= 'GTiff',overwrite = True)    
#        grass.run_command('v.out.ogr', input = 'str_finalcat',output = outputFolder  + 'str_finalcat.shp',format= 'ESRI_Shapefile',overwrite = True)


#####
        tempinfo = Dbf5(self.Path_allLakeply[:-3] + "dbf")
        allLakinfo = tempinfo.to_dataframe()
        landuseinfo = pd.read_csv(self.Path_Landuseinfo_in,sep=",",low_memory=False)
    
        finalcat_arr = garray.array(mapname="finalcat")
        acc_array = garray.array(mapname="acc_grass")
        dir_array = garray.array(mapname="dir_Arcgis")#ndir_Arcgis
        Lake1_arr = garray.array(mapname="SelectedLakes")    
        dem_array = garray.array(mapname="dem")
        rivlen_array = garray.array(mapname="Length")
        area_array = garray.array(mapname="Area")
        width_array = garray.array(mapname="width")
        depth_array = garray.array(mapname="depth")
        obs_array = garray.array(mapname="obs")
        Q_mean_array = garray.array(mapname="qmean")
        slope_array = garray.array(mapname="tanslopedegree")
        landuse_array = garray.array(mapname="landuse")
    
    
        temparray = garray.array()
        temparray[:,:] = -9999
    

        allcatid = np.unique(finalcat_arr)
        allcatid = allcatid[allcatid >= 0]
        catinfo2 = np.full((len(allcatid),18),-9999.00000)
    
        catinfodf = pd.DataFrame(catinfo2, columns = ['SubId', "DowSubId","Rivlen",'RivSlope','BasinSlope',
                            'BkfWidth','BkfDepth','IsLake','HyLakeId','LakeVol','LakeDepth',
                             'LakeArea','Laketype','IsObs','MeanElev','FloodP_n','Q_Mean','Ch_n'])
                                 
        catinfo,rivpath= Generatecatinfo(finalcat_arr,acc_array,dir_array,Lake1_arr,dem_array,
             area_array,catinfodf,allcatid,allLakinfo,width_array,depth_array,rivlen_array,obs_array,self.nrows,self.ncols,
             slope_array,landuse_array,landuseinfo,Q_mean_array)
        catinfo = Writecatinfotodbf(catinfo)
        
        catinfo.to_csv(self.Path_finalcatinfo, index = None, header=True)
    
    
        temparray[:,:] = rivpath[:,:]
        temparray.write(mapname="rivpath", overwrite=True)
        grass.run_command('r.null', map='rivpath',setnull=-9999)
        grass.run_command('db.in.ogr', input=self.Path_alllakeinfoinfo,output = 'alllakeinfo',overwrite = True)
        grass.run_command('v.db.join', map= 'SelectedLakes_F',column = 'value', other_table = 'alllakeinfo',other_column ='Hylak_id', overwrite = True)
        grass.run_command('db.in.ogr', input=self.Path_finalcatinfo,output = 'result',overwrite = True)
        grass.run_command('v.db.join', map= 'finalcat_F',column = 'Gridcode', other_table = 'result',other_column ='SubId', overwrite = True)
        grass.run_command('v.out.ogr', input = 'finalcat_F',output = os.path.join(self.tempfolder,'finalcat_info1.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
        
#        grass.run_command('v.out.ogr', input = 'finalcat_F',output = outputFolder  + 'finalcat_info.shp',format= 'ESRI_Shapefile',overwrite = True)
        PERMANENT.close
        
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        
        processing.run("native:centroids", {'INPUT':os.path.join(self.tempfolder,'finalcat_info1.shp'),'ALL_PARTS':False,'OUTPUT':os.path.join(self.tempfolder,'Centerpoints.shp')})
        processing.run("native:addxyfields", {'INPUT':os.path.join(self.tempfolder,'Centerpoints.shp'),'CRS':QgsCoordinateReferenceSystem(self.SpRef_in),'OUTPUT':os.path.join(self.tempfolder,'ctwithxy.shp')})
        processing.run("native:joinattributestable",{'INPUT':os.path.join(self.tempfolder,'finalcat_info1.shp'),'FIELD':'SubId','INPUT_2':os.path.join(self.tempfolder,'ctwithxy.shp'),'FIELD_2':'SubId',
                          'FIELDS_TO_COPY':['x','y'],'METHOD':0,'DISCARD_NONMATCHING':False,'PREFIX':'centroid_','OUTPUT':os.path.join(self.OutputFolder,'finalcat_info.shp')})
        Qgs.exitQgis()
        
###########################################################################3
    def Output_Clean(self,Out = 'Simple',clean = 'True',):
        
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        grass.run_command('g.region', raster='dem')
                
        if Out == 'Simple':
#            grass.run_command('v.out.ogr', input = 'finalcat_F',output = os.path.join(self.OutputFolder,'finalcat_info.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
            grass.run_command('v.out.ogr', input = 'SelectedLakes_F',output = os.path.join(self.OutputFolder,'SelectedLakes.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
            grass.run_command('v.out.ogr', input = 'str_finalcat',output = os.path.join(self.OutputFolder,'Channel.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
            grass.run_command('v.out.ogr', input = 'Hylake',output = os.path.join(self.OutputFolder,'AllLakes.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
        if Out == 'All':
            grass.run_command('r.out.gdal', input = 'SelectedLakes',output = os.path.join(self.OutputFolder,'SelectedLakes.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture')
        if clean == 'True':
            shutil.rmtree(self.grassdb,ignore_errors=True)
            shutil.rmtree(self.tempfolder,ignore_errors=True)


    def GenerateRavenInput(self,Finalcat_NM = 'finalcat_info',lenThres = 1,iscalmanningn = -1):
        
        if not os.path.exists(self.Raveinputsfolder):
            os.makedirs(self.Raveinputsfolder)

        finalcatchpath = os.path.join(self.OutputFolder,Finalcat_NM)
        
        tempinfo = Dbf5( finalcatchpath + ".dbf")#np.genfromtxt(hyinfocsv,delimiter=',')
        ncatinfo = tempinfo.to_dataframe().astype(float)
        ncatinfo2 = ncatinfo.drop_duplicates('SubId', keep='first')
        Writervhchanl(ncatinfo2,self.Raveinputsfolder,lenThres,iscalmanningn)
        writelake(ncatinfo2,self.Raveinputsfolder)
        print(self.Raveinputsfolder)
        
        