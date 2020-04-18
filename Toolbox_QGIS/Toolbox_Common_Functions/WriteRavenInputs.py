import numpy as np
import copy
import os


def WriteObsfiles(catinfo,outFolderraven,startyear,endyear):
    obsnms  = catinfo[['Obs_NM','SRC_obs']]
    obsnms  = obsnms.drop_duplicates('Obs_NM', keep='first')
    print(obsnms)


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



def writelake(catinfo,outFolderraven,nclakeinfo):
    f2 = open(os.path.join(outFolderraven,"TestLake.rvh"),"w")
    tab = '       '
    maxcatid = max(catinfo['SubId'].values)
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HyLakeId'] > 0:
            lakeid = int(catinfo.iloc[i]['HyLakeId'])
            catid = catinfo.iloc[i]['SubId']
            if float(catinfo.iloc[i]['BasArea'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LakeArea']):
                A = float(catinfo.iloc[i]['BasArea'])*0.95
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
    Includ_Nonc_Lakes = nclakeinfo[nclakeinfo['new_HruID'] > 0]
    for inclidx in range(0,len(Includ_Nonc_Lakes)):
            lakeid = int(Includ_Nonc_Lakes['Gridcode'].values[inclidx])
            catid = int(Includ_Nonc_Lakes['new_SubID'].values[inclidx])
            A = float(Includ_Nonc_Lakes['Lake_area'].values[inclidx])*1000*1000
#            A = catinfo.iloc[i]['LakeArea']*1000*1000
            h0 = float(Includ_Nonc_Lakes['Depth_avg'].values[inclidx])
            WeirCoe = 0.6
            hruid = int(Includ_Nonc_Lakes['new_HruID'].values[inclidx])
            Crewd = float(Includ_Nonc_Lakes['BkfWidth'].values[inclidx])
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
                    if float(catinfozone.iloc[i]['BasArea'])/1000.00/1000.00 <= float(catinfozone.iloc[i]['LakeArea']):
                        A = float(catinfozone.iloc[i]['BasArea'])*0.95
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

def Writervhchanl(ocatinfo,outFolder,lenThres,iscalmanningn,nclakeinfo):
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
    nclakeinfo['new_SubID'] = -9999
    nclakeinfo['new_HruID'] = -9999
    nclakeinfo['BkfWidth'] = -9999
    maxcatid = max(catinfo['SubId'].values)
    num_con_lake = len(catinfo.loc[catinfo['IsLake']>0])
    Non_con_lake_sid = 2*maxcatid + num_con_lake + 1000
    nonlake_jjj = 1
    
    for i in range(0,len(catinfo.index)):
        ### Get catchment width and dpeth
        catid = int(catinfo.iloc[i]['SubId'])
        temp = catinfo.iloc[i]['RivLength']
        
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
            chslope = 0.0001234
        
        if catinfo.iloc[i]['Ch_n'] > 0:
            nchn = catinfo.iloc[i]['Ch_n']
        else:
            nchn = 0.001234
            
        if catinfo.iloc[i]['FloodP_n'] > 0:
            floodn = catinfo.iloc[i]['FloodP_n']
        else:
            floodn = 0.035
            
        writechanel(pronam,max(catinfo.iloc[i]['BkfWidth'],1),max(catinfo.iloc[i]['BkfDepth'],1),
        chslope,ochn,catinfo.iloc[i]['MeanElev'],floodn,nchn,iscalmanningn)
        
        if catinfo.iloc[i]['IsObs'] > 0 :
            Guage = '1'
        else:
            Guage = '0'
        orvh.write("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage+"\n")
        
       ###write non connected lake subbasin 
        Non_con_lake_info = nclakeinfo.loc[nclakeinfo['SubId_cat'] == catid]
        for idx_ncl in range(0,len(Non_con_lake_info)):
            nclakeid = Non_con_lake_info['Gridcode'].values[idx_ncl]
            Strcat = str(Non_con_lake_sid + nonlake_jjj) ##### leave maxcatid to 2* maxcatid to connected lake hrus
            StrDid = str(catid)
            strRlen = '0'
            Guage = '0'
            orvh.write("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage+"\n")
            nclakeinfo.loc[nclakeinfo['Gridcode'] == nclakeid,'new_SubID'] = Non_con_lake_sid + nonlake_jjj
            nclakeinfo.loc[nclakeinfo['Gridcode'] == nclakeid,'BkfWidth'] = max(catinfo.iloc[i]['BkfWidth'],1)
            nonlake_jjj = nonlake_jjj + 1
        
     
    orvh.write(":EndSubBasins"+"\n")
    orvh.write("\n")
##########################################
    orvh.write(":HRUs"+"\n")
    orvh.write("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT"+"\n")
    orvh.write("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg"+"\n")
    maxcatid = max(catinfo['SubId'].values)
    maxcatidnclake = max(nclakeinfo['new_SubID'].values) + 1000
    nonlake_kkk = 1    

    
    for i in range(0,len(catinfo.index)):
        hruid = int(catinfo.iloc[i]['SubId'])
        catslope = catinfo.iloc[i]['BasSlope']
        cataspect= catinfo.iloc[i]['BasAspect']
        Non_con_lake_info = nclakeinfo.loc[nclakeinfo['SubId_cat'] == hruid]
        ### calculate Non conneced lake da
        if len(Non_con_lake_info) < 0:
            Non_con_L_DA_area = 0.0
        else:
            Non_con_L_DA_area = np.sum(Non_con_lake_info['Area_m'].values)/1000.00/1000.00
        
        if catinfo.iloc[i]['IsLake'] > 0:
            Con_L_area = float(catinfo.iloc[i]['LakeArea'])
        else:
            Con_L_area = 0.0
        
        landarea = 0.0
        total_water_area = Non_con_L_DA_area + Con_L_area
        
        if float(catinfo.iloc[i]['BasArea'])/1000.00/1000.00 <= total_water_area:
            catarea2 = float(catinfo.iloc[i]['BasArea'])*0.05/1000.00/1000.00
            landarea = float(catinfo.iloc[i]['BasArea'])*0.05/1000.00/1000.00
        else:  ### cat area larger than total water area
            catarea2 = float(catinfo.iloc[i]['BasArea'])/1000.00/1000.00 - total_water_area
            landarea =  float(catinfo.iloc[i]['BasArea'])/1000.00/1000.00 - total_water_area

        Restarea = float(catinfo.iloc[i]['BasArea'])/1000.00/1000.00 - landarea  ## in km2

        StrGid =  str(hruid)+tab
        catid = str(int(catinfo.iloc[i]['SubId']))+tab
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(catinfo.iloc[i]['MeanElev'])+tab
        lat = str(catinfo.iloc[i]['centroid_y'])+tab
        lon = str(catinfo.iloc[i]['centroid_x'])+tab
        LAND_USE_CLASS = 'FOREST'+tab
        VEG_CLASS = 'FOREST'+tab
        SOIL_PROFILE ='SOILPROF'+tab
        AQUIFER_PROFILE ='[NONE]'+tab
        TERRAIN_CLASS ='[NONE]'+tab
        SLOPE = str(catslope)+tab
        ASPECT = str(cataspect)+tab
        orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
        
        if catinfo.iloc[i]['IsLake'] > 0:
            hruid = int(catinfo.iloc[i]['SubId']) + int(maxcatid)
            catslope = catinfo.iloc[i]['BasSlope']
            cataspect= catinfo.iloc[i]['BasAspect']
            
            if Restarea <= float(catinfo.iloc[i]['LakeArea']):
                catarea2 = Restarea*0.95
                Con_L_area =  Restarea*0.95
            else:
                catarea2 = float(catinfo.iloc[i]['LakeArea'])
                Con_L_area =  float(catinfo.iloc[i]['LakeArea']) ### in km2

            StrGid =  str(hruid)+tab
            catid = str(int(catinfo.iloc[i]['SubId']))+tab
            StrGidarea = str(catarea2)+tab
            StrGidelev = str(catinfo.iloc[i]['MeanElev'])+tab
            lat = str(catinfo.iloc[i]['centroid_y'])+tab
            lon = str(catinfo.iloc[i]['centroid_x'])+tab
            LAND_USE_CLASS = 'WATER'+tab
            VEG_CLASS = 'WATER'+tab
            SOIL_PROFILE ='LAKE'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab
            SLOPE = str(catslope)+tab
            ASPECT = str(cataspect)+tab
            orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
            
        Leftarea = float(catinfo.iloc[i]['BasArea'])/1000.00/1000.00 - landarea - Con_L_area
        
        
        for idxnonclake in range(0,len(Non_con_lake_info)):
            total_non_c_lake_cat_a =  Leftarea * (Non_con_lake_info['Area_m'].values[idxnonclake]/1000.00/1000.00)/Non_con_L_DA_area
            nclakeid = Non_con_lake_info['Gridcode'].values[idxnonclake]
            
            if total_non_c_lake_cat_a <= Non_con_lake_info['Lake_area'].values[idxnonclake]:
                catarea2 = total_non_c_lake_cat_a*0.05
            else:  ### cat area larger than total water area
                catarea2 = total_non_c_lake_cat_a - Non_con_lake_info['Lake_area'].values[idxnonclake]
            
            StrGid =  str(int(Non_con_lake_info['new_SubID'].values[idxnonclake]))+tab
            catid = str(int(Non_con_lake_info['new_SubID'].values[idxnonclake]))+tab

            StrGidarea = str(catarea2)+tab
            StrGidelev = str(catinfo.iloc[i]['MeanElev'])+tab
            lat = str(catinfo.iloc[i]['centroid_y'])+tab
            lon = str(catinfo.iloc[i]['centroid_x'])+tab
            LAND_USE_CLASS = 'FOREST'+tab
            VEG_CLASS = 'FOREST'+tab
            SOIL_PROFILE ='SOILPROF'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab
            SLOPE = str(catslope)+tab
            ASPECT = str(cataspect)+tab
            orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")  

            StrGid =  str(int(maxcatidnclake + nonlake_kkk))+tab
#            print(Non_con_lake_info['new_SubID'].values[idxnonclake],maxcatidnclake,nclakeid)
            nclakeinfo.loc[nclakeinfo['Gridcode'] == nclakeid,'new_HruID'] = int(maxcatidnclake + nonlake_kkk)
            nonlake_kkk = nonlake_kkk + 1

            catid = str(int(Non_con_lake_info['new_SubID'].values[idxnonclake]))+tab
            StrGidarea = str(total_non_c_lake_cat_a - catarea2)+tab
            StrGidelev = str(catinfo.iloc[i]['MeanElev'])+tab
            lat = str(catinfo.iloc[i]['centroid_y'])+tab
            lon = str(catinfo.iloc[i]['centroid_x'])+tab
            LAND_USE_CLASS = 'WATER'+tab
            VEG_CLASS = 'WATER'+tab
            SOIL_PROFILE ='LAKE'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab
            SLOPE = str(catslope)+tab
            ASPECT = str(cataspect)+tab
            orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")              
            
        
    orvh.write(":EndHRUs"+"\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return nclakeinfo

