import numpy as np
import pandas as pd
import copy
import os
import sqlite3
import urllib

def DownloadStreamflowdata_CA(Station_NM,CA_HYDAT,StartYear,EndYear):
    obtaindata   = 1
    con          = sqlite3.connect(CA_HYDAT)
    ### obtain station info
    sqlstat      = "SELECT STATION_NUMBER, DRAINAGE_AREA_GROSS, DRAINAGE_AREA_EFFECT from STATIONS WHERE STATION_NUMBER=?"
    Station_info = pd.read_sql_query(sqlstat, con,params=[Station_NM])
    if len(Station_info) == 0:
        flowdata = -1
        obs_DA   = -9999
        obtaindata = -1 
        return flowdata,obs_DA,obtaindata
    
    DAS          = np.array([-1.2345,Station_info['DRAINAGE_AREA_GROSS'].values[0],Station_info['DRAINAGE_AREA_EFFECT'].values[0]]) 
    DAS          = DAS[DAS != None]
    if len(DAS) > 0:
        obs_DA       = np.nanmax(DAS)
    else:
        obs_DA       = -1.2345


    ## obtain streamflow data
    sqlstat = "select * from DLY_FLOWS WHERE STATION_NUMBER = ?"
    Readed_Streamflow = pd.read_sql_query(sqlstat, con, params=[Station_NM])
    
    Readed_Streamflow = Readed_Streamflow[Readed_Streamflow['YEAR'] >= StartYear]
    Readed_Streamflow = Readed_Streamflow[Readed_Streamflow['YEAR'] <= EndYear]
    ## Initial dataframe 
    if len(Readed_Streamflow) == 0:
        flowdata = -1
        obs_DA   = -9999
        obtaindata = -1 
        return flowdata,obs_DA,obtaindata
        
    year_ini  = Readed_Streamflow['YEAR'].values[0]
    mon_ini   = Readed_Streamflow['MONTH'].values[0]
    year_end  = Readed_Streamflow['YEAR'].values[len(Readed_Streamflow) - 1]
    mon_end   = Readed_Streamflow['MONTH'].values[len(Readed_Streamflow) - 1]
    ndays_end = Readed_Streamflow['NO_DAYS'].values[len(Readed_Streamflow) - 1]
    
    Date_ini  = str(year_ini)+'-'+str(mon_ini)+'-'+'01'
    Date_end  = str(year_end)+'-'+str(mon_end)+'-'+str(ndays_end)
    Date      = pd.date_range(start=Date_ini, end=Date_end, freq='D')
    flowdata  = pd.DataFrame(np.full((len(Date),2),-1.2345),columns = ['Flow','QC'],index = Date)
    
    ### loop read streamflow data 
    for index,row in Readed_Streamflow.iterrows():
        NDays = row['NO_DAYS']
        for iday in range(1,NDays+1):
            cdate = pd.to_datetime({'year': [row['YEAR']],'month': [row['MONTH']],'day': [iday]}).values
#            cdates = pd.to_datetime(str(row['YEAR'])+'-'+str(row['MONTH'])+'-'+str(iday))
            if row['FLOW'+str(iday)] != np.nan and row['FLOW'+str(iday)] != None and float(row['FLOW'+str(iday)]) > 0:
                flowdata.loc[cdate,'Flow'] = row['FLOW'+str(iday)]
                flowdata.loc[cdate,'QC']   = row['FLOW_SYMBOL'+str(iday)]
    return flowdata,obs_DA,obtaindata
            
    
def DownloadStreamflowdata_US(Station_NM,StartYear,EndYear):
    obtaindata   = 1
    #### Obtain station info 
    urlstlist    = 'https://waterdata.usgs.gov/nwis/inventory/?format=rdb&site_no='+str(int(Station_NM)).zfill(8)
    Reslist      =  urllib.request.urlopen(urlstlist)
    stlistdata   = Reslist.read()   
    stlistdata   = stlistdata.splitlines()
    station_info_name  = stlistdata[len(stlistdata) - 3].split()
    station_info_value = stlistdata[len(stlistdata) - 1].split()
     
    if (station_info_name[len(station_info_name) - 1].decode('utf-8') != 'contrib_drain_area_va'):
        obs_DA = -1.2345 
    else: 
        obs_DA = float(station_info_value[len(station_info_value) - 1].decode('utf-8'))*2.58999   #square miles to square km
    
    ## try to obtain data with in this period 
    Date_ini  = str(StartYear)+'-'+'01'+'-'+'01'
    Date_end  = str(EndYear)  +'-'+'12'+'-'+'31'
    urlstlist = 'https://waterdata.usgs.gov/nwis/dv?cb_00060=on&format=rdb&site_no='+str(int(Station_NM)).zfill(8)+'&referred_module=sw&begin_date='+Date_ini+'&end_date='+Date_end
#    print(urlstlist)
    Reslist =  urllib.request.urlopen(urlstlist)
    stlistdata = Reslist.read()   
    stlistdata = stlistdata.splitlines()
    
    ##obtain start of the data rows
    datarow    = -1
    for i in range(0,len(stlistdata)):
        istlistdata = stlistdata[i].split()
        if istlistdata[0] == '#' or len(istlistdata) != 5:
            continue   
        if istlistdata[1].decode('utf-8') == str(int(Station_NM)).zfill(8):
            datarow = i
            break
    Date_ini  = stlistdata[datarow].split()[2].decode('utf-8')
    Date_end  = stlistdata[len(stlistdata)-1].split()[2].decode('utf-8')
    Date      = pd.date_range(start=Date_ini, end=Date_end, freq='D')
    flowdata  = pd.DataFrame(np.full((len(Date),2),-1.2345),columns = ['Flow','QC'],index = Date)
    for i in range(datarow,len(stlistdata)):
        istlistdata = stlistdata[i].split()
        if len(istlistdata) < 5 or  istlistdata[3].decode('utf-8') == 'Ice':
            continue
        else:
            date  = istlistdata[2].decode('utf-8')
            cdate = pd.to_datetime({'year': [date[0:4]],'month': [date[5:7]],'day': [date[8:10]]}).values
            flowdata.loc[cdate,'Flow'] = float(istlistdata[3].decode('utf-8'))*0.0283168  # cubic feet per second to cubic meters per second
            flowdata.loc[cdate,'QC']   = istlistdata[4].decode('utf-8')     
    return flowdata,obs_DA,obtaindata
    

def Writeobsrvtfile(flowdata,obsnm,outObsfileFolder):
    toobsrvtfile = os.path.join(outObsfileFolder,obsnm['Obs_NM']+'_'+str(obsnm['SubId'])+'.rvt')
    f2 = open(toobsrvtfile, "w")
    f2.write(":ObservationData HYDROGRAPH "+str(obsnm['SubId'])+"   m3/s" + " \n")
    f2.write(flowdata.index[0].strftime('%Y-%m-%d') + "  " + '00:00:00  ' + '1     ' + str(len(flowdata))+ '\n')
    for id in range(0,len(flowdata)):
        f2.write('         ' + str(flowdata['Flow'].values[id])+ '\n')
    f2.write(':EndObservationData'+ '\n')
    f2.close()
        
def Modify_template_rvt(outFolderraven,outObsfileFolder,obsnm):
    toobsrvtfile = os.path.join(outFolderraven,'test.rvt')
    obsflodername = './' + os.path.split(outObsfileFolder)[1]+'/' 
    f2 = open(toobsrvtfile, "a")
    f2.write(":RedirectToFile    "+obsflodername+ obsnm['Obs_NM']+'_'+str(obsnm['SubId'])+'.rvt'+" \n")
    f2.close()
#    asdfadsfadsfadsf
        
def WriteObsfiles(catinfo,outFolderraven,outObsfileFolder,startyear,endyear,CA_HYDAT='#',Template_Folder='#'):
    obsnms    = catinfo[['Obs_NM','SRC_obs','SubId','DA']]
    obsnms    = obsnms.drop_duplicates('Obs_NM', keep='first')
    obsnms.loc[:,'DA']  = obsnms['DA'].values/1000/1000  # m2 to km2    
    index     = obsnms.index
    Date      = pd.date_range(start=str(startyear)+'-'+'01'+'-'+'01', end=str(endyear)+'-'+'12'+'-'+'31', freq='D')
    obsnms['Obs_DA_data']  = -1.2345
    obsnms['Missing_V']    = -1.2345
    
    for idx in index:
        obsnm    = obsnms.loc[idx,:]
        iobs_nm  = obsnms.loc[idx,'Obs_NM']
        iobs_src = obsnms.loc[idx,'SRC_obs']
        flowdata = pd.DataFrame(np.full((len(Date),2),-1.2345),columns = ['Flow','QC'],index = Date)
        if iobs_src[0] == '-':
            continue 
        elif iobs_src  == 'US':
            flowdata_read, DA_obs_data,Finddata    = DownloadStreamflowdata_US(Station_NM = iobs_nm,StartYear = startyear,EndYear = endyear)
        elif iobs_src  == 'CA':
            if CA_HYDAT == '#':
                Finddata = -1
                continue
            else:
                
                flowdata_read, DA_obs_data,Finddata = DownloadStreamflowdata_CA(Station_NM = iobs_nm,CA_HYDAT = CA_HYDAT,StartYear = startyear,EndYear = endyear)
        else:
            print("Country not included yet ")
            continue
        
        ####check if data are founded, and assign it to the output dataframes 
        if Finddata < 0:
            print("not find data")
        else:    
            flowdata.loc[flowdata.index.isin(flowdata_read.index), ['Flow', 'QC']] = flowdata_read[['Flow', 'QC']]
            obsnms.loc[idx,'Obs_DA_data']                                          = DA_obs_data
            obsnms.loc[idx,'Missing_V']                                            = len(flowdata[flowdata['Flow'] == -1.2345])
#        flowdata.loc[flowdata_read.index,'Flow'] = flowdata_read['Flow'].values
#        flowdata.loc[flowdata_read.index,'QC']   = flowdata_read['QC'].values
#        print(flowdata)
        Writeobsrvtfile(flowdata,obsnm,outObsfileFolder)
        if Template_Folder != '#':
            Modify_template_rvt(outFolderraven,outObsfileFolder,obsnm)
        ### Modify rvt file
        obsnms.to_csv(os.path.join(outObsfileFolder,'obsinfo.csv'),sep=',',index=False)
     


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



def writelake(catinfo,outFolderraven,HRU_ID_NM,HRU_Area_NM,Sub_ID_NM):
    f2 = open(os.path.join(outFolderraven,"TestLake.rvh"),"w")
    tab = '       '
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HyLakeId'] > 0 and catinfo.iloc[i]['HRU_Type'] > 0:
            lakeid = int(catinfo.iloc[i]['HyLakeId'])
            catid = catinfo.iloc[i][Sub_ID_NM]
            A     = catinfo.iloc[i][HRU_Area_NM]
            h0 = catinfo.iloc[i]['LakeDepth']
            WeirCoe = 0.6
            hruid = int(catinfo.iloc[i][HRU_ID_NM])
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

def Writervhchanl(ocatinfo,outFolder,lenThres,iscalmanningn,HRU_ID_NM,HRU_Area_NM,Sub_ID_NM):
    catinfo_hru = copy.copy(ocatinfo)
    catinfo     = copy.copy(ocatinfo)
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
    
    catinfo_sub = catinfo.drop_duplicates(Sub_ID_NM, keep='first')### remove duplicated subids, beacuse of including hrus in the dataframe 
    print('Total number of Subbasins are     ' +  str(int((len(catinfo_sub))))+'   '+Sub_ID_NM)
    for i in range(0,len(catinfo_sub)):
        ### Get catchment width and dpeth
        catid = int(catinfo_sub[Sub_ID_NM].values[i])
        temp = catinfo_sub['RivLength'].values[i]
        
        if (float(temp) > lenThres):
            catlen = float(temp)/1000 #### in km
            strRlen = str(catlen)
        else:
            catlen = -9999
            strRlen = 'ZERO-'
        if catinfo_sub['IsLake'].values[i] >= 0 :
            strRlen = 'ZERO-'
        #####################################################3
        Strcat = str(catid)
        if catid == catinfo_sub['DowSubId'].values[i]:
            StrDid = str(-1)
        else:
            StrDid = str(int(catinfo_sub['DowSubId'].values[i]))
            
        pronam = 'Chn_'+ Strcat

        chslope = max(catinfo_sub['RivSlope'].values[i],0.00001)
        
        if chslope < 0:
            chslope = 0.0001234
        
        if catinfo_sub['Ch_n'].values[i] > 0:
            nchn = catinfo_sub['Ch_n'].values[i]
        else:
            nchn = 0.001234
            
        if catinfo_sub['FloodP_n'].values[i] > 0:
            floodn = catinfo_sub['FloodP_n'].values[i]
        else:
            floodn = 0.035
            
        writechanel(pronam,max(catinfo_sub['BkfWidth'].values[i],1),max(catinfo_sub['BkfDepth'].values[i],1),
        chslope,ochn,catinfo_sub['MeanElev'].values[i],floodn,nchn,iscalmanningn)
        
        if catinfo_sub['IsObs'].values[i] > 0 :
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
    
    for i in range(0,len(catinfo_hru.index)):
        
        hruid = int(catinfo_hru[HRU_ID_NM].values[i])
        catslope = catinfo_hru['BasSlope'].values[i]
        cataspect= catinfo_hru['BasAspect'].values[i]
        
        catarea2 = catinfo_hru[HRU_Area_NM].values[i]
        
        StrGid =  str(hruid) #str( catinfo_hru.iloc[i][HRU_Area_NM])+tab
        
        catid = str(int(catinfo_hru[Sub_ID_NM].values[i]))+tab
        
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(catinfo_hru['MeanElev'].values[i])+tab
        lat = str(catinfo_hru['centroid_y'].values[i])+tab
        lon = str(catinfo_hru['centroid_x'].values[i])+tab
        if catinfo_hru['IsLake'].values[i] > 0:
            LAND_USE_CLASS = 'Lake_HRU'+tab
            VEG_CLASS = 'Lake_HRU'+tab
            SOIL_PROFILE ='SOILPROF_Lake'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab        
        else:    
            LAND_USE_CLASS = 'FOREST'+tab
            VEG_CLASS = 'FOREST'+tab
            SOIL_PROFILE ='SOILPROF'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab
            
        SLOPE = str(catslope)+tab
        ASPECT = str(cataspect)+tab
        orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
        
            
        
    orvh.write(":EndHRUs"+"\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return 

