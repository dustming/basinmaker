import numpy as np
import pandas as pd
import copy
import os
import sqlite3
import urllib

def DownloadStreamflowdata_CA(Station_NM,CA_HYDAT,StartYear,EndYear):
    """
    Function that used to obtain streamflow data of certain gauge from HYDAT database

    Inputs: 
    
        Station_NM     (string):     The name of the guage, "05PB019"
        CA_HYDAT       (string):     Path and filename of previously downloaded 
                                     external database containing streamflow observations, 
                                     e.g. HYDAT for Canada ("Hydat.sqlite3").
                                    
        Startyear      (integer):    Start year of simulation. Used to 
                                     read streamflow observations from external databases.
        EndYear        (integer):    End year of simulation. Used to 
                                     read streamflow observations from external databases.  
                                          
    Outputs: 
            
        None 

    Return:
    
        flowdata:      (Dataframe):  obtained streamflow observation dataframe between 
                                     Startyear and EndYear.
        obtaindata:    (integer):    1 indicate sucessfully obtain data, -1 indicate no data are founded
                                     for this gauge
        obs_DA:        (float):      The drainage area of this gauge readed from HYDAT database
           
    """
    
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
    """
    Function that used to obtain streamflow data of certain gauge from USGS website

    Inputs: 
    
        Station_NM     (string):     The name of the guage, "050012032"
                                    
        Startyear      (integer):    Start year of simulation. Used to 
                                     read streamflow observations from external databases.
        EndYear        (integer):    End year of simulation. Used to 
                                     read streamflow observations from external databases.  
                                          
    Outputs: 
            
        None 

    Return:
    
        flowdata:      (Dataframe):  obtained streamflow observation dataframe between 
                                     Startyear and EndYear.
        obtaindata:    (integer):    1 indicate sucessfully obtain data, -1 indicate no data are founded
                                     for this gauge
        obs_DA:        (float):      The drainage area of this gauge readed from HYDAT database
           
    """
    
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
    
    """
    Function that used to write xxx.rvt file for each gauge

    Inputs: 
    
        obsnm             (Dataframe):  Dataframe of observation gauge inculding all subbasin/HRUs
                                        information for this gauge as the station name of this gauge
        flowdata:         (Dataframe):  Obtained streamflow observation dataframe between 
                                        Startyear and EndYear.
                                    
        outObsfileFolder  (String):     Path and name of the output folder to save obervation rvt file
                                        of each gauge 
                                         
                                          
    Outputs: 
            
        xxx.rvt                     -   streamflow observation for each gauge xxx in shpfile 
                                        database will be automatically generagted in 
                                        folder outObsfileFolder. 

    Return:
    
       None
           
    """    
    toobsrvtfile = os.path.join(outObsfileFolder,obsnm['Obs_NM']+'_'+str(obsnm['SubId'])+'.rvt')
    f2 = open(toobsrvtfile, "w")
    f2.write(":ObservationData HYDROGRAPH "+str(obsnm['SubId'])+"   m3/s" + " \n")
    f2.write(flowdata.index[0].strftime('%Y-%m-%d') + "  " + '00:00:00  ' + '1     ' + str(len(flowdata))+ '\n')
    for id in range(0,len(flowdata)):
        f2.write('         ' + str(flowdata['Flow'].values[id])+ '\n')
    f2.write(':EndObservationData'+ '\n')
    f2.close()
        
def Modify_template_rvt(outFolderraven,outObsfileFolder,obsnm):
    
    """
    Function that used to modify raven model rvt file (test.rvt)
    Add  ":RedirectToFile    ./obs/xxx_obs.rvt"  for each gauge in the end of 
    model rvt file (test.rvt)
    
    Inputs: 
    
        obsnm             (Dataframe):  Dataframe of observation gauge inculding all subbasin/HRUs
                                        information for this gauge as the station name of this gauge
        outFolderraven    (String):     Path and name of the output folder of Raven input files
        outObsfileFolder  (String):     Path and name of the output folder to save obervation rvt file
                                        of each gauge 
                                         
                                          
    Outputs: 
            
        test.rvt                      - streamflow observation for each gauge xxx in shpfile 
                                        database will be automatically generagted in 
                                        folder outObsfileFolder. 
    Return:
    
       None
    
    todo:
       
       add name of the raven model base name as a paramter
    """ 
    toobsrvtfile = os.path.join(outFolderraven,'test.rvt')
    obsflodername = './' + os.path.split(outObsfileFolder)[1]+'/' 
    f2 = open(toobsrvtfile, "a")
    f2.write(":RedirectToFile    "+obsflodername+ obsnm['Obs_NM']+'_'+str(obsnm['SubId'])+'.rvt'+" \n")
    f2.close()
#    asdfadsfadsfadsf
        
def WriteObsfiles(catinfo,outFolderraven,outObsfileFolder,startyear,endyear,CA_HYDAT='#',Template_Folder='#'):
    
    """
    Function that used to generate Raven streamflow observation input files. All output will be stored in outFolderraven/obs/

    Inputs: 
        General
           outFolderraven   (string):     Folder path and name that save outputs
           catinfo       (DataFrame):     A dataframe includes all attribute for each HRU
                                          read from polygon shpefile generated by the toolbox   
           
        Parameters needed to define obs.rvt file:
           CA_HYDAT  (string):            (optional) path and filename of previously downloaded 
                                          external database containing streamflow observations, 
                                          e.g. HYDAT for Canada ("Hydat.sqlite3").
           Startyear (integer):           (optional) Start year of simulation. Used to 
                                          read streamflow observations from external databases.
           EndYear   (integer):           (optional) End year of simulation. Used to 
                                          read streamflow observations from external databases.  
           WarmUp    (integer):           (optional) The warmup time (in years) used after 
                                          startyear. Values in output file "obs/xxx.rvt" containing 
                                          observations will be set to NoData value "-1.2345".
               
        Input needed to copy raven template files:
           Template_Folder (string):      Folder name containing raven template files.If it is not 
                                          '#', the rvt file in the model folder will be modified.
                                          folowwling line will be added to the end of rvt file 
                                          for each gauge 
                                          ":RedirectToFile    ./obs/xxx_obs.rvt"    
                                          If Template_Folder is '#', nothing happened. 
                                          
    Outputs: 
            
        model.rvt         - (Optional) modified model rvt files which indicate the forcings 
                            and observation gauges.
        xxx.rvt           - streamflow observation for each gauge xxx in shpfile 
                            database will be automatically generagted in 
                            folder "<DataFolder>/Model/RavenInput/obs/". 
        obsinfo.csv       - information file generated reporting drainage area difference 
                            between observed in shpfile and standard database as well as 
                            number of missing values for each gauge
    Return:
        None
           
    """
        
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
    """
    Function that used to write channel paramter file channel rvt. 
        
    Inputs: 
    
        chname             (String):  the name of the channel for each SubBasins
                                      information for this gauge as the station name of this gauge
        chwd               (String):  channel width
        chdep              (String):  channel depth 
        chslope            (String):  channel slope 
        orchnl             (Object):  python file write object 
        elev               (String):  channel elevation 
        floodn             (String):  channel flood plain manning's coefficient 
        channeln           (String):  main channnel manning's coefficient 
        iscalmanningn      (Bool):    Ture use channeln or False use 0.035 as main channel
                                      manning's coefficient 
                                         
                                          
    Outputs: 
            
        modelchannel.rvp      - contains definition and parameters for channels 
        
    Return:
    
       None
    
    """ 
    
    ### Following SWAT instructions, assume a trapezoidal shape channel, with channel sides has depth and width ratio of 2. zch = 2
    zch = 2
    sidwd = zch * chdep ###river side width
    tab = "          "
    botwd = chwd - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*chwd
        sidwd = 0.5*0.5*chwd
        zch = (chwd - botwd)/2/chdep
    if iscalmanningn == True:
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
    """
    Function that used to generate Raven lake definition TestLake.rvh input files. All output will be stored in outFolderraven

    Inputs: 
        General
           outFolderraven   (string):     Folder path and name that save outputs
           catinfo       (DataFrame):     A dataframe includes all attribute for each HRU
                                          read from polygon shpefile generated by the toolbox   
           
        Parameters needed to define lake definition rvh file:
           HRU_ID_NM     (string):        Column name in Finalcat_NM that defines HRU ID
           HRU_Area_NM   (string):        Column name in Finalcat_NM that defines HRU area 
           Sub_ID_NM     (string):        Column name in Finalcat_NM that defines subbasin ID
           Lake_As_Gauge (integer):       If "1", all lake subbasins will labeled as gauged 
                                          subbasin such that Raven will export lake balance for 
                                          this lake. If "-1", lake subbasin will not be labeled 
                                          as gauge subbasin.
  
    Outputs: 
        TestLake.rvh          - contains definition and parameters of lakes
    Return:
       None
           
    """
        
    f2 = open(os.path.join(outFolderraven,"TestLake.rvh"),"w")
    tab = '       '
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HyLakeId'] > 0 and catinfo.iloc[i]['HRU_Type'] == 1: ## lake hru
            lakeid = int(catinfo.iloc[i]['HyLakeId'])
            catid = catinfo.iloc[i][Sub_ID_NM]
            A     = catinfo.iloc[i][HRU_Area_NM] ### in meters
            h0 = catinfo.iloc[i]['LakeDepth'] ## m
            WeirCoe = 0.6   
            hruid = int(catinfo.iloc[i][HRU_ID_NM])
            Crewd = catinfo.iloc[i]['BkfWidth']  ##3 m 
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

def Writervhchanl(ocatinfo,outFolder,lenThres,iscalmanningn,HRU_ID_NM,HRU_Area_NM,Sub_ID_NM,Lake_As_Gauge = 1):
    
    """
    Function that used to generate raven rvh file and channel rvp file. 
    Output will saved to outFolder 

    Inputs: 
        General
           outFolder     (string):     Folder path and name that save outputs
           ocatinfo      (DataFrame):  A dataframe includes all attribute for each HRU
                                       read from polygon shpefile generated by the toolbox       

        Parameters needed to define rvh file:
           lenThres      (float):      River length threshold; river length smaller than 
                                       this will write as zero in Raven rvh file
           iscalmanningn (integer):    If "1", use manning's coefficient in the shpfile table 
                                       and set to default value (0.035).
                                       If "-1", do not use manning's coefficients.
           HRU_ID_NM     (string):     Column name in Finalcat_NM that defines HRU ID
           HRU_Area_NM   (string)      Column name in Finalcat_NM that defines HRU area 
           Sub_ID_NM     (string):     Column name in Finalcat_NM that defines subbasin ID
           Lake_As_Gauge (integer):    If "1", all lake subbasins will labeled as gauged 
                                       subbasin such that Raven will export lake balance for 
                                       this lake. If "-1", lake subbasin will not be labeled 
                                       as gauge subbasin.     
               
    Outputs: 
        test.rvh              - contains subbasins and HRUs
        modelchannel.rvp      - contains definition and parameters for channels
        
    Return:
        None
    
    todo:
       
       Only two HRU type (lake and land) are defined here. need includes more hru types 
           
    """
            
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
        downcatid= int(catinfo_sub['DowSubId'].values[i])
        temp = catinfo_sub['RivLength'].values[i]
        
        if (float(temp) > lenThres):
            catlen = float(temp)/1000 #### in km
            strRlen = str(catlen)
        else:
            catlen = -9999
            strRlen = 'ZERO-'
        if catinfo_sub['IsLake'].values[i] >= 0:# and catinfo_sub['HRU_Type'].values[i] == 1:
            strRlen = 'ZERO-'
        #####################################################3
        Strcat = str(catid)
        if catid == downcatid:
            StrDid = str(-1)
        elif len(catinfo_sub.loc[catinfo_sub[Sub_ID_NM] == downcatid]) == 0:
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
        
        if catinfo_sub['IsObs'].values[i] > 0:
            Guage = '1'
        elif catinfo_sub['IsLake'].values[i] >= 0 and Lake_As_Gauge == True:# and catinfo_sub['HRU_Type'].values[i] == 1: 
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
        
        catarea2 = catinfo_hru[HRU_Area_NM].values[i]/1000/1000  ### in km2
        
        StrGid =  str(hruid) #str( catinfo_hru.iloc[i][HRU_Area_NM])+tab
        
        catid = str(int(catinfo_hru[Sub_ID_NM].values[i]))+tab
        
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(catinfo_hru['MeanElev'].values[i])+tab
        lat = str(catinfo_hru['centroid_y'].values[i])+tab
        lon = str(catinfo_hru['centroid_x'].values[i])+tab
        if catinfo_hru['IsLake'].values[i] > 0 and catinfo_hru['HRU_Type'].values[i] == 1:
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
    orvh.write(":PopulateHRUGroup Lake_HRUs With LANDUSE EQUALS Lake_HRU" + "\n")
    orvh.write(":PopulateHRUGroup Land_HRUs With LANDUSE EQUALS FOREST" + "\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return 

