import numpy as np
import pandas as pd
import copy
import os
import sqlite3
import urllib


def WriteStringToFile(Out_String,File_Path,WriteMethod):
    """Write String to a file
    
    Function that used to write Out_String to a file located at the File_Path.

    Parameters 
    ----------
    Out_String            : string
        The string that will be writed to the file located at File_Path
    File_Path             : string
        Path and filename of file that will be modified or created                            
    WriteMethod           : {'a','w'}
        If WriteMethod = "w", a new file will be created at the File_Path
        If WriteMethod = "a", the Out_String will be added to exist file  
                                          
    Notes
    ------    
        The file located at the File_Path will be modified or created  

    Returns
    -------
        None 
   
    Examples
    --------
    >>> from WriteRavenInputs import WriteStringToFile
    >>> Out_String = 'sometest at line 1\n some test at line 2\n some test at line 3\n'
    >>> File_Path  = 'C:/Path_to_the_Flie_with_file_name'
    >>> WriteStringToFile(Out_String = Out_String,File_Path = File_Path,WriteMethod = 'w')
    
    """
    
    if os.path.exists(File_Path): ### if file exist, we can either modify or overwrite it 
        with open(File_Path,WriteMethod) as f:
            f.write(Out_String)
    else: ## create a new file anyway, since file did not exist 
        with open(File_Path,'w') as f:
            f.write(Out_String)    
        
    
        
def DownloadStreamflowdata_CA(Station_NM,CA_HYDAT,StartYear,EndYear):
    """Return streamflow data from HYDAT
    
    Function that used to obtain streamflow data of certain gauge from HYDAT database

    Parameters 
    ----------
    Station_NM            : string
        The name of the gauge, "05PB019"
    CA_HYDAT              : string
        Path and filename of previously downloaded 
        external database containing streamflow observations, 
        e.g. HYDAT for Canada ("Hydat.sqlite3").                                
    Startyear             : integer
        Start year of simulation. Used to 
        read streamflow observations from external databases.
    EndYear               : integer  
        End year of simulation. Used to 
        read streamflow observations from external databases.  
                                          
    Notes
    ------    
        None 

    Returns
    -------
    flowdata              : data-type   
        obtained streamflow observation dataframe between 
        Startyear and EndYear.
    obtaindata            : bool
        True indicate successfully obtain data, False indicate no data are founded
        for this gauge
    obs_DA                : float
        The drainage area of this gauge read from HYDAT database
   
    Examples
    --------
    >>> from WriteRavenInputs import DownloadStreamflowdata_CA
    >>> Station_NM = '05PC019'
    >>> StartYear  = 2010
    >>> EndYear    = 2011
    >>> CA_HYDAT   = HYDAT_Path
    >>> flowdata,obs_DA,obtaindata = DownloadStreamflowdata_CA(Station_NM,CA_HYDAT,StartYear,EndYear)
    
    """
    
    obtaindata   = True
    con          = sqlite3.connect(CA_HYDAT)
    ### obtain station info
    sqlstat      = "SELECT STATION_NUMBER, DRAINAGE_AREA_GROSS, DRAINAGE_AREA_EFFECT from STATIONS WHERE STATION_NUMBER=?"
    Station_info = pd.read_sql_query(sqlstat, con,params=[Station_NM])
    if len(Station_info) == 0:
        flowdata = -1
        obs_DA   = -9999
        obtaindata = False 
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
        obtaindata = False 
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
    """Return streamflow data from USGS website 
    
    Function that used to obtain streamflow data of certain gauge from USGS website

    Parameters 
    ----------
    Station_NM            : string
        The name of the gauge, "05127000"                                
    Startyear             : integer
        Start year of simulation. Used to 
        read streamflow observations from external databases.
    EndYear               : integer  
        End year of simulation. Used to 
        read streamflow observations from external databases.  
                                          
    Notes
    ------    
        None 

    Returns
    -------
    flowdata              : data-type   
        obtained streamflow observation dataframe between 
        Startyear and EndYear.
    obtaindata            : bool
        True indicate successfully obtain data, False indicate no data are founded
        for this gauge
    obs_DA                : float
        The drainage area of this gauge read from HYDAT database
   
    Examples
    --------
    >>> from WriteRavenInputs import DownloadStreamflowdata_US
    >>> Station_NM = '05127000'
    >>> StartYear  = 2010
    >>> EndYear    = 2011
    >>> flowdata,obs_DA,obtaindata = DownloadStreamflowdata_CA(Station_NM,StartYear,EndYear)
    
    """
    
    
    obtaindata   = True
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
    

def Generate_Raven_Obsrvt_String(flowdata,obsnm,outObsfileFolder):   #Writeobsrvtfile(flowdata,obsnm,outObsfileFolder):
    
    """Generate a string in Raven observation rvt input file format 
    
    Function that is used to subbasin id and observation guage name from obsnm
    and reformat the streamflow observation data in flowdata 
    generate a string that follows the raven observation rvt input file format

    Parameters 
    ----------
    flowdata           : data-type 
        Obtained streamflow observation dataframe between 
        Startyear and EndYear. The index of the dataframe should be Date in 
        '%Y-%m-%d' format, and the streamflow observation data in m3/s should 
        in 'Flow' column      
    obsnm              : data-type 
        Dataframe of observation gauge information for this gauge including
        at least following two columns
        'Obs_NM': the name of the stream flow obsrvation gauge 
        'SubId' : the subbasin Id of this stremflow gauge located at.         
    outObsfileFolder   : string     
        Path and name of the output folder to save obervation rvt file
        of each gauge                                          
                                          
    Notes
    ------    
        None 

    Returns
    -------
    obs_rvt_file_path : string
        It is the file path inclding file names of the raven rvt input file
        for this gauge 
    output_string     : string
        It is the string that contains the content of the raven rvt input 
        file of this gauge   

    See Also
    --------
    DownloadStreamflowdata_US : Generate flowdata inputs needed by this function
    DownloadStreamflowdata_CA : Generate flowdata inputs needed by this function 

    Examples
    --------
    >>> from WriteRavenInputs import DownloadStreamflowdata_US,Generate_Raven_Obsrvt_String
    >>> import pandas as pd 
    >>> Station_NM  = '05127000'
    >>> StartYear   = 2010
    >>> EndYear     = 2011
    >>> Subbasin_ID = 1
    >>> flowdata_read, DA_obs_data,Finddata    = DownloadStreamflowdata_US(Station_NM = iobs_nm,StartYear = startyear,EndYear = endyear)
    >>> Date                                   = pd.date_range(start=str(startyear)+'-'+'01'+'-'+'01', end=str(endyear)+'-'+'12'+'-'+'31', freq='D')
    >>> flowdata    = pd.DataFrame(np.full((len(Date),2),-1.2345),columns = ['Flow','QC'],index = Date)
    >>> flowdata.loc[flowdata.index.isin(flowdata_read.index), ['Flow', 'QC']] = flowdata_read[['Flow', 'QC']]
    >>> obsnms = pd.DataFrame(data=[Subbasin_ID,Station_NM],columns=['SubId','Obs_NM'])
    >>> Outputfolderrvt = 'c:/some_folder_to_store_raven_rvt_file'
    >>> obs_rvt_file_path, output_string = Generate_Raven_Obsrvt_String(flowdata = flowdata,obsnm = obsnms,outObsfileFolder = Outputfolderrvt)
    
    """
    
    output_string_list = []    
    obs_rvt_file_path = os.path.join(outObsfileFolder,obsnm['Obs_NM']+'_'+str(obsnm['SubId'])+'.rvt')
    output_string_list.append(":ObservationData HYDROGRAPH "+str(obsnm['SubId'])+"   m3/s")
    output_string_list.append(flowdata.index[0].strftime('%Y-%m-%d') + "  " + '00:00:00  ' + '1     ' + str(len(flowdata)))
    for id in range(0,len(flowdata)):
        output_string_list.append('         ' + str(flowdata['Flow'].values[id]))
    output_string_list.append(':EndObservationData'+ '\n')
    output_string = "\n".join(output_string_list)
    
    return obs_rvt_file_path, output_string
    
def Generate_Raven_Timeseries_rvt_String(outFolderraven,outObsfileFolder,obsnm,Model_Name):  #Modify_template_rvt(outFolderraven,outObsfileFolder,obsnm):
    
    """Generate a string in Raven time series rvt input file format 
    
    Function that used to modify raven model timeseries rvt file (Model_Name.rvt)
    Add  ":RedirectToFile    ./obs/guagename_subbasinid.rvt" 
    for each gauge in the end of model rvt file (Model_Name.rvt)
    
    Parameters 
    ----------
    outFolderraven            : String
        Path and name of the output folder of Raven input files
    outObsfileFolder          : String
        Path and name of the output folder to save obervation rvt file
        of each gauge 
    obsnm                     : data-type 
        Dataframe of observation gauge information for this gauge including
        at least following two columns
        'Obs_NM': the name of the stream flow obsrvation gauge 
        'SubId' : the subbasin Id of this stremflow gauge located at. 
    Model_Name                : string
        The Raven model base name. File name of the raven input will be 
        Model_Name.xxx.                                                 
                                          
    Notes
    ------    
    None 

    See Also
    --------
    DownloadStreamflowdata_US            : Generate flowdata inputs 
                                           needed by this function
    DownloadStreamflowdata_CA            : Generate flowdata inputs 
                                           needed by this function 
    
    Returns
    -------
    output_string     : string
        It is the string that contains the content that will be used to 
        modify the raven time series rvt input file of this gauge   
    Examples
    --------
    >>> from WriteRavenInputs import Generate_Raven_Timeseries_rvt_String
    >>> outFolderraven    = 'c:/path_to_the_raven_input_folder/'
    >>> outObsfileFolder  = 'c:/path_to_the_raven_streamflow_observation gauge_folder/'
    >>> Subbasin_ID = 1
    >>> Station_NM  = '05127000'
    >>> obsnms = pd.DataFrame(data=[Subbasin_ID,Station_NM],columns=['SubId','Obs_NM'])
    >>> Model_Name = 'test'
    >>> output_string = Generate_Raven_Timeseries_rvt_String(outFolderraven,outObsfileFolder,obsnm,Model_Name)
    
    """ 
    
    toobsrvtfile = os.path.join(outFolderraven,Model_Name + '.rvt')
    obsflodername = './' + os.path.split(outObsfileFolder)[1]+'/' 
    output_string = "  \n"+":RedirectToFile    "+obsflodername+ obsnm['Obs_NM']+'_'+str(obsnm['SubId'])+'.rvt' +  "  \n"
    return output_string

        
def Generate_Raven_Obs_rvt_String(catinfo,outFolderraven,outObsfileFolder,startyear,endyear,CA_HYDAT='#',DownLoadObsData=True,Model_Name = 'test'):
    
    """ Generate Raven streamflow observation conent
     
    Function that used to generate content of Raven streamflow observation input file. 

    Parameters 
    ----------
    catinfo                   : data-type
        Dataframe of a routing structure that needes to define a Raven model. 
        Can be directly read from the database of the hru shpfile generated by
        the toolbox. At least include following columns: 
        'Obs_NM'  - the name of the stream flow obsrvation gauge 
        'SubId'   - the subbasin Id of this stremflow gauge located at. 
        'SRC_obs' - the country of the gauge located at 
        'DA'      - the drainage area controlled by the gauge obtained from 
                    the routing structure 
    outFolderraven            : String
        Path and name of the output folder of Raven input files
    outObsfileFolder          : String
        Path and name of the output folder to save obervation rvt file
        of each gauge 
    startyear                 : integer
        Start year of simulation. Used to 
        read streamflow observations from external databases.
    endyear                   : integer  
        End year of simulation. Used to 
        read streamflow observations from external databases.
    CA_HYDAT                  : string,  optional
        path and filename of previously downloaded 
        external database containing streamflow observations, 
        e.g. HYDAT for Canada ("Hydat.sqlite3").  
    DownLoadObsData           : Bool, optional 
        Input that used to indicate if the observation data will be Download
        from usgs website or read from hydat database for streamflow Gauge 
        in US or Canada,respectively. If this parameter is False, 
        while WriteObsrvt is True. The program will write the observation data
        file with "-1.2345" for each observation gauges.
    Model_Name                : string
        The Raven model base name. File name of the raven input will be 
        Model_Name.xxx.                                     
                                          
    Notes
    ------  
    None 
                            
    See Also
    --------
    DownloadStreamflowdata_US            : Generate flowdata inputs 
                                           needed by this function
    DownloadStreamflowdata_CA            : Generate flowdata inputs 
                                           needed by this function 
    Generate_Raven_Obsrvt_String         : Generate a string in Raven 
                                           observation rvt input file format
    Generate_Raven_Timeseries_rvt_String : Generate a string in Raven 
                                           time series rvt input file format 
    
    Returns
    -------
    obs_rvt_file_path_gauge_list                : string
        It is the list of string, each of them contains the content 
        that will be used to define the path of raven observation
        input file for one streamflow gauge 
    obs_rvt_file_string_gauge_list              : string
        It is the list of string, each of them define the content of 
        raven obaervation input file(xxx_obs.rvt) for one gauge
    Model_rvt_file_path                         : string
        It is the string that define the path of
        the raven model time series input file 
    Model_rvt_file_string_modify_gauge_list     : string
        It is the list of string, each of them define the content
        that needs to be added into the Raven model time series 
        file the path of
    obsnm                                       : DataFrame
        Dataframe of observation gauge information for all streamflow gauges  
        
    Examples
    --------
    >>> from WriteRavenInputs import Generate_Raven_Obs_rvt_String
    >>> outFolderraven    = 'c:/path_to_the_raven_input_folder/'
    >>> DataFolder = "C:/Path_to_foldr_of_example_dataset_provided_in_Github_wiki/"
    >>> Model_Folder     = os.path.join(DataFolder,'Model')
    >>> Raveinputsfolder = os.path.join(Model_Folder,'RavenInput')
    >>> Obs_Folder       = os.path.join(Raveinputsfolder,'obs')
    >>> finalcatchpath = os.path.join(DataFolder,'finalcat_hru_info.shp')
    >>> tempinfo = Dbf5(finalcatchpath[:-3] + "dbf")
    >>> ncatinfo = tempinfo.to_dataframe()
    >>> Model_Name = 'test'
    >>> Startyear = 2010
    >>> EndYear = 2017
    >>> CA_HYDAT = 'c/path_to_your_HYDAT_database/'
    >>> WarmUp = 1
    >>> DownLoadObsData = True
    >>> ncatinfo2 = ncatinfo.drop_duplicates('HRU_ID', keep='first')
    >>> obs_rvt_file_path_gauge_list,obs_rvt_file_string_gauge_list,Model_rvt_file_path,Model_rvt_file_string_modify_gauge_list,obsnms = Generate_Raven_Obs_rvt_String(ncatinfo2,Raveinputsfolder,Obs_Folder,
    ...                                                                                                                                                                Startyear + WarmUp,EndYear,CA_HYDAT,
    ...                                                                                                                                                                DownLoadObsData,
    ...                                                                                                                                                                Model_Name)
    >>>

            
    """
        
    obsnms    = catinfo[['Obs_NM','SRC_obs','SubId','DA']]
    obsnms    = obsnms.drop_duplicates('Obs_NM', keep='first')
    obsnms    = obsnms.loc[obsnms['Obs_NM'] != '-9999.0']
    
    obsnms.loc[:,'DA']  = obsnms['DA'].values/1000/1000  # m2 to km2    
    index     = obsnms.index
    Date      = pd.date_range(start=str(startyear)+'-'+'01'+'-'+'01', end=str(endyear)+'-'+'12'+'-'+'31', freq='D')
    obsnms['Obs_DA_data']  = -1.2345
    obsnms['Missing_V']    = -1.2345
    obs_rvt_file_path_gauge_list            = [] 
    obs_rvt_file_string_gauge_list          = []
    Model_rvt_file_path                     = os.path.join(outFolderraven,Model_Name+'.rvt')
    Model_rvt_file_string_modify_gauge_list = []    
    
    for idx in index:
        obsnm    = obsnms.loc[idx,:]
        iobs_nm  = obsnms.loc[idx,'Obs_NM']
        iobs_src = obsnms.loc[idx,'SRC_obs']
        flowdata = pd.DataFrame(np.full((len(Date),2),-1.2345),columns = ['Flow','QC'],index = Date)
        if iobs_src[0] == '-':
            Finddata = False
        elif iobs_src  == 'US':
            if DownLoadObsData == True:
                flowdata_read, DA_obs_data,Finddata    = DownloadStreamflowdata_US(Station_NM = iobs_nm,StartYear = startyear,EndYear = endyear)
            else:
                Finddata = False                
        elif iobs_src  == 'CA':
            if CA_HYDAT != '#' and DownLoadObsData == True:
                flowdata_read, DA_obs_data,Finddata = DownloadStreamflowdata_CA(Station_NM = iobs_nm,CA_HYDAT = CA_HYDAT,StartYear = startyear,EndYear = endyear)
            else:
                Finddata = False              
        else:
            Finddata = False
        
        ####check if data are founded, and assign it to the output dataframes 
        if Finddata == False:
            print("not find data")
        else:    
            flowdata.loc[flowdata.index.isin(flowdata_read.index), ['Flow', 'QC']] = flowdata_read[['Flow', 'QC']]
            obsnms.loc[idx,'Obs_DA_data']                                          = DA_obs_data
            obsnms.loc[idx,'Missing_V']                                            = len(flowdata[flowdata['Flow'] == -1.2345])


        obs_rvt_file_path, output_string = Generate_Raven_Obsrvt_String(flowdata,obsnm,outObsfileFolder)
#        WriteStringToFile(Out_String = output_string,File_Path = obs_rvt_file_path, WriteMethod = "w")

        obs_rvt_file_path_gauge_list.append(obs_rvt_file_path)
        obs_rvt_file_string_gauge_list.append(output_string)
                
        output_string = Generate_Raven_Timeseries_rvt_String(outFolderraven,outObsfileFolder,obsnm,Model_Name)
#        WriteStringToFile(Out_String = output_string,File_Path = os.path.join(outFolderraven,Model_Name+'.rvt'), WriteMethod = "a") 
        Model_rvt_file_string_modify_gauge_list.append(output_string)         
    
    return obs_rvt_file_path_gauge_list,obs_rvt_file_string_gauge_list,Model_rvt_file_path,Model_rvt_file_string_modify_gauge_list,obsnms
     


def Generate_Raven_Channel_rvp_string_sub(chname,chwd,chdep,chslope,elev,floodn,channeln,iscalmanningn): #writechanel(chname,chwd,chdep,chslope,orchnl,elev,floodn,channeln,iscalmanningn):
    """ Generate string of each subbasin for raven chennel rvp inputs
     
    Function that used to generate a string to define a channel profile 
    for each subbasin in Raven channel rvp input file format.

    Parameters 
    ----------
    chname             :String
        the name of the channel for each SubBasins
        information for this gauge as the station 
        name of this gauge
    chwd              :Float
        channel width
    chdep             :Float
        channel depth 
    chslope           :Float
        channel slope 
    elev              :Float
        channel elevation 
    floodn            :Float
        channel flood plain manning's coefficient 
    channeln          :Float
        main channnel manning's coefficient 
    iscalmanningn     :Bool
        True use channeln or False use 0.035 as main channel
        manning's coefficient 
                                               
    Notes
    ------    
    None 
    
    Returns
    -------
    output_string     : string
        It is the string that contains the content that will be used to 
        to define a channel profile for given subbasin in 
        Raven channel rvp input file format.  
    
    """ 
    output_string_list = []      
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
    output_string_list.append(Channame) #orchnl.write(Channame+"\n")
    Chanslop = "  :Bedslope"+tab+str(chslope)
    output_string_list.append(Chanslop) #orchnl.write(Chanslop+"\n")
    output_string_list.append("  :SurveyPoints") #orchnl.write("  :SurveyPoints"+"\n")
    output_string_list.append("    0"+tab+str(zfld)) #orchnl.write("    0"+tab+str(zfld)+"\n")
    output_string_list.append("    "+str(sidwdfp)+tab+str(elev)) #orchnl.write("    "+str(sidwdfp)+tab+str(elev)+"\n")
    output_string_list.append("    "+str(sidwdfp + 2*chwd)+tab+str(elev)) #orchnl.write("    "+str(sidwdfp + 2*chwd)+tab+str(elev)+"\n")
    output_string_list.append("    "+str(sidwdfp + 2*chwd + sidwd)+tab+str(zbot)) #orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd)+tab+str(zbot)+"\n")
    output_string_list.append("    "+str(sidwdfp + 2*chwd + sidwd + botwd)+tab+str(zbot)) #orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd + botwd)+tab+str(zbot)+"\n")
    output_string_list.append("    "+str(sidwdfp + 2*chwd + 2*sidwd + botwd)+tab+str(elev)) #orchnl.write("    "+str(sidwdfp + 2*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    output_string_list.append("    "+str(sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(elev)) #orchnl.write("    "+str(sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    output_string_list.append("    "+str(2*sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(zfld)) #orchnl.write("    "+str(2*sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(zfld)+"\n")
    output_string_list.append("  :EndSurveyPoints") #orchnl.write("  :EndSurveyPoints"+"\n")
    output_string_list.append("  :RoughnessZones") #orchnl.write("  :RoughnessZones"+"\n")
    output_string_list.append("    0" + tab + str(floodn)) #orchnl.write("    0" + tab + str(floodn) +"\n")
    output_string_list.append("    " + str(sidwdfp + 2*chwd)+ tab + mann) #orchnl.write("    " + str(sidwdfp + 2*chwd)+ tab + mann +"\n")
    output_string_list.append("    " + str(sidwdfp + 2*chwd + 2*sidwd + botwd)+ tab + str(floodn)) #orchnl.write("    " + str(sidwdfp + 2*chwd + 2*sidwd + botwd)+ tab + str(floodn) +"\n")
    output_string_list.append("  :EndRoughnessZones") #orchnl.write("  :EndRoughnessZones"+"\n")
    output_string_list.append(":EndChannelProfile") #orchnl.write(":EndChannelProfile"+"\n")
    output_string_list.append("\n") #orchnl.write("\n")
    output_string_list.append("##############new channel ##############################") #orchnl.write("##############new channel ##############################\n")
    output_string = "\n".join(output_string_list)
    
    return output_string
#########################################################################################################33



def Generate_Raven_Lake_rvh_String(catinfo,Raveinputsfolder,Model_Name):
    """Generate string of raven lake rvh input 
    
    Function that used to generate the content for 
    Raven lake definition Model_Name_Lake.rvh input file.
    
    Parameters 
    ----------    
    catinfo              : DataFrame
        A dataframe includes all attribute for each HRU
        read from polygon shpfile generated by the toolbox 
    Raveinputsfolder     : string
        Folder path and name that save outputs  
        
    Notes
    ------    
    None 
    
    See Also
    --------
    None 
    
    Returns
    -------
    Lake_rvh_string       : string
        It is the string that contains the content that will be used to 
        to define lake parameters for all lakes in 
        Raven lake rvh input file format.
    Lake_rvh_file_path    : string
        It is the string that define the path of
        the raven channel rvp input file.
          
    Examples
    --------
    >>> from WriteRavenInputs import Generate_Raven_Lake_rvh_String
    >>> outFolderraven    = 'c:/path_to_the_raven_input_folder/'
    >>> DataFolder = "C:/Path_to_foldr_of_example_dataset_provided_in_Github_wiki/"
    >>> Model_Folder     = os.path.join(DataFolder,'Model')
    >>> Raveinputsfolder = os.path.join(Model_Folder,'RavenInput')
    >>> finalcatchpath = os.path.join(DataFolder,'finalcat_hru_info.shp')
    >>> tempinfo = Dbf5(finalcatchpath[:-3] + "dbf")
    >>> ncatinfo = tempinfo.to_dataframe()
    >>> Model_Name = 'test'
    >>> ncatinfo2 = ncatinfo.drop_duplicates('HRU_ID', keep='first')
    >>> Lake_rvh_string, Lake_rvh_file_path= Generate_Raven_Lake_rvh_String(ncatinfo2,Raveinputsfolder,lenThres,Model_Name)             
           
    """
    Lake_rvh_file_path   = os.path.join(Raveinputsfolder,Model_Name + "_Lake.rvh")
    Lake_rvh_string_list = []
    tab = '       '
    
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HRU_IsLake'] > 0: ## lake hru
            lakeid = int(catinfo.iloc[i]['HyLakeId'])
            catid = catinfo.iloc[i]['SubId']
            A     = catinfo.iloc[i]['HRU_Area'] ### in meters
            h0 = catinfo.iloc[i]['LakeDepth'] ## m
            WeirCoe = 0.6   
            hruid = int(catinfo.iloc[i]['HRU_ID'])
            Crewd = catinfo.iloc[i]['BkfWidth']  ##3 m 
#            if slakeinfo.iloc[0]['Wshd_area'] < 6000 and slakeinfo.iloc[0]['Wshd_area'] > 0:
        ######write lake information to file
            Lake_rvh_string_list.append(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## ")#f2.write(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## " +"\n")
            Lake_rvh_string_list.append("  :SubBasinID  "+str(int(catid)))#f2.write("  :SubBasinID  "+str(int(catid))+ "\n")
            Lake_rvh_string_list.append("  :HRUID   "+str(int(hruid)))#f2.write("  :HRUID   "+str(int(hruid))+ "\n")
            Lake_rvh_string_list.append("  :Type RESROUTE_STANDARD   ")#f2.write("  :Type RESROUTE_STANDARD   "+"\n")
            Lake_rvh_string_list.append("  :WeirCoefficient  "+str(WeirCoe))#f2.write("  :WeirCoefficient  "+str(WeirCoe)+ "\n")
            Lake_rvh_string_list.append("  :CrestWidth "+str(Crewd))#f2.write("  :CrestWidth "+str(Crewd)+ "\n")
            Lake_rvh_string_list.append("  :MaxDepth "+str(h0))#f2.write("  :MaxDepth "+str(h0)+ "\n")
            Lake_rvh_string_list.append("  :LakeArea    "+str(A))#f2.write("  :LakeArea    "+str(A)+ "\n")
            Lake_rvh_string_list.append(":EndReservoir   ")#f2.write(":EndReservoir   "+"\n")
            Lake_rvh_string_list.append("#############################################")#f2.write("#############################################"+"\n")
            Lake_rvh_string_list.append("###New Lake starts") #f2.write("###New Lake starts"+"\n")
    
    Lake_rvh_string = "\n".join(Lake_rvh_string_list)
    return Lake_rvh_string,Lake_rvh_file_path
    #### write lake input files for different lake zone


def Return_Group_Name_Based_On_Value(value,GroupNames,Group_Thresthold_Values):
    """Return group name 
    It is a function to return group name in GroupNames, based on value and
    Group_Thresthold_Values
    
    Parameters 
    ----------
    value                   : float 
        it is a value that used to identify which group this value belong to
    GroupNames              : list 
        it is a list contain group names 
    Group_Thresthold_Values : list 
        it is a list contains thresthold to define different groups, the value
        of Group_Thresthold_Values should increase from 0 to end.
    Notes
    ------    
    The dimension of the GroupNames should equal to 1 or 
    len(Group_Thresthold_Values) + 1 
    
    See Also
    --------
    None 
    
    Returns
    -------
    GroupName  : String
        the group name determined by value 
    """
    
    ### only one group 
    if len(GroupNames) == 1: 
        GroupName = GroupNames[0]
    elif len(GroupNames) > 1:
        for i in range(0,len(Group_Thresthold_Values)):
            
            ### 
            if value < Group_Thresthold_Values[i]:
                GroupName = GroupNames[i]
                break    
            ## value larger than all thresthold value
            elif i == len(Group_Thresthold_Values) - 1:
                GroupName = GroupNames[i + 1]
    return GroupName
                      
     
def Generate_Raven_Channel_rvp_rvh_String(ocatinfo,Raveinputsfolder,lenThres,iscalmanningn,
                                          Lake_As_Gauge,Model_Name,
                                          SubBasinGroup_NM_Lake,SubBasinGroup_Area_Lake,
                                          SubBasinGroup_NM_Channel,SubBasinGroup_Length_Channel):                                       #Writervhchanl(ocatinfo,Raveinputsfolder,lenThres,iscalmanningn,HRU_ID_NM,HRU_Area_NM,Sub_ID_NM,Lake_As_Gauge = False,Model_Name = 'test'):
    """Generate string of raven chennel rvp input and rvh input
    
    Function that used to generate the content of raven rvh file and 
    channel rvp file. 
    
    Parameters 
    ----------    
    ocatinfo             : DataFrame
        A dataframe includes all attribute for each HRU
        read from polygon shpfile generated by the toolbox 
    Raveinputsfolder     : string
        Folder path and name that save outputs
    lenThres             : float      
        River length threshold; river length smaller than 
        this will write as zero in Raven rvh file
    iscalmanningn        : integer
        If "1", use manning's coefficient in the shpfile table 
        and set to default value (0.035).
        If "-1", do not use manning's coefficients.
    Lake_As_Gauge       : Bool
        If "True", all lake subbasins will labeled as gauged 
        subbasin such that Raven will export lake balance for 
        this lake. If "False", lake subbasin will not be labeled 
        as gauge subbasin.   
    Model_Name          : string
        The Raven model base name. File name of the raven input will be 
        Model_Name.xxx.              
    SubBasinGroup_NM_Channel       : List 
        It is a list of names for subbasin groups, which are grouped based
        on channel length of each subbsin. Should at least has one name
    SubBasinGroup_Length_Channel   : List
        It is a list of float channel length thresthold in meter, to divide
        subbasin into different groups. for example, [1,10,20] will divide
        subbasins into four groups, group 1 with channel length (0,1];
        group 2 with channel length (1,10],
        group 3 with channel length (10,20],
        group 4 with channel length (20,Max channel length].
    SubBasinGroup_NM_Lake          : List 
        It is a list of names for subbasin groups, which are grouped based
        on Lake area of each subbsin. Should at least has one name
    SubBasinGroup_Area_Lake        : List
        It is a list of float lake area thresthold in m2, to divide
        subbasin into different groups. for example, [1,10,20] will divide
        subbasins into four groups, group 1 with lake area (0,1];
        group 2 with lake are (1,10],
        group 3 with lake are (10,20],
        group 4 with lake are (20,Max channel length].
                    
    Notes
    ------    
    None 
    
    See Also
    --------
    Generate_Raven_Channel_rvp_string_sub  : Generate a string to define channel 
                                             for given subbasin in Raven channel
                                             rvp  input format.                                
    Returns
    -------
    Channel_rvp_string       : string
        It is the string that contains the content that will be used to 
        to define channel profiles for all subbasin in 
        Raven channel rvp input file format.
    Channel_rvp_file_path    : string
        It is the string that define the path of
        the raven channel rvp input file.
    Model_rvh_string         : string
        It is the string that contains the content that will be used to 
        to define routing and hru inputs for all subbasin in 
        Raven channel rvh input file format.
    Model_rvh_file_path      : string
        It is the string that define the path of
        the raven channel rvh input file. 
    Model_rvp_string_modify  : string
        It is the string that contains the content that will be used to 
        to modify model rvp file. 
    Model_rvp_file_path      : string
        It is the string that define the path of
        the raven channel rvp input file. 
          
    Examples
    --------
    >>> from WriteRavenInputs import Generate_Raven_Channel_rvp_rvh_String
    >>> outFolderraven    = 'c:/path_to_the_raven_input_folder/'
    >>> DataFolder = "C:/Path_to_foldr_of_example_dataset_provided_in_Github_wiki/"
    >>> Model_Folder     = os.path.join(DataFolder,'Model')
    >>> Raveinputsfolder = os.path.join(Model_Folder,'RavenInput')
    >>> finalcatchpath = os.path.join(DataFolder,'finalcat_hru_info.shp')
    >>> tempinfo = Dbf5(finalcatchpath[:-3] + "dbf")
    >>> ncatinfo = tempinfo.to_dataframe()
    >>> Model_Name = 'test'
    >>> lenThres = 1
    >>> iscalmanningn = -1
    >>> Lake_As_Gauge = -1
    >>> ncatinfo2 = ncatinfo.drop_duplicates('HRU_ID', keep='first')
    >>> Channel_rvp_file_path,Channel_rvp_string,Model_rvh_file_path,Model_rvh_string,Model_rvp_file_path,Model_rvp_string_modify = Generate_Raven_Channel_rvp_rvh_String(ncatinfo2,Raveinputsfolder,lenThres,
    ...                                                                                                                                                                   iscalmanningn,Lake_As_Gauge,Model_Name
    ...                                                                                                                                                                   )
    >>>
        
    """
    Channel_rvp_file_path   = os.path.join(Raveinputsfolder,Model_Name + "_channel.rvp")
    Channel_rvp_string_list = []
    Model_rvh_file_path     = os.path.join(Raveinputsfolder,Model_Name + ".rvh")
    Model_rvh_string_list   = []   
    Model_rvp_file_path     = os.path.join(Raveinputsfolder,Model_Name + ".rvp")
    Model_rvp_string_modify ="\n" +":RedirectToFile " + Model_Name+"_channel.rvp" + "\n"
            
    tab = "     "
                
    catinfo_hru = copy.copy(ocatinfo)
    catinfo     = copy.copy(ocatinfo)
#    print int(catinfo.iloc[0]['SUBID']),len(catinfo.index)

##################3
    Model_rvh_string_list.append("# --------------------------------------------")#orvh.write("# --------------------------------------------"+"\n")
    Model_rvh_string_list.append("# Raven HRU Input file")#orvh.write("# Raven HRU Input file"+"\n")
    Model_rvh_string_list.append("#  lake catchment emulation")#orvh.write("#  lake catchment emulation"+"\n")
    Model_rvh_string_list.append("# --------------------------------------------")#orvh.write("# --------------------------------------------"+"\n")
    Model_rvh_string_list.append(":SubBasins")#orvh.write(":SubBasins"+"\n")
    Model_rvh_string_list.append("  :Attributes   NAME  DOWNSTREAM_ID       PROFILE REACH_LENGTH  GAUGED")#orvh.write("  :Attributes   NAME  DOWNSTREAM_ID       PROFILE REACH_LENGTH  GAUGED"+"\n")
    Model_rvh_string_list.append("  :Units        none           none          none           km    none")#orvh.write("  :Units        none           none          none           km    none"+"\n")
    
    
    catinfo_sub = catinfo.drop_duplicates('SubId', keep='first')### remove duplicated subids, beacuse of including hrus in the dataframe 
    
    SubBasin_Group_Channel = pd.DataFrame(data=np.full((len(catinfo_sub),2),np.nan),columns = ['SubId','SubBasin_Group_NM'])
    SubBasin_Group_Lake    = pd.DataFrame(data=np.full((len(catinfo_sub),2),np.nan),columns = ['SubId','SubBasin_Group_NM'])
    
#    print(catinfo_sub[['SubId','DowSubId']])
    print('Total number of Subbasins are     ' +  str(int((len(catinfo_sub))))+'   '+'SubId')
    for i in range(0,len(catinfo_sub)):
        ### Get catchment width and dpeth
        catid = int(catinfo_sub['SubId'].values[i])
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
#        print(catid,downcatid,len(catinfo_sub.loc[catinfo_sub['SubId'] == downcatid]))
        if catid == downcatid:
            StrDid = str(-1)
        elif len(catinfo_sub.loc[catinfo_sub['SubId'] == downcatid]) == 0:
            StrDid = str(-1)
        else:
            StrDid = str(int(catinfo_sub['DowSubId'].values[i]))
            
        GroupName = Return_Group_Name_Based_On_Value(catinfo_sub['RivLength'].values[i],
                                                     SubBasinGroup_NM_Channel,
                                                     SubBasinGroup_Length_Channel)
        SubBasin_Group_Channel.loc[i,'SubId'] = catid 
        SubBasin_Group_Channel.loc[i,'SubBasin_Group_NM'] = GroupName  
        
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
            
        output_string_chn_rvp_sub = Generate_Raven_Channel_rvp_string_sub(pronam,max(catinfo_sub['BkfWidth'].values[i],1),
                                                max(catinfo_sub['BkfDepth'].values[i],1),
                                                chslope,catinfo_sub['MeanElev'].values[i],
                                                floodn,nchn,iscalmanningn)
                                                
        Channel_rvp_string_list.append(output_string_chn_rvp_sub)
        
        if catinfo_sub['IsObs'].values[i] > 0:
            Guage = '1'
        elif catinfo_sub['IsLake'].values[i] >= 0 and Lake_As_Gauge == True:# and catinfo_sub['HRU_Type'].values[i] == 1: 
            Guage = '1'
        else:
            Guage = '0'
        Model_rvh_string_list.append("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage) #orvh.write("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage+"\n")
        
    Model_rvh_string_list.append(":EndSubBasins")#orvh.write(":EndSubBasins"+"\n")
    Model_rvh_string_list.append("\n")#orvh.write("\n")
##########################################
    Model_rvh_string_list.append(":HRUs")#orvh.write(":HRUs"+"\n")
    Model_rvh_string_list.append("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT")#orvh.write("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT"+"\n")
    Model_rvh_string_list.append("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg")#orvh.write("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg"+"\n")
    
    for i in range(0,len(catinfo_hru.index)):
        
        hruid = int(catinfo_hru['HRU_ID'].values[i])
        catslope = catinfo_hru['HRU_S_mean'].values[i]
        cataspect= catinfo_hru['HRU_A_mean'].values[i]
        
        catarea2 = max(0.0001, catinfo_hru['HRU_Area'].values[i]/1000/1000)  ### in km2
        
        StrGid =  str(hruid) #str( catinfo_hru.iloc[i][HRU_Area_NM])+tab
        
        catid = str(int(catinfo_hru['SubId'].values[i]))+tab
        
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(catinfo_hru['HRU_E_mean'].values[i])+tab
        lat = str(catinfo_hru['HRU_CenY'].values[i])+tab
        lon = str(catinfo_hru['HRU_CenX'].values[i])+tab
        LAND_USE_CLASS = catinfo_hru['LAND_USE_C'].values[i]+tab
        VEG_CLASS = catinfo_hru['VEG_C'].values[i]+tab
        SOIL_PROFILE =catinfo_hru['SOIL_PROF'].values[i]+tab
        AQUIFER_PROFILE ='[NONE]'+tab
        TERRAIN_CLASS ='[NONE]'+tab        

            
        SLOPE = str(catslope)+tab
        ASPECT = str(360 - cataspect)+tab
        Model_rvh_string_list.append("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT) #orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
       
        if catinfo_hru['HRU_IsLake'].values[i] == 1:        
            GroupName = Return_Group_Name_Based_On_Value(catinfo_hru['HRU_Area'].values[i],
                                                         SubBasinGroup_NM_Lake,
                                                         SubBasinGroup_Area_Lake)
            SubBasin_Group_Lake.loc[i,'SubId'] = catinfo_hru['SubId'].values[i] 
            SubBasin_Group_Lake.loc[i,'SubBasin_Group_NM'] = GroupName     
        
    Model_rvh_string_list.append(":EndHRUs")#orvh.write(":EndHRUs"+"\n")
    Model_rvh_string_list.append(":PopulateHRUGroup Lake_HRUs With LANDUSE EQUALS Lake_HRU")#orvh.write(":PopulateHRUGroup Lake_HRUs With LANDUSE EQUALS Lake_HRU" + "\n")
    Model_rvh_string_list.append(":PopulateHRUGroup Land_HRUs With LANDUSE EQUALS FOREST")#orvh.write(":PopulateHRUGroup Land_HRUs With LANDUSE EQUALS FOREST" + "\n")
    Model_rvh_string_list.append(":RedirectToFile " + Model_Name+"_Lake.rvh")#orvh.write(":RedirectToFile TestLake.rvh")

    for i in range(0,len(SubBasinGroup_NM_Channel)):
        Model_rvh_string_list.append(":SubBasinGroup   " + SubBasinGroup_NM_Channel[i])
        SubBasin_Group_Channel_i = SubBasin_Group_Channel.loc[SubBasin_Group_Channel['SubBasin_Group_NM'] == SubBasinGroup_NM_Channel[i]]
        SubIDs_In_Group = SubBasin_Group_Channel_i['SubId'].values
        SubIDs_In_Group_Str_list = ['    ']
        for j in range(0,len(SubIDs_In_Group)):
            SubIDs_In_Group_Str_list.append(str(int(SubIDs_In_Group[i])))
        SubIDs_In_Group_Str = "   ".join(SubIDs_In_Group_Str_list)
         
        Model_rvh_string_list.append(SubIDs_In_Group_Str)
        Model_rvh_string_list.append(":EndSubBasinGroup   ")
        
        Model_rvh_string_list.append("# :SBGroupPropertyOverride " + SubBasinGroup_NM_Channel[i] + "MANNINGS_N 0.001")
        Model_rvh_string_list.append("# :SBGroupPropertyMultiplier "+ SubBasinGroup_NM_Channel[i] + "MANNINGS_N 3")

    for i in range(0,len(SubBasinGroup_NM_Lake)):
        Model_rvh_string_list.append(":SubBasinGroup   " + SubBasinGroup_NM_Lake[i])
        SubBasin_Group_Lake_i = SubBasin_Group_Lake.loc[SubBasin_Group_Lake['SubBasin_Group_NM'] == SubBasinGroup_NM_Lake[i]]
        SubIDs_In_Group = SubBasin_Group_Lake_i['SubId'].values
        SubIDs_In_Group_Str_list = ['    ']
        for j in range(0,len(SubIDs_In_Group)):
            SubIDs_In_Group_Str_list.append(str(int(SubIDs_In_Group[i])))
        SubIDs_In_Group_Str = "   ".join(SubIDs_In_Group_Str_list)
         
        Model_rvh_string_list.append(SubIDs_In_Group_Str)
        Model_rvh_string_list.append(":EndSubBasinGroup   ")
        
        Model_rvh_string_list.append("# :SBGroupPropertyOverride " + SubBasinGroup_NM_Lake[i] + "RESERVOIR_CREST_WIDTH 12.0")
        Model_rvh_string_list.append("# :SBGroupPropertyMultiplier "+ SubBasinGroup_NM_Lake[i] + "RESERVOIR_CREST_WIDTH 1.0")
        
        
    Channel_rvp_string = "\n".join(Channel_rvp_string_list)
    Model_rvh_string   = "\n".join(Model_rvh_string_list)


    return Channel_rvp_file_path,Channel_rvp_string,Model_rvh_file_path,Model_rvh_string,Model_rvp_file_path,Model_rvp_string_modify

