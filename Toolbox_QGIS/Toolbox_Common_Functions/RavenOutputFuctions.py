import pandas as pd 
import os 
import numpy as np

def plotGuagelineobs(scenario,data,outfilename):
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rc('font', family='serif')
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(figsize=(7.48031, 3))
    ax = fig.add_subplot(1, 1, 1)
    colors = np.array(['b','k','g','y','c'])
    for i in range(0,len(scenario)):
        ax.plot(data.index, data[scenario[i]], color=colors[i], 
                ls='solid', linewidth=1,label=scenario[i])
                
    ax.scatter(data.index, data['Obs'],color = 'grey',
               s=0.1,label='Observation')
    plt.legend(loc='upper left',frameon=False,ncol=2,prop={'size': 6})
    plt.ylim(0,max(data[scenario[i]].values)+max(data[scenario[i]].values)*0.1)
    plt.xlabel('Model Time')
    plt.ylabel('Discharge (m3/s)')
    plt.savefig(outfilename, bbox_inches='tight',dpi=300)
    plt.close()  

def plotGuageerror(basename,scenario,data,Diagno):
    for i in range(0,len(scenario)):
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(1, 1, 1)
        results = []
        for k in range(0,len(data)):
            if not math.isnan(data[basename +'_'+scenario[0] +'_obs'].values[k]):
                results.append((data[basename +'_'+scenario[0] +'_obs'].values[k],
                                data[basename +'_'+scenario[i] +'_sim'].values[k],
                                data[basename +'_'+scenario[i] +'_sim'].values[k] -
                                data[basename +'_'+scenario[0] +'_obs'].values[k]))
        results = np.array(results)
        if len(results) <= 0:
            print(scenario[i])
            plt.close()
            continue
        plt.hist(results[:,2], bins='auto')
        plt.savefig('./Figures/'+basename+scenario[i]+'_errhist.pdf', bbox_inches='tight',dpi=300)
        plt.close()
#######
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(1, 1, 1)
        ax.scatter(results[0:len(results)-2,2], results[1:len(results)-1,2],color = 'grey',s=0.1)
        plt.savefig('./Figures/'+basename+scenario[i]+'_errt1t2.pdf', bbox_inches='tight',dpi=300)
        plt.close()
########
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(figsize=(3, 3))
        ax = fig.add_subplot(1, 1, 1)
        plot_acf(results[:,2])
        plt.savefig('./Figures/'+basename+scenario[i]+'_erracf.pdf', bbox_inches='tight',dpi=300)
        plt.close() 

def plotGuageerrorquainy(basename,scenario,data,Diagno,errpars):
    for i in range(0,len(scenario)):
        results = []
        for k in range(0,len(data)):
            if not math.isnan(data[basename +'_'+scenario[0] +'_obs'].values[k]):
                results.append((data[basename +'_'+scenario[0] +'_obs'].values[k],
                                data[basename +'_'+scenario[i] +'_sim'].values[k],
                                data[basename +'_'+scenario[i] +'_sim'].values[k] -
                                data[basename +'_'+scenario[0] +'_obs'].values[k]))
        results = np.array(results)
        if len(results) <= 0:
            continue
        plt.rc('font', family='serif')
        plt.rc('xtick', labelsize='x-small')
        plt.rc('ytick', labelsize='x-small')
        fig = plt.figure(figsize=(7.48031, 3))
        ax = fig.add_subplot(1, 1, 1)
        colors = np.array(['b','k','g','y','c'])
        for j in range(0,len(results)):
            jsim = results[j,1]
            jerror = results[j,2]/(errpars['a'].values[i]*jsim+errpars['b'].values[i])
            nlarge = results[results[:,1] <=jsim]
            ax.scatter(float(len(nlarge))/float(len(results)),jerror,color = 'grey',s=0.1)
        plt.savefig('./Figures/'+basename+scenario[i]+'_errqutiv.pdf', bbox_inches='tight',dpi=300)
        plt.close() 
def NSE(obs,sim):
    import numpy as np
    obsmean     = np.mean(obs)
    diffobsmean = (sim-obsmean)**2
    diffsim     = (sim - obs)**2
    NSE         = 1-     np.sum(diffsim)/np.sum(diffobsmean)  
    return NSE
     
def PlotHydrography_Raven_alone(Path_rvt_Folder = '#',Path_Hydrographs_output_file=['#'],Scenario_NM = ['#','#'],OutputFolder = './'):

    Obs_rvt_NMS = []
    ###obtain obs rvt file name 
    for file in os.listdir(Path_rvt_Folder):
        if file.endswith(".rvt"):
            Obs_rvt_NMS.append(file)
                
    Metric_OUT = pd.DataFrame(Obs_rvt_NMS,columns = ['Obs_NM'])
    Metric_OUT['SubId'] = np.nan
    Metric_OUT['NSE']   = np.nan
    print(Metric_OUT)
    Obs_subids = []                 
    for i in range(0,len(Obs_rvt_NMS)):
        ###find subID
        obs_nm =  Obs_rvt_NMS[i]
        ifilepath = os.path.join(Path_rvt_Folder,obs_nm)
        f = open(ifilepath, "r")
#        print(ifilepath)    
        for line in f:
            firstline_info = line.split()
#            print(firstline_info)
            if firstline_info[0] == ':ObservationData':
                obssubid  = int(firstline_info[2])
                break  ### only read first line
            else:
                obssubid  = -1.2345 
                break     ### only read first line
#        print(obssubid)        
        Obs_subids.append(obssubid)
        Metric_OUT.loc[i,'SubId'] = obssubid
        ## this is not a observation rvt file 
        if obssubid ==-1.2345 :
            continue 
        ####assign column name in the hydrography.csv
        
        colnm_obs = 'sub'+str(obssubid)+' (observed) [m3/s]'
        colnm_sim = 'sub'+str(obssubid)+' [m3/s]'
        colnm_Date = 'date'
        colnm_hr   = 'hour'
            
        ##obtain data from all provided hydrograpy csv output files each hydrograpy csv need has a coorespond scenario name
        Initial_data_frame = 1
        readed_data_correc = 1
        data_len = []
#        print(Metric_OUT)
        for j in range(0,len(Path_Hydrographs_output_file)):
            
            Path_Hydrographs_output_file_j = Path_Hydrographs_output_file[j]
#            print(Path_Hydrographs_output_file_j)#
            i_simresult = pd.read_csv(Path_Hydrographs_output_file[j],sep=',')
            colnames = i_simresult.columns
            ## check if obs name exist in the hydrograpy csv output files
            if colnm_obs in colnames:
                    
                ## Initial lize the reaed in data frame 
                if Initial_data_frame == 1:
                    Readed_Data = i_simresult[[colnm_Date,colnm_hr]]
                    Readed_Data['Obs']  = i_simresult[colnm_obs] 
                    Readed_Data[Scenario_NM[j]] = i_simresult[colnm_sim] 
                    Readed_Data['Date'] = pd.to_datetime(i_simresult[colnm_Date] + ' ' + i_simresult[colnm_hr])
                    Initial_data_frame = -1 
                    data_len.append(len(Readed_Data))
                else:
                    tempdata = i_simresult[[colnm_Date,colnm_hr]]
                    tempdata['Date'] = pd.to_datetime(i_simresult[colnm_Date] + ' ' + i_simresult[colnm_hr])
                    rowmask    = Readed_Data['Date'].isin(tempdata['Date'].values)
                    Readed_Data.loc[rowmask,Scenario_NM[j]] = i_simresult[colnm_sim].values
                    data_len.append(len(tempdata))
            else:
                readed_data_correc = -1
                continue 
        print(readed_data_correc)            
        if readed_data_correc == -1:
            continue 
                
        datalen = min(data_len)
#        Readed_Data = Readed_Data.drop(columns=[colnm_Date,colnm_hr]) 
        Readed_Data = Readed_Data.head(datalen)
        Readed_Data = Readed_Data.set_index('Date') 
            
        Readed_Data = Readed_Data.resample('D').sum()
        Readed_Data['ModelTime'] = Readed_Data.index.strftime('%Y-%m-%d')
    
        Data_NSE    = Readed_Data[Readed_Data['Obs'] > 0]
        Metric_OUT.loc[i,'NSE']   = NSE(Data_NSE['Obs'].values,Data_NSE[Scenario_NM[0]].values)
        print("adfadsfadsfadsf")
        print(Metric_OUT)
        
#        plotGuagelineobs(Scenario_NM,Readed_Data,os.path.join(OutputFolder,obs_nm + '.pdf'))        
def Caluculate_Lake_Active_Depth_and_Lake_Evap(Path_Finalcat_info = '#',Path_ReservoirStages = '#', Path_ReservoirMassBalance = '#',Output_Folder= '#'):
    import pandas as pd 
    import numpy as np 
    from simpledbf import Dbf5 
        
    hyinfocsv         = Path_Finalcat_info[:-3] + "dbf"
    tempinfo          = Dbf5(hyinfocsv)
    finalcat_info     = tempinfo.to_dataframe()
        
    Res_Stage_info    = pd.read_csv(Path_ReservoirStages,sep=',',header = 0)
    Res_MB_info       = pd.read_csv(Path_ReservoirMassBalance,sep=',',header = 0)
    
    Res_Stage_info['Date_2'] = pd.to_datetime(Res_Stage_info['date'] + ' ' + Res_Stage_info['hour'])
    Res_Stage_info = Res_Stage_info.set_index('Date_2')
    
    Res_MB_info['Date_2'] = pd.to_datetime(Res_MB_info['date'] + ' ' + Res_MB_info['hour'])
    Res_MB_info = Res_MB_info.set_index('Date_2')    
    
    
    Res_Wat_Ave_Dep_Vol  = Res_Stage_info.copy(deep=True)
    Res_Wat_Ave_Dep_Vol['Lake_Area'] = np.nan
    Res_Wat_Ave_Dep_Vol['Lake_Stage_Ave'] = np.nan
    Res_Wat_Ave_Dep_Vol['Lake_Vol'] = np.nan
    
    Res_Wat_Ave_Dep_Vol['Lake_Stage_Below_Crest'] = np.nan
    Res_Wat_Ave_Dep_Vol['Lake_Vol_Below_Crest'] = np.nan
    Res_Wat_Ave_Dep_Vol['Lake_Area_Below_Crest'] = np.nan
    
    Col_NMS_Stage  = list(Res_Stage_info.columns)  
    Col_NMS_MB     = list(Res_MB_info.columns)   
    
    finalcat_info_lake_hru     = finalcat_info.loc[(finalcat_info['IsLake'] > 0) & (finalcat_info['HRU_Type'] == 1)]
    
    
    ####
    stage_out_NM  = ['Lake_Id','Lake_Area','Lake_DA','Lake_SubId','#_Day_Active_stage','Min_stage','Max_stage','Ave_stage']
    data_stage = np.full((len(Col_NMS_Stage),8),np.nan)
    Stage_Statis = pd.DataFrame(data = data_stage,columns = stage_out_NM)
    istage = 0
    
    ###
    Evp_out_NM  = ['Year','Total_Lake_Evp_Loss']
    data_evp    = np.full((50,2),0.00000)
    MB_Statis = pd.DataFrame(data = data_evp,columns = Evp_out_NM)
    
    Year_Begin = Res_MB_info.index.min().year
    Year_end = Res_MB_info.index.max().year
    imb = 0
    for iyr in range(Year_Begin,Year_end + 1):
        MB_Statis.loc[imb,'Year'] = iyr
        imb = imb + 1


    for i in range(0,len(Res_Wat_Ave_Dep_Vol)):
        idate = Res_Wat_Ave_Dep_Vol.index[i]
        
        Res_Wat_Ave_Vol_iday = 0.0
        Res_Wat_Ave_Dep_Vol_iday = 0.0
        Res_Wat_Lake_Area = 0.0
        
        Res_Wat_Lake_Area_Below = 0.0 
        Res_Wat_Ave_Vol_Below_iday = 0.0
        Res_Wat_Ave_Dep_Vol_Below_iday = 0.0
        
        for j in range(0,len(finalcat_info_lake_hru)):
            LakeId     = finalcat_info_lake_hru['HyLakeId'].values[j]
            Lake_Area  = finalcat_info_lake_hru['HRU_Area'].values[j] 
            Lake_Subid = finalcat_info_lake_hru['SubId'].values[j] 
            Lake_DA    = finalcat_info_lake_hru['DA'].values[j]             
            Mb_Col_NM       = 'sub'+str(int(Lake_Subid))+' losses [m3]'  
            Stage_Col_NM    = 'sub'+str(int(Lake_Subid)) + ' '  
            if Mb_Col_NM in Col_NMS_MB and Stage_Col_NM in Col_NMS_Stage:
                Res_Wat_Ave_Dep_Vol_iday  = Res_Wat_Ave_Dep_Vol[Stage_Col_NM].values[i] * Lake_Area + Res_Wat_Ave_Dep_Vol_iday
                Res_Wat_Ave_Dep_Lake_Area = Lake_Area 
                
                if Res_Wat_Ave_Dep_Vol[Stage_Col_NM].values[i] < 0:
                    Res_Wat_Ave_Vol_Below_iday = Res_Wat_Ave_Dep_Vol[Stage_Col_NM].values[i] * Lake_Area + Res_Wat_Ave_Vol_Below_iday
                    Res_Wat_Lake_Area_Below =  Lake_Area
#                    print(Res_Wat_Ave_Vol_Below_iday,Res_Wat_Lake_Area_Below,Res_Wat_Ave_Dep_Vol[Stage_Col_NM].values[i])
                
        if Res_Wat_Ave_Dep_Lake_Area > 0:
            Res_Wat_Ave_Dep_iday =  Res_Wat_Ave_Dep_Vol_iday/Res_Wat_Ave_Dep_Lake_Area
        else:
            Res_Wat_Ave_Dep_iday = np.nan 

        if Res_Wat_Lake_Area_Below > 0:
            Res_Wat_Ave_Dep_Vol_Below_iday =  Res_Wat_Ave_Vol_Below_iday/Res_Wat_Lake_Area_Below
        else:
            Res_Wat_Ave_Dep_Vol_Below_iday = np.nan
            
                    
        Res_Wat_Ave_Dep_Vol.loc[idate,'Lake_Area'] = Res_Wat_Ave_Dep_Lake_Area
        Res_Wat_Ave_Dep_Vol.loc[idate,'Lake_Stage_Ave'] = Res_Wat_Ave_Dep_iday
        Res_Wat_Ave_Dep_Vol.loc[idate,'Lake_Vol'] = Res_Wat_Ave_Dep_Vol_iday 
    
        Res_Wat_Ave_Dep_Vol.loc[idate,'Lake_Area_Below_Crest'] = Res_Wat_Lake_Area_Below
        Res_Wat_Ave_Dep_Vol.loc[idate,'Lake_Stage_Below_Crest'] = Res_Wat_Ave_Dep_Vol_Below_iday
        Res_Wat_Ave_Dep_Vol.loc[idate,'Lake_Vol_Below_Crest'] = Res_Wat_Ave_Vol_Below_iday 
        
                     
    
    Res_Wat_Ave_Dep_Vol = Res_Wat_Ave_Dep_Vol[['date','hour','Lake_Area','Lake_Stage_Ave','Lake_Vol','Lake_Area_Below_Crest','Lake_Stage_Below_Crest','Lake_Vol_Below_Crest']]
    
    Res_Wat_Ave_Dep_Vol.to_csv(os.path.join(Output_Folder,'Lake_Wat_Ave_Depth_Vol.csv'),sep = ',')
                
    for i in range(0,len(finalcat_info_lake_hru)):
        LakeId     = finalcat_info_lake_hru['HyLakeId'].values[i]
        Lake_Area  = finalcat_info_lake_hru['HRU_Area'].values[i] 
        Lake_Subid = finalcat_info_lake_hru['SubId'].values[i] 
        Lake_DA    = finalcat_info_lake_hru['DA'].values[i] 
        
        ####
        stage_idx =  Stage_Statis.index[istage]
        Stage_Statis.loc[stage_idx,'Lake_Id']  = LakeId
        Stage_Statis.loc[stage_idx,'Lake_Area']  = Lake_Area
        Stage_Statis.loc[stage_idx,'Lake_SubId'] = Lake_Subid
        Stage_Statis.loc[stage_idx,'Lake_DA'] = Lake_DA
        istage = istage + 1 
    
        ###            
        Mb_Col_NM       = 'sub'+str(int(Lake_Subid))+' losses [m3]'  
        Stage_Col_NM    = 'sub'+str(int(Lake_Subid)) + ' '
        
        if Mb_Col_NM not in Col_NMS_MB or Stage_Col_NM not in Col_NMS_Stage:
            print(Mb_Col_NM in Col_NMS_MB,Mb_Col_NM[0] == Col_NMS_MB[7][0],len(Mb_Col_NM),len(Mb_Col_NM[7]),Mb_Col_NM,Col_NMS_MB[7])
            print(Stage_Col_NM in Col_NMS_Stage,Stage_Col_NM[0] == Col_NMS_Stage[7][0],len(Stage_Col_NM),len(Col_NMS_Stage[7]),Stage_Col_NM,Col_NMS_Stage[7])
        
        if Mb_Col_NM in Col_NMS_MB and Stage_Col_NM in Col_NMS_Stage:
            Mb_info_lake     = Res_MB_info[Mb_Col_NM]
            Stage_info_lake  = Res_Stage_info[Stage_Col_NM]
            
            ### For stage
            Stage_Statis = Calulate_Yearly_Reservior_stage_statistics(Stage_info_lake,Stage_Statis,stage_idx,Stage_Col_NM)
            
            ### for Mass 
            
            for iyr in range(Year_Begin,Year_end+1):
                Mb_info_lake_iyr        = Mb_info_lake.loc[Stage_info_lake.index.year == iyr].values
                MB_Statis.loc[MB_Statis['Year'] == iyr,'Total_Lake_Evp_Loss'] = MB_Statis.loc[MB_Statis['Year'] == iyr,'Total_Lake_Evp_Loss'] + sum(Mb_info_lake_iyr)
                
    Stage_Statis.to_csv(os.path.join(Output_Folder,'Lake_Stage_Yearly_statistics.csv'),sep = ',')
    MB_Statis.to_csv(os.path.join(Output_Folder,'Lake_MB_Yearly_statistics.csv'),sep = ',')
            
        #### 
def Calulate_Yearly_Reservior_stage_statistics(Stage_info_lake,Stage_Statis,stage_idx,Stage_Col_NM):
    Year_Begin = Stage_info_lake.index.min().year
    Year_end = Stage_info_lake.index.max().year
    
    
    Num_Day_Active_stage_sum = 0 
    Min_stage_sum = 0
    Max_stage_sum = 0
    Ave_stage_sum = 0
    
    nyear = 0
    for iyr in range(Year_Begin,Year_end+1):
        nyear = nyear + 1
        Stage_info_lake_iyr        = Stage_info_lake.loc[Stage_info_lake.index.year == iyr].values
        Active_Stage_info_lake_iyr = Stage_info_lake_iyr[Stage_info_lake_iyr < 0]        
        if(len(Active_Stage_info_lake_iyr) > 0):
            Num_Day_Active_stage_sum = Num_Day_Active_stage_sum + len(Active_Stage_info_lake_iyr)
            Min_stage_sum      = Min_stage_sum + min(Stage_info_lake_iyr)
            Max_stage_sum      = Max_stage_sum + max(Stage_info_lake_iyr)
            Ave_stage_sum      = Ave_stage_sum + np.average(Stage_info_lake_iyr)
        
      
        
    Stage_Statis.loc[stage_idx,'#_Day_Active_stage'] = Num_Day_Active_stage_sum/nyear
    Stage_Statis.loc[stage_idx,'Min_stage'] = Min_stage_sum/nyear
    Stage_Statis.loc[stage_idx,'Max_stage'] = Max_stage_sum/nyear
    Stage_Statis.loc[stage_idx,'Ave_stage'] = Ave_stage_sum/nyear
    
    return Stage_Statis    
        

#        'sub800883 losses [m3]'  'sub800883 losses [m3]'