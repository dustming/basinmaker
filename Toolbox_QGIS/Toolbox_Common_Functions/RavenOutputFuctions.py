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
    import pandas as pd 
    import os 
    import numpy as np
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
#            print(Path_Hydrographs_output_file_j)
            i_simresult = pd.read_csv(Path_Hydrographs_output_file[j],sep=',')
            colnames = i_simresult.columns
            print(colnames)
            print(colnm_obs)
            print(colnm_obs in colnames)
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
#    def Caluculate_Lake_Active_Depth_and_Lake_Evap(Path_Finalcat_info = '#',Path_ReservoirStages = '#', Path_ReservoirMassBalance = '#'):
        
        