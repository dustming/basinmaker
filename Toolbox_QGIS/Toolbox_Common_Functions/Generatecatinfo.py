
import numpy as np
from GetBasinoutlet import Getbasinoutlet,Nextcell
from Calculate_River_Len_Slope import Getcatrivlenslope_hydroshed
import copy



def Streamorderanddrainagearea(catinfo):
    catlist = np.full((len(catinfo)), -9)
    icat = 0
    iseg = 1
    ### find first segments of all reaches, no upstream reaches 
    for i in range(0,len(catinfo)):
        idx = catinfo.index[i]
        if catinfo['SubId'].values[i] == catinfo['DowSubId'].values[i]:
            catinfo.loc[idx,'DowSubId'] = -1
        catid = catinfo['SubId'].values[i]
        if len(catinfo[catinfo['DowSubId'] == catid]) == 0: ### the river seg has no upstream segment 
            catlist[icat] = int(catinfo['DowSubId'].values[i])   #### store next reach segment 
            catinfo.loc[idx,'DA'] = catinfo['BasArea'].values[i]
            catinfo.loc[idx,'Strahler'] = 1
            catinfo.loc[idx,'Seg_order'] = 1
            catinfo.loc[idx,'Seg_ID'] = iseg
            icat = icat + 1
            iseg = iseg +1
            
    catlist = np.unique(catlist)
    catlist = catlist[catlist > 0]
#    print(catlist)
    ### Loop for each first reach, until go to reach intersection 
    newcatlist = np.full((len(catinfo)), -9)
    inewstart = 0
    
    for i in range(0,len(catlist)):
        catid = catlist[i]
        F_intersect = 1 
#        print("new start            ",i,catid)
        while F_intersect == 1 and catid > 0:
            Up_Reaches_info = catinfo[catinfo['DowSubId'] == catid]
            cur_Reach_info = catinfo[catinfo['SubId'] == catid]
            curcat_idx = catinfo['SubId'] == catid

            if len(Up_Reaches_info) == 1:   ### only have one upstream 
                catinfo.loc[curcat_idx,'DA'] = cur_Reach_info['BasArea'].values[0] + Up_Reaches_info['DA'].values[0]
                catinfo.loc[curcat_idx,'Strahler'] = Up_Reaches_info['Strahler'].values[0]
                catinfo.loc[curcat_idx,'Seg_order'] = Up_Reaches_info['Seg_order'].values[0] + 1
                catinfo.loc[curcat_idx,'Seg_ID'] = Up_Reaches_info['Seg_ID'].values[0]
#                print('1',catid,catinfo.loc[curcat_idx,'DA'].values,catinfo.loc[curcat_idx,'Strahler'].values,catinfo.loc[curcat_idx,'Sub_order'].values)
                catid =  int(cur_Reach_info['DowSubId'].values[0])
            else:  ### has mutiple upstram 
                if np.min(Up_Reaches_info['Strahler'].values) > 0: ### all upstream has been processed 
                    catinfo.loc[catinfo['SubId'] == catid,'DA'] = cur_Reach_info['BasArea'].values[0] + np.sum(Up_Reaches_info['DA'].values)
                    if np.min(Up_Reaches_info['Strahler'].values) == np.max(Up_Reaches_info['Strahler'].values): ### two reach has the same order
                        catinfo.loc[curcat_idx,'Strahler'] = Up_Reaches_info['Strahler'].values[0] + 1
                        catinfo.loc[curcat_idx,'Seg_order'] = 1
                        catinfo.loc[curcat_idx,'Seg_ID'] = iseg +1
                        iseg = iseg +1
#                        print('2',catid,catinfo.loc[catinfo['SubId'] == catid,'DA'].values,catinfo.loc[catinfo['SubId'] == catid,'Strahler'].values,catinfo.loc[catinfo['SubId'] == catid,'Sub_order'].values)
                    else:
                        max_straorder = np.max(Up_Reaches_info['Strahler'].values)
                        catinfo.loc[curcat_idx,'Strahler']  = max_straorder
                        catinfo.loc[curcat_idx,'Seg_order'] = 1
                        catinfo.loc[curcat_idx,'Seg_ID'] = iseg +1
                        iseg = iseg +1
#                        print('3',catid,catinfo.loc[catinfo['SubId'] == catid,'DA'].values,catinfo.loc[catinfo['SubId'] == catid,'Strahler'].values,catinfo.loc[catinfo['SubId'] == catid,'Sub_order'].values) 
                    catid =  int(cur_Reach_info['DowSubId'].values[0])
                else:  ## there are some reach has not been processed, save id to the list and wait for another loob 
                    newcatlist[inewstart]  =  int(catid)
                    inewstart = inewstart + 1  
                    F_intersect = 0               

    newcatlist = np.unique(newcatlist)
    newcatlist = newcatlist[newcatlist>0]   
#####################################
    riv_segs = catinfo['Seg_ID'].values
    riv_segs = np.unique(riv_segs)
    riv_segs = riv_segs[riv_segs>0]
    for i in range(0,len(riv_segs)):
        iriv_seg = riv_segs[i]
        catinfo_saseg = catinfo[catinfo['Seg_ID'] == iriv_seg]
        segs_nobkf_idx = catinfo_saseg['BkfWidth'] < 0  ### index that do not have bankfull width data 
        segs_bkf_idx = catinfo_saseg['BkfWidth'] > 0 
        segs_norivslop_idx = catinfo_saseg['RivSlope'] < 0  ### index that do not have bankfull width data 
        if len(catinfo_saseg[segs_nobkf_idx]) > 0 or len(catinfo_saseg[segs_norivslop_idx]) > 0:      ####bankfulll width data not avaiable for some of 
            nobkfdatarives = catinfo_saseg[segs_nobkf_idx]
            bkfdatarives = catinfo_saseg[segs_bkf_idx]
            norivsloperives = catinfo_saseg[segs_norivslop_idx]

            seg_rivslope_ave = -1.2345
            seg_bkfwidth_ave = -1.2345
            seg_bkfdepth_ave = -1.2345
            seg_bkfqmean_ave = -1.2345

            seg_max_dems = catinfo_saseg['Max_DEM'].values
            seg_min_dems = catinfo_saseg['Min_DEM'].values
            if len(seg_max_dems[seg_max_dems > -1000]) > 0:
                seg_max_dem = np.max(seg_max_dems[seg_max_dems > -1000])
                seg_min_dem = np.max(seg_min_dems[seg_min_dems > -1000])
            else:
                continue
            
            if (seg_max_dem > seg_min_dem and np.sum(catinfo_saseg['RivLength'].values > 0)):
                seg_rivslope_ave = (seg_max_dem - seg_min_dem)/np.sum(catinfo_saseg['RivLength'].values)
            else:
                seg_rivslope_ave = 0.00012345

            if len(bkfdatarives) > 0:
                seg_bkfwidth_ave = np.average(bkfdatarives['BkfWidth'].values,weights = bkfdatarives['RivLength'].values)
                seg_bkfdepth_ave = np.average(bkfdatarives['BkfDepth'].values,weights = bkfdatarives['RivLength'].values)
                seg_bkfqmean_ave = np.average(bkfdatarives['Q_Mean'].values,weights = bkfdatarives['RivLength'].values)
            ### calcuate river 
            for i in range(0,len(nobkfdatarives)):
                icat = nobkfdatarives['SubId'].vaues[i]
                icat_idx = catinfo['SubId'] == icat
                catinfo.loc[icat_idx,'BkfWidth']  = seg_bkfwidth_ave
                catinfo.loc[icat_idx,'BkfDepth']  = seg_bkfdepth_ave
                catinfo.loc[icat_idx,'Q_Mean']  = seg_bkfqmean_ave

            for i in range(0,len(norivsloperives)):
                icat = norivsloperives['SubId'].values[i]
                icat_idx = catinfo['SubId'] == icat
                catinfo.loc[icat_idx,'RivSlope']  = seg_rivslope_ave           
                
    ### calcuate channel manning's coefficient     
    for i in range(0,len(catinfo)):
        idx =  catinfo.index[i]
        if catinfo['BkfWidth'].values[i] > 0 and catinfo['RivSlope'].values[i] > 0 :
            catinfo.loc[idx,'Ch_n'] = calculateChannaln(catinfo['BkfWidth'].values[i],catinfo['BkfDepth'].values[i],
                              catinfo['Q_Mean'].values[i],catinfo['RivSlope'].values[i])
        if catinfo['IsObs'].values[i] > 0:
            if catinfo['DA_Obs'].values[i] >0:
                catinfo.loc[idx,'DA_error'] = (catinfo['DA'].values[i]/1000.0/1000.0 - catinfo['DA_Obs'].values[i])/catinfo['DA_Obs'].values[i]
    return catinfo 



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
    if (V > 0):
        n = (Rch**(2.0/3.0))*(slope**(1.0/2.0))/V
    else:
        n = 0.0012345
    return n



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



def Generatecatinfo(Watseds,fac,fdir,lake,dem,area,catinfo,allcatid,lakeinfo,width,depth,
                    rivlen,obs,nrows,ncols,slope,landuse,landuseinfo,Q_Mean):
    finalcat = copy.copy(Watseds)
    rivpath = copy.copy(Watseds)
    rivpath[:,:] = -9999
    for i in range(0,len(allcatid)):
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
        catrivlen,catrivslp,catrivslp2,rivpath = Getcatrivlenslope_hydroshed(rowcol[:,0],rowcol[:,1],rivlen,dem,fac,fdir,finalcat,
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
            catinfo.loc[i,'RivSlope'] = 0.0001234
        else:
            catinfo.loc[i,'RivSlope'] = catrivslp
#        catinfo[i,27] = catrivslp2
        if len(slopet) < 1:
            catinfo.loc[i,'RivSlope'] = 0.0001234
            catinfo.loc[i,'BasinSlope'] = 0.0001234
        catinfo.loc[i,'FloodP_n'] = Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo)
        catinfo.loc[i,'Q_Mean'] = catQ
#    writeraster(outputFolder + "/" + "rivpath.asc",rivpath,OutputFolder + "/" + "dir")
    return catinfo,rivpath


     

def Generatecatinfo_riv(Watseds,fac,fdir,lake,dem,catinfo,allcatid,width,depth,
                    obs,slope,aspect,landuse,slop_deg,Q_Mean,netcat,landuseinfo,lakeinfo,
                    nrows,ncols,leninfo,areainfo,obsinfo,NonConcLakeInfo,NonCL_array):
    finalcat = copy.copy(Watseds)
    for i in range(0,len(allcatid)):
        catid = allcatid[i].astype(int)
        catmask2 = netcat ==catid
        catinfo.loc[i,'SubId'] = catid
        catmask = Watseds == catid
        trow,tcol = Getbasinoutlet(catid,finalcat,fac,fdir,nrows,ncols)
        k = 1
        ttrow,ttcol = trow,tcol
        while catinfo['DowSubId'].values[i] < 0 and k < 20:
            nrow,ncol = Nextcell(fdir,ttrow,ttcol)### get the downstream catchment id
            if nrow < 0 or ncol < 0:
                catinfo.loc[i,'DowSubId'] = -1
                break;
            elif nrow >= nrows or ncol >= ncols:
                catinfo.loc[i,'DowSubId'] = -1
                break;
            elif finalcat[nrow,ncol] <= 0:
                catinfo.loc[i,'DowSubId'] = -1
            else:
                catinfo.loc[i,'DowSubId'] = finalcat[nrow,ncol]
            k = k + 1
            ttrow = nrow
            ttcol = ncol 
################################## Get lake information        
        lakeinriv = lake[catmask]
        lakeids = np.unique(lakeinriv)
        lakeids = lakeids[lakeids > 0]
        if len(lakeids) == 1:
            lakeid = lakeids[0]
        elif len(lakeids) > 1:
            print('Warning:  stream    ',catid,'connected with ',len(lakeids),'   Lakes')
            lakeid = lakeids[0]
            for j in range(1,len(lakeids)):
                if len(np.argwhere(lakeinriv == lakeid)) < len(np.argwhere(lakeinriv == lakeids[j])):
                    lakeid = lakeids[j]
        else:
            lakeid = -1 
            
        if lakeid > 0:
            slakeinfo = lakeinfo.loc[lakeinfo['Hylak_id'] == lakeid]
            catinfo.loc[i,'IsLake'] = 1
            catinfo.loc[i,'HyLakeId'] = lakeid
            catinfo.loc[i,'LakeVol'] = slakeinfo.iloc[0]['Vol_total']
            catinfo.loc[i,'LakeArea']= slakeinfo.iloc[0]['Lake_area']
            catinfo.loc[i,'LakeDepth']= slakeinfo.iloc[0]['Depth_avg']
            catinfo.loc[i,'Laketype'] = slakeinfo.iloc[0]['Lake_type']
########Check if it is observation points
#        print(catid,obs[trow,tcol],fac[trow,tcol],fac[nrow,ncol],finalcat[trow,tcol],finalcat[nrow,ncol])
        if obs[trow,tcol]  > 0:
#            arcpy.AddMessage(str(catid)+"      "+str(obs[trow,tcol]))
            catinfo.loc[i,'IsObs'] =  obs[trow,tcol]
            obsid = float(obs[trow,tcol])
            catinfo.loc[i,'DA_Obs']  = obsinfo.loc[obsinfo['Obs_ID'] == obsid]['DA_obs'].values[0]
            catinfo.loc[i,'Obs_NM']  = obsinfo.loc[obsinfo['Obs_ID'] == obsid]['STATION_NU'].values[0]
            catinfo.loc[i,'SRC_obs'] = obsinfo.loc[obsinfo['Obs_ID'] == obsid]['SRC_obs'].values[0]
#            print(obsinfo.loc[obsinfo['Obs_ID'] == obsid]['SRC_obs'])
########Slopes slope,aspect,landuse,slop_deg
        slopeinriv = slope[catmask2]
        aspectinriv = aspect[catmask2]
        slop_deginriv = slop_deg[catmask2]
        deminriv = dem[catmask]
        deminriv2 = dem[catmask2]
        
        if(len(slop_deginriv[slop_deginriv >= 0])) > 0:
            slop_deginriv[slop_deginriv <0] = np.NaN  
            catinfo.loc[i,'BasSlope'] = np.maximum(np.nanmean(slop_deginriv),0.1)
        else:
            catinfo.loc[i,'BasSlope'] = -1.2345      

        if(len(aspectinriv[aspectinriv >= 0])) > 0:
            aspectinriv[aspectinriv <0] = np.NaN  
            catinfo.loc[i,'BasAspect'] = np.maximum(np.nanmean(aspectinriv),0.1)
        else:
            catinfo.loc[i,'BasAspect'] = -1.2345  
        
        if(len(deminriv[deminriv > 0])) > 0:
            deminriv[deminriv <=0] = np.NaN
            maxdem = np.nanmax(deminriv)
            mindem = np.nanmin(deminriv)
            catinfo.loc[i,'Min_DEM'] = mindem
            catinfo.loc[i,'Max_DEM'] = maxdem
        else:
            maxdem = -1.2345
            mindem = -1.2345

        if(len(deminriv2[deminriv2 > 0])) > 0:
            deminriv2[deminriv2 <=0] = np.NaN  
            catinfo.loc[i,'MeanElev'] = np.nanmean(deminriv2)
        else:
            catinfo.loc[i,'MeanElev'] =-1.2345
                    
            
        rivlen = np.unique(leninfo.loc[leninfo['Gridcode'] == catid]['Length_m'].values)  #'Area_m'
        if len(rivlen) == 1:
            catinfo.loc[i,'RivLength'] = rivlen
            if rivlen >= 0:
                if max(0,float((maxdem - mindem))/float(rivlen)) == 0:
                    catinfo.loc[i,'RivSlope'] =-1.2345
                else:
                    catinfo.loc[i,'RivSlope'] = max(0,float((maxdem - mindem))/float(rivlen))
            else:
                catinfo.loc[i,'RivSlope'] = -1.2345
        else:
            print("Warning  river length of stream  " , catid, "   need check   ", len(rivlen) )
            catinfo.loc[i,'RivLength'] = -1.2345
            catinfo.loc[i,'RivSlope'] = -1.2345
            
            
        catarea = np.unique(areainfo.loc[areainfo['Gridcode'] == catid]['Area_m'].values)  #'Area_m'
        if len(catarea) == 1:
            catinfo.loc[i,'BasArea'] = catarea
        else:
            print("Warning  basin area of stream  " , catid, "   need check   ", len(catarea) )
            catinfo.loc[i,'BasArea'] = -1.2345

##########
        Landtypes = landuse[catmask]
        Landtypeid = np.unique(Landtypes)
        Landtypeid1 = Landtypeid[Landtypeid >= 0]
        Landtypeid2 = Landtypeid1[Landtypeid1 > 0]
        Landtypes = Landtypes[Landtypes > 0]
        if len(Landtypes) > 0 and float(len(Landtypeid2))/float(len(Landtypeid1)) >= 0.1:
            sum = 0.0
            for j in range(0,len(Landtypeid2)):
                iid = Landtypeid2[j]
                sum = sum + landuseinfo[landuseinfo['RasterV'] == iid]['MannV'].values*len(np.argwhere(Landtypes == iid))
            floodn = sum/(len(Landtypes))
        else:
            floodn = 0.035
        catinfo.loc[i,'FloodP_n'] = floodn
            
########Got basin width and depth
        widthinriv = width[catmask]
        depthinriv = depth[catmask]
        Q_Meaninriv = Q_Mean[catmask]
        
        widthids = np.unique(widthinriv)
        widthids = widthids[widthids > 0]

        if(len(widthids)) > 0:
            widthinriv[widthinriv <=0] = np.NaN
            depthinriv[depthinriv <=0] = np.NaN
            Q_Meaninriv[Q_Meaninriv <=0] = np.NaN
            catinfo.loc[i,'BkfWidth'] = np.nanmean(widthinriv)
            catinfo.loc[i,'BkfDepth'] = np.nanmean(depthinriv)
            catinfo.loc[i,'Q_Mean'] = np.nanmean(Q_Meaninriv)
        else:
            catinfo.loc[i,'BkfWidth'] = -1.2345
            catinfo.loc[i,'BkfDepth'] = -1.2345
            catinfo.loc[i,'Q_Mean'] =  -1.2345   
             
######## Non connected lakes
        catinfo.loc[i,'NonLDArea'] =  0.0
        Nonlakes = NonCL_array[catmask2]
        Nonlakes = Nonlakes[Nonlakes > 0]
        Nonlakes = np.unique(Nonlakes)
        Daarea = 0.0
        if len(Nonlakes) > 0:
            for inl in range(0,len(Nonlakes)):
                NLid = int(Nonlakes[inl])
                nlakeinfoidx = NonConcLakeInfo['Gridcode'] == NLid
                SubID_NL = -1
                ### assign NONconnectlake attributes 
                slakeinfo = lakeinfo.loc[lakeinfo['Hylak_id'] == NLid]
                NonConcLakeInfo.loc[nlakeinfoidx,'HyLakeId']  = NLid
                NonConcLakeInfo.loc[nlakeinfoidx,'LakeVol']   = slakeinfo.iloc[0]['Vol_total']
                NonConcLakeInfo.loc[nlakeinfoidx,'LakeArea']  = slakeinfo.iloc[0]['Lake_area']
                NonConcLakeInfo.loc[nlakeinfoidx,'LakeDepth'] = slakeinfo.iloc[0]['Depth_avg']
                NonConcLakeInfo.loc[nlakeinfoidx,'Laketype']  = slakeinfo.iloc[0]['Lake_type']
               ### find non connect lake drainage subid 
                trow,tcol = Getbasinoutlet(NLid,NonCL_array,fac,fdir,nrows,ncols)
                k = 1
                ttrow,ttcol = trow,tcol
                while SubID_NL < 0 and k < 20:
                    nrow,ncol = Nextcell(fdir,ttrow,ttcol)### get the downstream catchment id
                    if nrow < 0 or ncol < 0:
                        SubID_NL = -1
                        break;
                    elif nrow >= nrows or ncol >= ncols:
                        SubID_NL = -1
                        break;
                    elif NonCL_array[nrow,ncol] == NLid:
                        SubID_NL = -1
                    else:
                        SubID_NL = netcat[nrow,ncol]
                    k = k + 1
                    ttrow = nrow
                    ttcol = ncol 
                if SubID_NL == catid: ### the non connected lake drainage to this catid 
                    Daarea = Daarea + NonConcLakeInfo.loc[nlakeinfoidx]['Area_m'].values[0]
                NonConcLakeInfo.loc[nlakeinfoidx,'SubId_riv'] = SubID_NL
            catinfo.loc[i,'NonLDArea'] =  Daarea
    catinfo = Writecatinfotodbf(catinfo)
    return catinfo,NonConcLakeInfo



    
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

