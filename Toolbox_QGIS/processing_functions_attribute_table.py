import numpy as np
import pandas as pd 
import copy 

def Calculate_Longest_flowpath(mainriv_merg_info):
    mainriv_merg_info_sort = mainriv_merg_info.sort_values(["DA"], ascending = (False))
#    print(mainriv_merg_info_sort[['SubId','DowSubId','DA','Strahler','RivLength']])
    longest_flow_pathes = np.full(100,0)
#    print(longest_flow_pathes)
    npath = 1

    #### loop upstream to find longest flow path
    Pathid  = np.full(1000,-1)
    subid = mainriv_merg_info_sort['SubId'].values[0]
    npath_current = 1
    Pathid[npath - 1] = subid
#    print('####################################################################################')
    while len(Pathid[Pathid > 0]) > 0:
        nPathid  =  np.full(1000,-1)
        npath    = npath_current

#        print('###################################')
#        print(npath,Pathid[0:npath])
#        print('###################################')
        for ipath in range(0,npath_current):
            c_subid_ipath = Pathid[ipath]

            if c_subid_ipath < 0:  ### means this path has been closed due to no more subbasin within the lake domain
                continue

            longest_flow_pathes[ipath] = mainriv_merg_info_sort.loc[mainriv_merg_info_sort['SubId'] == c_subid_ipath,'RivLength'] +longest_flow_pathes[ipath] ## add river length to current path
            Strahler_order_ipath = mainriv_merg_info_sort.loc[mainriv_merg_info_sort['SubId'] == c_subid_ipath,'Strahler'].values[0]

            upstream_sub_infos = mainriv_merg_info_sort.loc[mainriv_merg_info_sort['DowSubId'] == c_subid_ipath] ## get upstream info

            if len(upstream_sub_infos) <= 0: ## no more upstream catchment within the domain of the lake
#                print("path        closed        ",ipath)
                continue

            ## look for upstream catchment has the same upstream_sub_infos_eq_Strahler first
#            print(Strahler_order_ipath)
#            print(upstream_sub_infos['Strahler'])
            upstream_sub_infos_eq_Strahler = upstream_sub_infos.loc[upstream_sub_infos['Strahler'] == Strahler_order_ipath]

            if len(upstream_sub_infos_eq_Strahler) > 0: ### has a upstream river has the saem strahler id, no new path will be added
                nPathid[ipath] = upstream_sub_infos_eq_Strahler['SubId'].values[0] ### add this upstream id to nPathid
                continue
            else:
                upstream_sub_infos_eq_Strahler_1 = upstream_sub_infos.loc[upstream_sub_infos['Strahler'] == Strahler_order_ipath - 1]

                for inpath in range(0,len(upstream_sub_infos_eq_Strahler_1)):
                    ### this brance sperate into two or several reaches, the starting river length for all of them are the same
                    if inpath == 0:
                        nPathid[ipath] = upstream_sub_infos_eq_Strahler_1['SubId'].values[inpath]
#                        print(nPathid[ipath],ipath,upstream_sub_infos_eq_Strahler_1['SubId'].values[inpath],'aaaaa',range(0,len(upstream_sub_infos_eq_Strahler_1)))
                    else:
                        nPathid[npath + 1 - 1] = upstream_sub_infos_eq_Strahler_1['SubId'].values[inpath]
                        longest_flow_pathes[npath + 1 - 1] = longest_flow_pathes[ipath]
#                        print(npath + 1 - 1,longest_flow_pathes[npath + 1 - 1],nPathid[npath + 1 - 1],'bbbbb',range(0,len(upstream_sub_infos_eq_Strahler_1)))
                        npath = npath + 1



        Pathid = nPathid
        npath_current = npath
    Longestpath = max(longest_flow_pathes)

    return Longestpath
    
    
    
def New_SubId_To_Dissolve(subid,catchmentinfo,mapoldnew_info,upsubid = -1,ismodifids = -1,modifiidin = [-1],mainriv = [-1],Islake = -1,seg_order = -1):
    sub_colnm = 'SubId'
    routing_info      = catchmentinfo[['SubId','DowSubId']].astype('float').values
    if ismodifids < 0:
        Modify_subids1            = Defcat(routing_info,subid)   ### find all subids drainage to this subid
        if upsubid > 0:
            Modify_subids2        = Defcat(routing_info,upsubid)
            mask = np.in1d(Modify_subids1, Modify_subids2)
            Modify_subids = Modify_subids1[np.logical_not(mask)]
        else:
            Modify_subids =  Modify_subids1

    else:
        Modify_subids =  modifiidin
#    print("##########################################")
#    print(subid)
#    print(Modify_subids)
    cbranch                  = catchmentinfo[catchmentinfo[sub_colnm].isin(Modify_subids)].copy()
    tarinfo                  = catchmentinfo[catchmentinfo[sub_colnm] == subid].copy()   ### define these subs attributes
    ### average river slope info

    mainriv_merg_info = mainriv.loc[mainriv['SubId'].isin(Modify_subids)].copy()
    mainriv_merg_info = mainriv_merg_info.loc[mainriv_merg_info['RivLength'] > 0].copy()
    idx = tarinfo.index[0]
#    print(tarinfo.loc[idx,'BasArea'],"1")
    if len(mainriv_merg_info) > 0:
        tarinfo.loc[idx,'RivLength'] = np.sum(mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'RivSlope']  = np.average(mainriv_merg_info['RivSlope'].values ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'FloodP_n']  = np.average(mainriv_merg_info['FloodP_n'].values ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'Q_Mean']    = np.average(mainriv_merg_info['Q_Mean'].values   ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'Ch_n']      = np.average(mainriv_merg_info['Ch_n'].values     ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'BkfWidth']  = np.max(mainriv_merg_info['BkfWidth'].values)
        tarinfo.loc[idx,'BkfDepth']  = np.max(mainriv_merg_info['BkfDepth'].values)

    tarinfo.loc[idx,'BasArea']       = np.sum(cbranch['BasArea'].values)
#    tarinfo.loc[idx,'NonLDArea']     = np.sum(cbranch['NonLDArea'].values)
    if len(cbranch) > 0:
        tarinfo.loc[idx,'BasSlope']      = np.average(cbranch['BasSlope'].values,  weights = cbranch['BasArea'].values)
        tarinfo.loc[idx,'MeanElev']      = np.average(cbranch['MeanElev'].values,  weights = cbranch['BasArea'].values)
        tarinfo.loc[idx,'BasAspect']     = np.average(cbranch['BasAspect'].values, weights = cbranch['BasArea'].values)

        tarinfo.loc[idx,'Max_DEM']       = np.max(cbranch['Max_DEM'].values)
        tarinfo.loc[idx,'Min_DEM']       = np.min(cbranch['Min_DEM'].values)
#    print(tarinfo.loc[idx,'BasArea'],"2")
    if Islake == 1:   ## Meger subbasin covered by lakes, Keep lake outlet catchment  DA, stream order info
        Longestpath = Calculate_Longest_flowpath(mainriv_merg_info)
        tarinfo.loc[idx,'RivLength'] = Longestpath

    elif Islake == 2:
        tarinfo.loc[idx,'RivLength'] = 0.0
        tarinfo.loc[idx,'IsLake']    = 2
    elif Islake < 0:
#        tarinfo.loc[idx,'Strahler']      = -1.2345
#        tarinfo.loc[idx,'Seg_ID']        = -1.2345
#        tarinfo.loc[idx,'Seg_order']     = -1.2345
#        tarinfo.loc[idx,'DA']            = -1.2345
        tarinfo.loc[idx,'HyLakeId']      = -1.2345
        tarinfo.loc[idx,'LakeVol']       = -1.2345
        tarinfo.loc[idx,'LakeArea']      = -1.2345
        tarinfo.loc[idx,'LakeDepth']     = -1.2345
        tarinfo.loc[idx,'Laketype']      = -1.2345
        tarinfo.loc[idx,'IsLake']        = -1.2345

    tarinfo.loc[idx,'centroid_x']    = -1.2345
    tarinfo.loc[idx,'centroid_y']    = -1.2345

    if seg_order >0 :
        tarinfo.loc[idx,'Seg_order']      = seg_order

    mask1 = mapoldnew_info['SubId'].isin(Modify_subids) ## catchment newly determined to be merged to target catchment
    mask2 = mapoldnew_info['nsubid'] == subid ###for subbasin already processed to drainage into this target catchment
    mask  = np.logical_or(mask1,mask2)


    ### the old downsub id of the dissolved polygon is stored in DowSubId
    for col in tarinfo.columns:
        if col == 'SubId':
#            print(tarinfo[col].values[0])
            mapoldnew_info.loc[mask,'nsubid']     = tarinfo[col].values[0]
#            print(mapoldnew_info.loc[mask,'nsubid'])
        elif col == 'nsubid' or col == 'ndownsubid' or col == 'Old_SubId' or col == 'Old_DowSubId':
            continue
        else:
            mapoldnew_info.loc[mask, col]         = tarinfo[col].values[0]
#    print(mapoldnew_info.loc[mask1,['BasArea','nsubid','SubId']])
#    print(mapoldnew_info.loc[mapoldnew_info['nsubid'] == subid,['BasArea','nsubid','SubId']])
    return mapoldnew_info

def Evaluate_Two_Dataframes(Expected,Result,Check_Col_NM = 'SubId'):
    ## check if two have the same column names 
    if (Expected.columns != Result.columns).all():
        return False 
    neql = 0
    ## check for each column two dataframe has the same value 
    for col in Expected.columns:
        Array_expect = Expected[col].values
        Array_result = Result[col].values
        if (Array_expect != Array_result).all():
            print(col)
            mask = Array_expect !=Array_result
            print(Expected[Check_Col_NM].values[mask])
        else:
            neql = neql + 1 
            
    return neql == len(Expected.columns) 
            


def UpdateTopology(mapoldnew_info,UpdateStreamorder = 1,UpdateSubId = 1):
    """ Functions will update subid,downsubid, calcuate stream order and 
        update drainage area in the attribute table mapoldnew_info
    ----------

    Notes
    -------

    Returns:
    -------
        mapoldnew_info
    """

    idx = mapoldnew_info.index

    if UpdateSubId > 0:
        for i in range(0,len(idx)):
            nsubid     = mapoldnew_info.loc[idx[i],'nsubid']
            subid      = mapoldnew_info.loc[idx[i],'SubId']
            odownsubid = mapoldnew_info.loc[idx[i],'DowSubId']

            donsubidinfo = mapoldnew_info.loc[mapoldnew_info['SubId'] == odownsubid].copy()

            if (len(donsubidinfo) >0):
                mapoldnew_info.loc[idx[i],'ndownsubid'] = donsubidinfo['nsubid'].values[0]
            else:
                mapoldnew_info.loc[idx[i],'ndownsubid'] = -1

        mapoldnew_info['Old_SubId']    = mapoldnew_info['SubId'].values
        mapoldnew_info['Old_DowSubId'] = mapoldnew_info['DowSubId'].values
        mapoldnew_info['SubId']        = mapoldnew_info['nsubid'].values

        mapoldnew_info['DowSubId'] = mapoldnew_info['ndownsubid'].values

    if UpdateStreamorder < 0:
        return mapoldnew_info

    mapoldnew_info_unique      = mapoldnew_info.drop_duplicates('SubId', keep='first')

    mapoldnew_info_unique      = Streamorderanddrainagearea(mapoldnew_info_unique)

    for i in range(0,len(mapoldnew_info_unique)):
        isubid    =  mapoldnew_info_unique['SubId'].values[i]
        mapoldnew_info.loc[mapoldnew_info['SubId'] == isubid,'Strahler']  = mapoldnew_info_unique['Strahler'].values[i]
        mapoldnew_info.loc[mapoldnew_info['SubId'] == isubid,'Seg_ID']    = mapoldnew_info_unique['Seg_ID'].values[i]
        mapoldnew_info.loc[mapoldnew_info['SubId'] == isubid,'Seg_order'] = mapoldnew_info_unique['Seg_order'].values[i]
        mapoldnew_info.loc[mapoldnew_info['SubId'] == isubid,'DA']        = mapoldnew_info_unique['DA'].values[i]

    return mapoldnew_info



def Streamorderanddrainagearea(catinfoall):
    """ Functions will  calcuate stream order and 
        update drainage area in the attribute table catinfoall
    ----------

    Notes
    -------

    Returns:
    -------
        catinfoall
    """
    catinfo                 = catinfoall.loc[catinfoall['IsLake'] != 2].copy()  ### remove none connected lake catchments, which do not connected to the river system
    catinfo_ncl             = catinfoall.loc[catinfoall['IsLake'] == 2].copy()
    routing_ncl             = catinfo_ncl[['SubId','DowSubId']].astype('float').values

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

            #### calculate DA of head watershed include None connected lakes
            if len(routing_ncl) == 0:
                 DA_ncl = 0.0
            else:
                Upstreamcats      = Defcat(routing_ncl,catid)     ### alll subuds
                Up_cat_info       = catinfo_ncl.loc[catinfo_ncl['SubId'].isin(Upstreamcats)].copy()

                if len(Up_cat_info) > 0:
                    DA_ncl            = sum(Up_cat_info['BasArea'].values)
                else:
                    DA_ncl            = 0.0

            catinfo.loc[idx,'DA'] = DA_ncl + catinfo['BasArea'].values[i]
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
            Up_Reaches_info = catinfo.loc[catinfo['DowSubId'] == catid].copy()
            cur_Reach_info = catinfo.loc[catinfo['SubId'] == catid].copy()
            curcat_idx = catinfo['SubId'] == catid

            #### calculate DA of None connected lakes
            if len(routing_ncl) == 0:
                DA_ncl = 0.0
            else:
#                print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
#                print(catid)
#                print(routing_ncl)
                Upstreamcats      = Defcat(routing_ncl,catid)     ### alll subuds
                Up_cat_info       = catinfo_ncl.loc[catinfo_ncl['SubId'].isin(Upstreamcats)].copy()
                if len(Up_cat_info) > 0:
                    DA_ncl            = sum(Up_cat_info['BasArea'].values)
                else:
                    DA_ncl            = 0.0


            if(len(cur_Reach_info) <= 0):  ### reach the most downstream of the watersheds
                break

            if len(Up_Reaches_info) == 1:   ### only have one upstream
                catinfo.loc[curcat_idx,'DA'] = cur_Reach_info['BasArea'].values[0] + Up_Reaches_info['DA'].values[0] + DA_ncl
                catinfo.loc[curcat_idx,'Strahler'] = Up_Reaches_info['Strahler'].values[0]
                catinfo.loc[curcat_idx,'Seg_order'] = Up_Reaches_info['Seg_order'].values[0] + 1
                catinfo.loc[curcat_idx,'Seg_ID'] = Up_Reaches_info['Seg_ID'].values[0]
#                print('1',catid,catinfo.loc[curcat_idx,'DA'].values,catinfo.loc[curcat_idx,'Strahler'].values,catinfo.loc[curcat_idx,'Sub_order'].values)
                catid =  int(cur_Reach_info['DowSubId'].values[0])
            else:  ### has mutiple upstram
                if np.min(Up_Reaches_info['Strahler'].values) > 0: ### all upstream has been processed
                    catinfo.loc[catinfo['SubId'] == catid,'DA'] = cur_Reach_info['BasArea'].values[0] + np.sum(Up_Reaches_info['DA'].values) + DA_ncl
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

    mask = catinfoall['SubId'].isin(catinfo['SubId'].values)
    catinfoall.loc[mask,'Seg_ID']    = catinfo['Seg_ID'].values
    catinfoall.loc[mask,'Seg_order'] = catinfo['Seg_order'].values
    catinfoall.loc[mask,'Strahler']  = catinfo['Strahler'].values
    catinfoall.loc[mask,'Seg_ID']    = catinfo['Seg_ID'].values
    catinfoall.loc[mask,'DA']        = catinfo['DA'].values

    ### calcuate channel manning's coefficient
    for i in range(0,len(catinfoall)):
        idx =  catinfoall.index[i]
        # if catinfo['BkfWidth'].values[i] > 0 and catinfo['RivSlope'].values[i] > 0 :
        #     catinfo.loc[idx,'Ch_n'] = calculateChannaln(catinfo['BkfWidth'].values[i],catinfo['BkfDepth'].values[i],
        #                       catinfo['Q_Mean'].values[i],catinfo['RivSlope'].values[i])
        if catinfoall['IsObs'].values[i] > 0:
            if catinfoall['DA_Obs'].values[i] >0:
                catinfoall.loc[idx,'DA_error'] = (catinfoall['DA'].values[i]/1000.0/1000.0 - catinfoall['DA_Obs'].values[i])/catinfoall['DA_Obs'].values[i]
    return catinfoall


def Defcat(out,outletid):
    """ Functions will return upstream ids in out taht drainage 
        to outletid
    ----------

    Notes
    -------

    Returns:
    -------
        Shedid
    """
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((len(out)+1,1),-99999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((len(out)+1,1),-99999999999)
        poshdid = 0
#        print("################################################a")
        for i in range(0,len(otsheds)):
#            print(otsheds)
#            print(psid,outletid)
            Shedid[psid] = otsheds[i]
#            print(Shedid[psid],otsheds[i])
#            print("##################################################b")
            psid = psid + 1
            irow = np.argwhere(rout[:,1]==otsheds[i]).astype(int)
#            print(len(irow))
            for j in range(0,len(irow)):
                #### if the catchment id already processed skip
                if rout[irow[j],0] in Shedid:
                    continue
                noutshd[poshdid] = rout[irow[j],0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd>=0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid>=0]
    return Shedid
    
            
                    