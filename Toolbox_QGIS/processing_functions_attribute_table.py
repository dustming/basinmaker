import numpy as np
import pandas as pd 

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
    