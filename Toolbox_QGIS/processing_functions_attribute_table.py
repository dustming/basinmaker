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