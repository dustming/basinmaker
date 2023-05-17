import copy

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from basinmaker.utilities.utilities import *
import numbers
from joblib import Parallel, delayed
import tempfile
pd.options.mode.chained_assignment = None


def update_selected_subid_using_sec_downsubid(sec_down_subinfo, upstream_subs, cat_ply, hyshdinfo):
    is_sec_down_subid_in_selected_subid = False
    is_subid_of_sec_downsubid_in_selected_subid = False
    update_downsubids_using_sec_downsubid = False
    update_topology = False
    # identify sec_down_subinfo subbains that downsubid is not in the subid list
    sec_down_subinfo_mostdown = sec_down_subinfo[~sec_down_subinfo['Sec_DowSubId'].isin(
        sec_down_subinfo['SubId'].values)]
    sec_down_subid_in_selected_subid = sec_down_subinfo_mostdown['Sec_DowSubId'].isin(
        upstream_subs).copy(deep=True)

    # check if downsubid in the secondary downsubid list exist in the selected subbasin
    # if not means the secondary downsubid has no impact on this area
    # return the input directly
    sec_down_subinfo_downsub_in_selected = sec_down_subinfo_mostdown[
        sec_down_subid_in_selected_subid]

    if len(sec_down_subinfo_downsub_in_selected) <= 0:
        return upstream_subs, cat_ply, update_topology

    # check subid of secondary downsub id that are in the selected subids
    # is in the selected subid or not
    # if they all inlucded in the selected subid
    # means the secondary downsubid has no impact on this area
    # return the input directly
    missing_subid_in_sec_table = sec_down_subinfo_downsub_in_selected[~sec_down_subinfo_downsub_in_selected['SubId'].isin(
        upstream_subs)].copy(deep=True)

    if len(missing_subid_in_sec_table) <= 0:
        return upstream_subs, cat_ply, update_topology

    # obtain routing topology for secondary table
    hyshdinfo_sec = sec_down_subinfo[[
        'SubId', 'Sec_DowSubId']].astype("int32").values
    for subid_sec in missing_subid_in_sec_table['SubId'].values:
        # obtain subids in the sec_down_subinfo drainage to this  subid_sec
        upstream_subs_sec = defcat(hyshdinfo_sec, subid_sec)
        sec_down_subinfoselect = sec_down_subinfo[sec_down_subinfo['SubId'].isin(
            upstream_subs_sec)]
        # update the DownSubId in the cat_ply using using Sec_DowSubId for upstream_subs
        cat_ply = cat_ply.merge(sec_down_subinfoselect, on='SubId', how='left')
        mask = cat_ply['SubId'].isin(upstream_subs_sec)
        cat_ply.loc[mask, 'DowSubId'] = cat_ply.loc[mask,
                                                    'Sec_DowSubId'].values
        # update the routing topology and get the subbasin drainage to this subid_sec
        hyshdinfo = cat_ply[['SubId', 'DowSubId']].astype("int32").values
        upstream_subs_new_sub_sec = defcat(hyshdinfo, subid_sec)
        # add the new subids into the eixisting upstream subids
        upstream_subs = np.concatenate(
            (upstream_subs, upstream_subs_new_sub_sec), axis=0)
        cat_ply = cat_ply.drop(columns='Sec_DowSubId')
    update_topology = True
    return upstream_subs, cat_ply, update_topology


def return_extracted_subids(cat_ply, mostdownid, mostupstreamid, sec_down_subinfo):

    # flags for sec down subid
    has_sec_downsub = False
    update_downsubids_using_sec_downsubid = False

    hyshdinfo = cat_ply[['SubId', 'DowSubId']].astype("int32").values
    sum_update_topology = 0
    update_topology = 0
    # find all subid control by this subid
    for i_down in range(0, len(mostdownid)):
        # Loop for each downstream id
        tar_subid = mostdownid[i_down]

        upstream_subs, cat_ply, update_topology = return_subids_drainage_to_subid(
            tar_subid, hyshdinfo, sec_down_subinfo, cat_ply)
        sum_update_topology = sum_update_topology + update_topology
        if i_down == 0:
            selected_subs = upstream_subs
        else:
            selected_subs = np.concatenate(
                (selected_subs, upstream_subs), axis=0)

    selected_subs = np.unique(selected_subs)

    # find all subid control by this subid
    remove_subs = np.empty(0, dtype=int)

    for i_up in range(0, len(mostupstreamid)):
        # Loop for each downstream id
        tar_subid = mostupstreamid[i_up]

        if tar_subid < 0:
            continue

        upstream_subs, cat_ply, update_topology = return_subids_drainage_to_subid(
            tar_subid, hyshdinfo, sec_down_subinfo, cat_ply)

        if i_up == 0:
            remove_subs = upstream_subs
        else:
            remove_subs = np.concatenate((remove_subs, upstream_subs), axis=0)

    if len(remove_subs) > 0:
        remove_subs = np.unique(remove_subs)

        mask = ~np.in1d(selected_subs, remove_subs)

        selected_subs = selected_subs[mask]
    if sum_update_topology >= 1:
        update_topology = True
    else:
        update_topology = False
    return selected_subs, cat_ply, update_topology


def return_subids_drainage_to_subid(tar_subid, hyshdinfo, sec_down_subinfo, cat_ply):

    # find all subid control by this subid
    upstream_subs = defcat(hyshdinfo, tar_subid)
    update_topology = 0
    # check if has sencondary down subid
    if len(sec_down_subinfo) > 0:
        upstream_subs, cat_ply, update_topology = update_selected_subid_using_sec_downsubid(
            sec_down_subinfo, upstream_subs, cat_ply, hyshdinfo)

    return upstream_subs, cat_ply, update_topology


def remove_landuse_type_input_based_on_area(landuse_thres, hruinfo, sub_area, Landuse_ID):

    if landuse_thres <= 0:
        return hruinfo
    # calculate the landuse area of each landuse group in each subbasin
    subinfo_lu = hruinfo[['SubId', Landuse_ID, 'HRU_Area']].copy(deep=True)
    subinfo_lu = subinfo_lu.rename(columns={"HRU_Area": "Input_A_G"})
    subinfo_lu = subinfo_lu.groupby(
        ['SubId', Landuse_ID], as_index=False).sum()

    # calcuate landuse area ratio
    subinfo_lu = pd.merge(subinfo_lu, sub_area, on='SubId')
    subinfo_lu['Area_ratio'] = subinfo_lu["Input_A_G"]/subinfo_lu['Bas_A_G']

    # obtain dominated landuse ID
    subinfo_lu = subinfo_lu.sort_values(
        by=['SubId', 'Input_A_G'], ascending=False)
    subinfo_lu_dominated = subinfo_lu.drop_duplicates(
        subset=['SubId'], keep='first').copy(deep=True)

    # find sub
    subinfo_lu_need_change = subinfo_lu[subinfo_lu['Area_ratio'] < landuse_thres][[
        'SubId', Landuse_ID]].copy(deep=True)

    for i in range(0, len(subinfo_lu_need_change)):
        subid = subinfo_lu_need_change['SubId'].values[i]
        landuse = subinfo_lu_need_change[Landuse_ID].values[i]
        mask1 = hruinfo['SubId'] == subid
        mask2 = hruinfo[Landuse_ID] == landuse
        mask = np.logical_and(mask1, mask2)
        hruinfo.loc[mask, Landuse_ID] = subinfo_lu_dominated.loc[subinfo_lu_dominated['SubId']
                                                                 == subid][Landuse_ID].values[0]

    return hruinfo


def simplify_hrus_method2(area_ratio_thresholds, hruinfo, Landuse_ID,
                          Soil_ID, Veg_ID, Other_Ply_ID_1, Other_Ply_ID_2):

    hru_lake_info = hruinfo.loc[hruinfo['HRU_IsLake'] > 0].copy()
    hru_land_info = hruinfo.loc[hruinfo['HRU_IsLake'] <= 0].copy()

    sub_area = hruinfo[['SubId', 'HRU_Area']].copy(deep=True)
    sub_area = sub_area.rename(columns={"HRU_Area": "Bas_A_G"})
    sub_area = sub_area.groupby(['SubId'], as_index=False).sum()
#    hruinfo = pd.merge(hruinfo, subinfo, on='SubId')

    landuse_thres = area_ratio_thresholds[0]
    list = [Landuse_ID, Soil_ID, Other_Ply_ID_1]
    for i in range(0, len(list)):
        Item = list[i]
        landuse_thres = area_ratio_thresholds[i]
        hru_land_info = remove_landuse_type_input_based_on_area(
            landuse_thres, hru_land_info, sub_area, Item)
    hru_land_info[Veg_ID] = hru_land_info[Landuse_ID]

    hruinfo = hru_lake_info.append(hru_land_info)

    return hruinfo


def simplidfy_hrus(min_hru_pct_sub_area, hruinfo, importance_order):

    hruinfo['HRU_ID_New2'] = hruinfo['HRU_ID_New']

    subids = np.unique(hruinfo['SubId'].values)

    # loop for each subbasin
    for i in range(0, len(subids)):
        subid = subids[i]
        # hrus in this subbasin
        sub_hru_info = hruinfo.loc[hruinfo['SubId'] == subid].copy(deep=True)
        # calculate area and subbasin minimum hur ares
        subasrea = np.sum(sub_hru_info['HRU_Area'].values)
        subarea_thrs = min_hru_pct_sub_area * subasrea

        # obtain hru need to be removed
        need_remove_hrus = sub_hru_info.loc[sub_hru_info['HRU_Area'] < subarea_thrs].copy(
        )
        # obtain hrus that will be keepted
        good_hrus = sub_hru_info.loc[sub_hru_info['HRU_Area']
                                     >= subarea_thrs].copy()

        hru_columns = good_hrus.columns
        # do not modify this columns
        hru_columns = hru_columns[hru_columns != 'HRU_ID_New2']
        # check if the need_remove_hru can merge togeher for importance order
        # 1
        colnm1 = importance_order[0]
        import1_Area_need_remove_hru = need_remove_hrus.copy(deep=True)
        unique_import1 = np.unique(import1_Area_need_remove_hru[colnm1].values)
        # if subid == 116:
        #     print(subarea_thrs,min_hru_pct_sub_area,subasrea)
        #     print(need_remove_hrus[['HRU_ID_New2','HRU_Area',colnm1]])

        for i_import in range(0, len(unique_import1)):
            i_colnm1 = unique_import1[i_import]
            i_import1_Area_need_remove_hru = import1_Area_need_remove_hru.loc[import1_Area_need_remove_hru[colnm1] == i_colnm1].copy(
                deep=True)
            total_area_i_import_in_need_remove_hrus = np.sum(
                i_import1_Area_need_remove_hru['HRU_Area'].values)
            # check total area of the most important column in the importance list
            # if it is larger than area threthold,
            # then merge them together

            if total_area_i_import_in_need_remove_hrus >= subarea_thrs:

                # merge to the hru with largest area in the list
                i_import1_Area_need_remove_hru = i_import1_Area_need_remove_hru.sort_values(
                    by=['HRU_Area'], ascending=False)
                for i_nm in range(0, len(hru_columns)):
                    columnname = hru_columns[i_nm]
                    if columnname == 'SHAPE':
                        continue
                    # modify the hru attributes to good hru attribute
                    hruinfo.loc[hruinfo['HRU_ID_New2'].isin(
                        i_import1_Area_need_remove_hru['HRU_ID_New2'].values), columnname] = i_import1_Area_need_remove_hru[columnname].values[0]

                # remove modified HRU from need removed hru list
                need_remove_hrus = need_remove_hrus[~need_remove_hrus['HRU_ID_New2'].isin(
                    i_import1_Area_need_remove_hru['HRU_ID_New2'].values)]
        # loop for each
        indexes = need_remove_hrus.index
        # if subid == 116:
        #     print(need_remove_hrus[['HRU_ID_New2','HRU_Area',colnm1]])

        for j in range(0, len(indexes)):
            idx = indexes[j]
            hruid = need_remove_hrus.loc[idx, 'HRU_ID_New2']
            # find a target hru from good_hrus, and merge them by change attribute

            # loop the importance_order,
            for k in range(0, len(importance_order)):
                colnm = importance_order[k]
                # find the attribute value of the problem hru for this column
                attri_remove_hru = need_remove_hrus.loc[idx, colnm]

                # check if there is good hrus has the same attribute value
                good_hrus_k = good_hrus.loc[good_hrus[colnm]
                                            == attri_remove_hru].copy()

                if len(good_hrus_k) > 0:
                    # sort by hru areas
                    good_hrus_k = good_hrus_k.sort_values(
                        by='HRU_Area', ascending=False)
                    for i_nm in range(0, len(hru_columns)):
                        columnname = hru_columns[i_nm]
                        if columnname == 'SHAPE':
                            continue
                        # modify the hru attributes to good hru attribute
                        hruinfo.loc[hruinfo['HRU_ID_New2'] == hruid,
                                    columnname] = good_hrus_k[columnname].values[0]
                else:
                    continue
    hruinfo = hruinfo.drop(columns=['HRU_ID_New2'])
    return hruinfo


def update_non_connected_catchment_info(catinfo):
    routing_info = catinfo[["SubId", "DowSubId"]].astype("float").values
    catinfo_non_connected = catinfo.loc[catinfo["Lake_Cat"] == 2].copy()

    catids_nc = catinfo_non_connected["SubId"].copy()
    max_seg_id = catinfo["Seg_ID"].max()
    # catinfo.loc[
    #     catinfo["SubId"].isin(catids_nc), "RivLength"
    # ] = 0.0  ## no reiver length since not connected.

    for i in range(0, len(catinfo_non_connected)):
        c_subid = catinfo_non_connected["SubId"].values[i]
        d_subid = catinfo_non_connected["DowSubId"].values[i]
        d_sub_info = catinfo.loc[catinfo["SubId"] == d_subid].copy()

        lc_subid = d_subid

        # Upstreamcats = defcat(routing_info, c_subid)  ### alll subuds
        #
        # Up_cat_info = catinfo.loc[catinfo["SubId"].isin(Upstreamcats)].copy()
        #
        # DA = sum(Up_cat_info["BasArea"].values)
        #
        # catinfo.loc[catinfo["SubId"] == c_subid, "DrainArea"] = DA

        if len(d_sub_info) < 1:
            catinfo.loc[catinfo["SubId"] == c_subid, "Strahler"] = 1
            catinfo.loc[catinfo["SubId"] == c_subid,
                        "Seg_ID"] = max_seg_id + i + 1
            catinfo.loc[catinfo["SubId"] == c_subid, "Seg_order"] = 1
            continue

        # add nonconnected lake catchment area to downsubbasin drinage area
        #        if d_sub_info['Lake_Cat'].values[0]  != 2:
        #            catinfo.loc[catinfo['SubId'] == d_subid,'DA'] = d_sub_info['DA'].values[0] + DA

        while d_sub_info["Lake_Cat"].values[0] == 2:
            lc_subid_info = catinfo.loc[catinfo["SubId"] == lc_subid].copy()
            d_subid = lc_subid_info["DowSubId"].values[0]
            d_sub_info = catinfo.loc[catinfo["SubId"] == d_subid].copy()
            if len(d_sub_info) < 1:
                lc_subid = -1
                break
            if lc_subid == d_sub_info['DowSubId'].values[0]:
                print(lc_subid, d_subid)
                lc_subid = -1
                break
            lc_subid = d_subid

        if lc_subid == -1:
            continue

        # catinfo.loc[catinfo["SubId"] == c_subid, "RivSlope"] = -1.2345
        # catinfo.loc[catinfo["SubId"] == c_subid, "Ch_n"] = -1.2345
        # catinfo.loc[catinfo["SubId"] == c_subid, "RivLength"] = -1.2345
        # catinfo.loc[catinfo["SubId"] == c_subid, "Min_DEM"] = -1.2345
        # catinfo.loc[catinfo["SubId"] == c_subid, "Max_DEM"] = -1.2345
        # catinfo.loc[catinfo["SubId"] == c_subid, "FloodP_n"] = -1.2345
        # if "DA_Chn_L" in catinfo.columns:
        #     catinfo.loc[catinfo["SubId"] == c_subid, "DA_Chn_L"] = -1.2345
        #     catinfo.loc[catinfo["SubId"] == c_subid, "DA_Chn_Slp"] = -1.2345

        # catinfo.loc[catinfo["SubId"] == c_subid, "Q_Mean"] = d_sub_info[
        #     "Q_Mean"
        # ].values[0]
        # catinfo.loc[catinfo["SubId"] == c_subid, "BkfWidth"] = d_sub_info[
        #     "BkfWidth"
        # ].values[0]
        # catinfo.loc[catinfo["SubId"] == c_subid, "BkfDepth"] = d_sub_info[
        #     "BkfDepth"
        # ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Strahler"] = d_sub_info[
            "Strahler"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Seg_ID"] = d_sub_info[
            "Seg_ID"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Seg_order"] = d_sub_info[
            "Seg_order"
        ].values[0]

    return catinfo


def Calculate_Longest_flowpath(mainriv_merg_info):
    mainriv_merg_info_sort = mainriv_merg_info.sort_values(
        ["DrainArea"], ascending=(False))
    #    print(mainriv_merg_info_sort[['SubId','DowSubId','DA','Strahler','RivLength']])
    longest_flow_pathes = np.full(100, 0)
    #    print(longest_flow_pathes)
    npath = 1

    # loop upstream to find longest flow path
    Pathid = np.full(1000, -1)
    subid = mainriv_merg_info_sort["SubId"].values[0]
    npath_current = 1
    Pathid[npath - 1] = subid
    #    print('####################################################################################')
    while len(Pathid[Pathid > 0]) > 0:
        nPathid = np.full(1000, -1)
        npath = npath_current

        #        print('###################################')
        #        print(npath,Pathid[0:npath])
        #        print('###################################')
        for ipath in range(0, npath_current):
            c_subid_ipath = Pathid[ipath]

            if (
                c_subid_ipath < 0
            ):  # means this path has been closed due to no more subbasin within the lake domain
                continue

            longest_flow_pathes[ipath] = (
                mainriv_merg_info_sort.loc[
                    mainriv_merg_info_sort["SubId"] == c_subid_ipath, "RivLength"
                ]
                + longest_flow_pathes[ipath]
            )  # add river length to current path
            Strahler_order_ipath = mainriv_merg_info_sort.loc[
                mainriv_merg_info_sort["SubId"] == c_subid_ipath, "Strahler"
            ].values[0]

            upstream_sub_infos = mainriv_merg_info_sort.loc[
                mainriv_merg_info_sort["DowSubId"] == c_subid_ipath
            ]  # get upstream info

            if (
                len(upstream_sub_infos) <= 0
            ):  # no more upstream catchment within the domain of the lake
                #                print("path        closed        ",ipath)
                continue

            # look for upstream catchment has the same upstream_sub_infos_eq_Strahler first
            #            print(Strahler_order_ipath)
            #            print(upstream_sub_infos['Strahler'])
            upstream_sub_infos_eq_Strahler = upstream_sub_infos.loc[
                upstream_sub_infos["Strahler"] == Strahler_order_ipath
            ]

            if (
                len(upstream_sub_infos_eq_Strahler) > 0
            ):  # has a upstream river has the saem strahler id, no new path will be added
                nPathid[ipath] = upstream_sub_infos_eq_Strahler["SubId"].values[
                    0
                ]  # add this upstream id to nPathid
                continue
            else:
                upstream_sub_infos_eq_Strahler_1 = upstream_sub_infos.loc[
                    upstream_sub_infos["Strahler"] == Strahler_order_ipath - 1
                ]

                for inpath in range(0, len(upstream_sub_infos_eq_Strahler_1)):
                    # this brance sperate into two or several reaches, the starting river length for all of them are the same
                    if inpath == 0:
                        nPathid[ipath] = upstream_sub_infos_eq_Strahler_1[
                            "SubId"
                        ].values[inpath]
                    #                        print(nPathid[ipath],ipath,upstream_sub_infos_eq_Strahler_1['SubId'].values[inpath],'aaaaa',range(0,len(upstream_sub_infos_eq_Strahler_1)))
                    else:
                        nPathid[npath + 1 - 1] = upstream_sub_infos_eq_Strahler_1[
                            "SubId"
                        ].values[inpath]
                        longest_flow_pathes[npath + 1 -
                                            1] = longest_flow_pathes[ipath]
                        #                        print(npath + 1 - 1,longest_flow_pathes[npath + 1 - 1],nPathid[npath + 1 - 1],'bbbbb',range(0,len(upstream_sub_infos_eq_Strahler_1)))
                        npath = npath + 1

        Pathid = nPathid
        npath_current = npath
    Longestpath = max(longest_flow_pathes)

    return Longestpath


def remove_possible_small_subbasins(mapoldnew_info, area_thresthold=50, length_thresthold=15):
    mapoldnew_info_new = mapoldnew_info.copy(deep=True)
    # check the gauge column name, in case it is v2.1 using Has_Gauge
    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in mapoldnew_info.columns:
        Gauge_col_Name = "Has_Gauge"

    # get small subbasin that is not lake
    area_filter = mapoldnew_info['BasArea'] / 1000/1000 < area_thresthold
    length_filter = mapoldnew_info['RivLength'] / 1000 < length_thresthold
    small_sub_non_lake = mapoldnew_info[np.logical_or(
        area_filter, length_filter)].copy(deep=True)
    small_sub_non_lake = small_sub_non_lake[small_sub_non_lake['Lake_Cat'] == 0].copy(
        deep=True)
    small_sub_non_lake = small_sub_non_lake[small_sub_non_lake[Gauge_col_Name] == 0].copy(
        deep=True)
    small_sub_non_lake = small_sub_non_lake.sort_values(['Strahler','DrainArea'], ascending=True)
    # process connected lakes  merge polygons
    for i in range(0, len(small_sub_non_lake)):
        small_sub_id = small_sub_non_lake['SubId'].values[i]
        small_downsub_id = small_sub_non_lake['DowSubId'].values[i]
        small_sub_seg_id = small_sub_non_lake['Seg_ID'].values[i]

        small_sub_is_not_Lake = small_sub_non_lake['Lake_Cat'].values[i] == 0
        small_sub_is_not_gauge = small_sub_non_lake[Gauge_col_Name].values[i] <= 0

        down_sub_info = mapoldnew_info[mapoldnew_info['SubId'] == small_downsub_id].copy(
            deep=True)
        upstream_sub_info = mapoldnew_info[mapoldnew_info['DowSubId'] == small_sub_id].copy(
            deep=True)

        # check if it has the same segment id with downstream subbasin
        if len(down_sub_info) > 0:
            has_down_sub = True
            if down_sub_info['Seg_ID'].values[0] == small_sub_seg_id:
                down_sub_has_same_seg_id = True
            else:
                down_sub_has_same_seg_id = False
            if down_sub_info['Lake_Cat'].values[0] > 0:
                down_sub_is_lake = True
            else:
                down_sub_is_lake = False
        else:
            has_down_sub = False
            down_sub_has_same_seg_id = False
            down_sub_is_lake = False
        case = -1

        # if down stream sub is a lake sub, the seg_id was changed to the lake outlet subid
        if has_down_sub and down_sub_has_same_seg_id and small_sub_is_not_Lake and small_sub_is_not_gauge:
            #            tarinfo = down_sub_info
            # merge to downstream subbasin id
            tarinfo = mapoldnew_info_new[mapoldnew_info_new['SubId']
                                         == down_sub_info['SubId'].values[0]].copy(deep=True)
            modify = True
            down_sub_info_2 = mapoldnew_info_new[mapoldnew_info_new['SubId']
                                                 == down_sub_info['SubId'].values[0]].copy(deep=True)
            ndown_subid = down_sub_info_2['DowSubId'].values[0]
            update_river = True
            case = 1

        elif down_sub_is_lake and has_down_sub and case == -1:
            # merge to downstream subbasin id
            #            tarinfo = down_sub_info
            tarinfo = mapoldnew_info_new[mapoldnew_info_new['SubId']
                                         == down_sub_info['SubId'].values[0]].copy(deep=True)
            modify = True
            down_sub_info_2 = mapoldnew_info_new[mapoldnew_info_new['SubId']
                                                 == down_sub_info['SubId'].values[0]].copy(deep=True)
            ndown_subid = down_sub_info_2['DowSubId'].values[0]
            update_river = False

        elif has_down_sub and case == -1:
            # merge to downstream subbasin id
            #            tarinfo = down_sub_info
            tarinfo = mapoldnew_info_new[mapoldnew_info_new['SubId']
                                         == down_sub_info['SubId'].values[0]].copy(deep=True)
            down_sub_info_2 = mapoldnew_info_new[mapoldnew_info_new['SubId']
                                                 == down_sub_info['SubId'].values[0]].copy(deep=True)
            ndown_subid = down_sub_info_2['DowSubId'].values[0]

            target_upstream_info = mapoldnew_info[mapoldnew_info['DowSubId'] == small_downsub_id].copy(
                deep=True)
            upstream_subid_tar = np.unique(
                target_upstream_info['SubId'].values)

            common_elements = np.intersect1d(
                upstream_subid_tar, small_sub_non_lake['SubId'].values)

            if len(common_elements) == 1:
                update_river = True
                modify = True
            else:
                upstream_subid_need_modify = target_upstream_info[target_upstream_info['SubId'].isin(
                    common_elements)].sort_values('DrainArea', ascending=False)
                if small_sub_id != upstream_subid_need_modify['SubId'].values[0]:
                    update_river = False
                else:
                    update_river = True
                modify = True

        else:
            modify = False
            if small_downsub_id > 0:
                print(
                    small_sub_id,
                    has_down_sub,
                    down_sub_has_same_seg_id,
                    small_sub_is_not_Lake,
                    small_sub_is_not_gauge,
                    has_down_sub,
                    len(mapoldnew_info[mapoldnew_info['Seg_ID']
                        == small_sub_seg_id])
                )
            tarinfo = []
            continue

        if modify:
            mask1 = mapoldnew_info_new['SubId'] == small_sub_id
            mask3 = mapoldnew_info_new['nsubid'] == small_sub_id
            # for subbasin already processed to drainage into this target catchment
            mask2 = (mapoldnew_info_new["nsubid"]
                     == tarinfo["nsubid"].values[0])
            masktemp = np.logical_or(mask1, mask2)
            mask = np.logical_or(masktemp, mask3)

            if len(tarinfo) > 1 and len(np.unique(tarinfo['HyLakeId'].values)) > 1:
                print(small_sub_id, small_sub_seg_id)
                print(tarinfo[['SubId', 'DowSubId', 'Seg_ID', 'HyLakeId']])

            attributes = mapoldnew_info_new[mapoldnew_info_new['SubId'].isin(
                [small_sub_id, down_sub_info['SubId'].values[0]])]
            mapoldnew_info_new.loc[mask, 'BasArea'] = np.sum(
                attributes["BasArea"].values)
            mapoldnew_info_new.loc[mask, "BasSlope"] = np.average(
                attributes["BasSlope"].values, weights=attributes["BasArea"].values
            )
            mapoldnew_info_new.loc[mask, "MeanElev"] = np.average(
                attributes["MeanElev"].values, weights=attributes["BasArea"].values
            )
            mapoldnew_info_new.loc[mask, "BasAspect"] = np.average(
                attributes["BasAspect"].values, weights=attributes["BasArea"].values
            )
            for col in tarinfo.columns:
                if col == "RivLength" and update_river:
                    mapoldnew_info_new.loc[mask, "RivLength"] = np.sum(
                        attributes["RivLength"].values)
                elif col in ["RivSlope", "FloodP_n", "Q_Mean", "Ch_n"] and update_river:
                    mapoldnew_info_new.loc[mask, col] = np.average(
                        attributes[col].values,
                        weights=attributes[col].values,
                    )
                elif col in ["Max_DEM"] and update_river:
                    mapoldnew_info_new.loc[mask, col] = np.max(
                        attributes[col].values
                    )
                elif col in ["Min_DEM"] and update_river:
                    mapoldnew_info_new.loc[mask, col] = np.min(
                        attributes[col].values
                    )

                elif col == 'SubId':

                    mapoldnew_info_new.loc[mask,
                                           "nsubid"] = tarinfo["nsubid"].values[0]

                elif col == 'DowSubId':
                    mapoldnew_info_new.loc[mask, col] = ndown_subid
                elif (
                    col == "nsubid"
                    or col == "ndownsubid"
                    or col == "Old_SubId"
                    or col == "Old_DowSubId"
                    or col == "SHAPE"
                    or col == "geometry"
                    or col == 'BasArea'
                    or col == 'BasSlope'
                    or col == 'MeanElev'
                    or col == 'BasAspect'
                ):
                    continue
                else:
                    mapoldnew_info_new.loc[mask, col] = tarinfo[col].values[0]

    return mapoldnew_info_new


def New_SubId_To_Dissolve(
    subid,
    catchmentinfo,
    mapoldnew_info,
    upsubid=-1,
    ismodifids=-1,
    modifiidin=[-1],
    mainriv=[-1],
    Lake_Cat=-1,
    seg_order=-1,
):
    sub_colnm = "SubId"
    routing_info = catchmentinfo[["SubId", "DowSubId"]].astype("int32").values
    if ismodifids < 0:
        Modify_subids1 = defcat(
            routing_info, subid
        )  # find all subids drainage to this subid
        if upsubid > 0:
            Modify_subids2 = defcat(routing_info, upsubid)
            mask = np.in1d(Modify_subids1, Modify_subids2)
            Modify_subids = Modify_subids1[np.logical_not(mask)]
        else:
            Modify_subids = Modify_subids1

    else:
        Modify_subids = modifiidin
    #    print("##########################################")
    #    print(subid)
    #    print(Modify_subids)
    cbranch = catchmentinfo[catchmentinfo[sub_colnm].isin(
        Modify_subids)].copy()
    tarinfo = catchmentinfo[
        catchmentinfo[sub_colnm] == subid
    ].copy()  # define these subs attributes
    # average river slope info
    mainriv_merg_info = mainriv.loc[mainriv["SubId"].isin(
        Modify_subids)].copy()
    mainriv_merg_info = mainriv_merg_info.loc[mainriv_merg_info["RivLength"] > 0].copy(
    )
    idx = tarinfo.index[0]
    #    print(tarinfo.loc[idx,'BasArea'],"1")
    if len(mainriv_merg_info) > 0:
        tarinfo.loc[idx, "RivLength"] = np.sum(
            mainriv_merg_info["RivLength"].values)
        tarinfo.loc[idx, "RivSlope"] = np.average(
            mainriv_merg_info["RivSlope"].values,
            weights=mainriv_merg_info["RivLength"].values,
        )
        tarinfo.loc[idx, "FloodP_n"] = np.average(
            mainriv_merg_info["FloodP_n"].values,
            weights=mainriv_merg_info["RivLength"].values,
        )
        tarinfo.loc[idx, "Q_Mean"] = np.average(
            mainriv_merg_info["Q_Mean"].values,
            weights=mainriv_merg_info["RivLength"].values,
        )
        tarinfo.loc[idx, "Ch_n"] = np.average(
            mainriv_merg_info["Ch_n"].values,
            weights=mainriv_merg_info["RivLength"].values,
        )
        tarinfo.loc[idx, "BkfWidth"] = np.max(
            mainriv_merg_info["BkfWidth"].values)
        tarinfo.loc[idx, "BkfDepth"] = np.max(
            mainriv_merg_info["BkfDepth"].values)

        tarinfo.loc[idx, "Max_DEM"] = np.max(
            mainriv_merg_info["Max_DEM"].values)
        tarinfo.loc[idx, "Min_DEM"] = np.min(
            mainriv_merg_info["Min_DEM"].values)

    tarinfo.loc[idx, "BasArea"] = np.sum(cbranch["BasArea"].values)
    #    tarinfo.loc[idx,'NonLDArea']     = np.sum(cbranch['NonLDArea'].values)
    if len(cbranch) > 0:
        tarinfo.loc[idx, "BasSlope"] = np.average(
            cbranch["BasSlope"].values, weights=cbranch["BasArea"].values
        )
        tarinfo.loc[idx, "MeanElev"] = np.average(
            cbranch["MeanElev"].values, weights=cbranch["BasArea"].values
        )
        tarinfo.loc[idx, "BasAspect"] = np.average(
            cbranch["BasAspect"].values, weights=cbranch["BasArea"].values
        )
    #    print(tarinfo.loc[idx,'BasArea'],"2")
    if (
        Lake_Cat == 1
    ):  # Meger subbasin covered by lakes, Keep lake outlet catchment  DA, stream order info
        #        Longestpath = Calculate_Longest_flowpath(mainriv_merg_info)
        # tarinfo.loc[idx, "RivSlope"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "RivLength"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "FloodP_n"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "Ch_n"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "Min_DEM"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "Max_DEM"] = -1.2345  # Longestpath
        # if "DA_Chn_L" in tarinfo.columns:
        #     tarinfo.loc[idx, "DA_Chn_L"] = -1.2345  # Longestpath
        #     tarinfo.loc[idx, "DA_Chn_Slp"] = -1.2345  # Longestpath
        do_nothing = 1

    elif Lake_Cat == 2:
        # tarinfo.loc[idx, "RivSlope"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "RivLength"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "FloodP_n"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "Ch_n"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "Min_DEM"] = -1.2345  # Longestpath
        # tarinfo.loc[idx, "Max_DEM"] = -1.2345  # Longestpath
        # if "DA_Chn_L" in tarinfo.columns:
        #     tarinfo.loc[idx, "DA_Chn_L"] = -1.2345  # Longestpath
        #     tarinfo.loc[idx, "DA_Chn_Slp"] = -1.2345  # Longestpath
        tarinfo.loc[idx, "Lake_Cat"] = 2
    elif Lake_Cat <= 0:
        #        tarinfo.loc[idx,'Strahler']      = -1.2345
        #        tarinfo.loc[idx,'Seg_ID']        = -1.2345
        #        tarinfo.loc[idx,'Seg_order']     = -1.2345
        #        tarinfo.loc[idx,'DA']            = -1.2345
        tarinfo.loc[idx, "HyLakeId"] = 0
        tarinfo.loc[idx, "LakeVol"] = 0
        tarinfo.loc[idx, "LakeArea"] = 0
        tarinfo.loc[idx, "LakeDepth"] = 0
        tarinfo.loc[idx, "Laketype"] = 0
        tarinfo.loc[idx, "Lake_Cat"] = 0

    tarinfo.loc[idx, "centroid_x"] = -1.2345
    tarinfo.loc[idx, "centroid_y"] = -1.2345

    if seg_order > 0:
        tarinfo.loc[idx, "Seg_order"] = seg_order

    mask1 = mapoldnew_info["SubId"].isin(
        Modify_subids
    )  # catchment newly determined to be merged to target catchment
    mask2 = (
        mapoldnew_info["nsubid"] == subid
    )  # for subbasin already processed to drainage into this target catchment
    mask = np.logical_or(mask1, mask2)

    # the old downsub id of the dissolved polygon is stored in DowSubId
    for col in tarinfo.columns:
        if col == "SubId":
            #            print(tarinfo[col].values[0])
            mapoldnew_info.loc[mask, "nsubid"] = tarinfo[col].values[0]
        #            print(mapoldnew_info.loc[mask,'nsubid'])
        elif (
            col == "nsubid"
            or col == "ndownsubid"
            or col == "Old_SubId"
            or col == "Old_DowSubId"
            or col == "SHAPE"
            or col == "geometry"
        ):
            continue
        else:
            mapoldnew_info.loc[mask, col] = tarinfo[col].values[0]
    #    print(mapoldnew_info.loc[mask1,['BasArea','nsubid','SubId']])
    #    print(mapoldnew_info.loc[mapoldnew_info['nsubid'] == subid,['BasArea','nsubid','SubId']])
    return mapoldnew_info


def Evaluate_Two_Dataframes(Expected, Result, Check_Col_NM="SubId"):
    # check if two have the same column names
    if (Expected.columns != Result.columns).all():
        print(Expected.columns)
        print(Result.columns)
        return False
    neql = 0
    # check for each column two dataframe has the same value
    for col in Expected.columns:
        if col == "HRU_ID":
            neql = neql + 1
            continue
        Array_expect = Expected[col].values
        Array_result = Result[col].values
        if (Array_expect != Array_result).all():
            print(col)
            mask = Array_expect != Array_result
            print(Expected[Check_Col_NM].values[mask])
        else:
            neql = neql + 1

    return neql == len(Expected.columns)


def UpdateTopology(mapoldnew_info, UpdateStreamorder=1, UpdateSubId=1):
    """Functions will update subid,downsubid, calcuate stream order and
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
        for i in range(0, len(idx)):
            nsubid = mapoldnew_info.loc[idx[i], "nsubid"]
            subid = mapoldnew_info.loc[idx[i], "SubId"]
            odownsubid = mapoldnew_info.loc[idx[i], "DowSubId"]

            donsubidinfo = mapoldnew_info.loc[
                mapoldnew_info["SubId"] == odownsubid
            ].copy()

            if len(donsubidinfo) > 0:
                mapoldnew_info.loc[idx[i], "ndownsubid"] = donsubidinfo[
                    "nsubid"
                ].values[0]
            else:
                mapoldnew_info.loc[idx[i], "ndownsubid"] = -1

        mapoldnew_info["Old_SubId"] = mapoldnew_info["SubId"].values
        mapoldnew_info["Old_DowSubId"] = mapoldnew_info["DowSubId"].values
        mapoldnew_info["SubId"] = mapoldnew_info["nsubid"].values

        mapoldnew_info["DowSubId"] = mapoldnew_info["ndownsubid"].values

    if UpdateStreamorder < 0:
        return mapoldnew_info

    mapoldnew_info_unique = mapoldnew_info.drop_duplicates(
        "SubId", keep="first")

    mapoldnew_info_unique = streamorderanddrainagearea(mapoldnew_info_unique)

    for i in range(0, len(mapoldnew_info_unique)):
        isubid = mapoldnew_info_unique["SubId"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "Strahler"
        ] = mapoldnew_info_unique["Strahler"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "Seg_ID"
        ] = mapoldnew_info_unique["Seg_ID"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "Seg_order"
        ] = mapoldnew_info_unique["Seg_order"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "DrainArea"
        ] = mapoldnew_info_unique["DrainArea"].values[i]

    return mapoldnew_info


def calculate_Tc(DrainArea, DA_Chn_L, DA_Chn_Slp):

    TC_1 = 0.675*(DrainArea/1000/1000)**0.5
    TC_2 = 0.2426 * (DA_Chn_L/1000)*((DrainArea/1000/1000)
                                     ** (-0.1)) * (DA_Chn_Slp**(-0.2))
    TC_3 = 0.3 * ((DA_Chn_L/1000)**0.76) * (DA_Chn_Slp**(-0.19))
    TC_4 = (1/0.6)*2.8*((DA_Chn_L/1000) / (DA_Chn_Slp**(0.5)))**0.47
    TC_5 = (1/0.6)*0.000326*((DA_Chn_L) / (DA_Chn_Slp**(0.5)))**0.79
    return TC_1, TC_2, TC_3, TC_4, TC_5


def streamorderanddrainagearea(catinfoall):
    """Functions will  calcuate stream order and
        update drainage area in the attribute table catinfoall
    ----------

    Notes
    -------

    Returns:
    -------
        catinfoall
    """
    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in catinfoall.columns:
        Gauge_col_Name = "Has_Gauge"

    # remove none connected lake catchments, which do not connected to the river system
    catinfo = catinfoall.loc[catinfoall["Lake_Cat"] != 2].copy()
    catinfo_ncl = catinfoall.loc[catinfoall["Lake_Cat"] == 2].copy()
    routing_ncl = catinfo_ncl[["SubId", "DowSubId"]].astype("float").values

    catlist = np.full((len(catinfo)), -9)
    icat = 0
    iseg = 1
    # find first segments of all reaches, no upstream reaches
    for i in range(0, len(catinfo)):
        idx = catinfo.index[i]
        if catinfo["SubId"].values[i] == catinfo["DowSubId"].values[i]:
            catinfo.loc[idx, "DowSubId"] = -1
        catid = catinfo["SubId"].values[i]
        # the river seg has no upstream segment
        if (len(catinfo[catinfo["DowSubId"] == catid]) == 0):
            # store next reach segment
            catlist[icat] = int(catinfo["DowSubId"].values[i])

            # calculate DA of head watershed include None connected lakes
            if len(routing_ncl) == 0:
                DA_ncl = 0.0
                slp_ncl = 0.0
            else:
                Upstreamcats = defcat(routing_ncl, catid)  # alll subuds
                Up_cat_info = catinfo_ncl.loc[
                    catinfo_ncl["SubId"].isin(Upstreamcats)
                ].copy()

                if len(Up_cat_info) > 0:
                    DA_ncl = sum(Up_cat_info["BasArea"].values)
                    slp_ncl = np.average(
                        Up_cat_info["BasSlope"].values, weights=Up_cat_info["BasArea"].values)
                else:
                    DA_ncl = 0.0
                    slp_ncl = 0.0

            catinfo.loc[idx, "DrainArea"] = DA_ncl + \
                catinfo["BasArea"].values[i]
            catinfo.loc[idx, "Strahler"] = 1
            catinfo.loc[idx, "Seg_order"] = 1
            catinfo.loc[idx, "Seg_ID"] = iseg
            # head watershed
            # if "DA_Chn_L" in catinfo.columns:
            #     catinfo.loc[idx, "DA_Chn_L"] = -1.2345
            #     catinfo.loc[idx, "DA_Chn_Slp"] = -1.2345
            # catinfo.loc[idx, "RivLength"] = -1.2345
            # catinfo.loc[idx, "RivSlope"] = -1.2345
            # catinfo.loc[idx, "Min_DEM"] = -1.2345
            # catinfo.loc[idx, "FloodP_n"] = -1.2345
            # catinfo.loc[idx, "Ch_n"] = -1.2345
            # catinfo.loc[idx, "Max_DEM"] = -1.2345
            catinfo.loc[idx, "DA_Slope"] = (
                catinfo["BasSlope"].values[i]*catinfo["BasArea"].values[i]+DA_ncl*slp_ncl)/catinfo.loc[idx, "DrainArea"]

            # TC_1,TC_2,TC_3,TC_4,TC_5 = calculate_Tc(catinfo.loc[idx, "DrainArea"],catinfo.loc[idx, "DA_Chn_L"],catinfo.loc[idx, "DA_Chn_Slp"])
            #
            # catinfo.loc[idx, "Tc_1"] = TC_1
            # catinfo.loc[idx, "Tc_2"] = TC_2
            # catinfo.loc[idx, "Tc_3"] = TC_3
            # catinfo.loc[idx, "Tc_4"] = TC_4
            # catinfo.loc[idx, "Tc_5"] = TC_5

            icat = icat + 1
            iseg = iseg + 1

    catlist = np.unique(catlist)
    catlist = catlist[catlist > 0]
    #    print(catlist)
    # Loop for each first reach, until go to reach intersection
    newcatlist = np.full((len(catinfo)), -9)
    inewstart = 0

    for i in range(0, len(catlist)):
        catid = catlist[i]
        F_intersect = 1
        #        print("new start            ",i,catid)
        while F_intersect == 1 and catid > 0:
            Up_Reaches_info = catinfo.loc[catinfo["DowSubId"] == catid].copy()
            cur_Reach_info = catinfo.loc[catinfo["SubId"] == catid].copy()
            curcat_idx = catinfo["SubId"] == catid

            # calculate DA of None connected lakes
            if len(routing_ncl) == 0:
                DA_ncl = 0.0
                slp_ncl = 0.0
            else:
                #                print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
                #                print(catid)
                #                print(routing_ncl)
                Upstreamcats = defcat(routing_ncl, catid)  # alll subuds
                Up_cat_info = catinfo_ncl.loc[
                    catinfo_ncl["SubId"].isin(Upstreamcats)
                ].copy()
                if len(Up_cat_info) > 0:
                    DA_ncl = sum(Up_cat_info["BasArea"].values)
                    slp_ncl = np.average(
                        Up_cat_info["BasSlope"].values, weights=Up_cat_info["BasArea"].values)
                else:
                    DA_ncl = 0.0
                    slp_ncl = 0.0

            if (
                len(cur_Reach_info) <= 0
            ):  # reach the most downstream of the watersheds
                break

            if len(Up_Reaches_info) == 1:  # only have one upstream
                catinfo.loc[curcat_idx, "DrainArea"] = (
                    cur_Reach_info["BasArea"].values[0]
                    + Up_Reaches_info["DrainArea"].values[0]
                    + DA_ncl
                )
                catinfo.loc[curcat_idx, "Strahler"] = Up_Reaches_info[
                    "Strahler"
                ].values[0]
                catinfo.loc[curcat_idx, "Seg_order"] = (
                    Up_Reaches_info["Seg_order"].values[0] + 1
                )

                catinfo.loc[curcat_idx,
                            "Seg_ID"] = Up_Reaches_info["Seg_ID"].values[0]

                if "DA_Chn_L" in catinfo.columns:
                    catinfo.loc[curcat_idx, "DA_Chn_L"] = np.maximum(
                        Up_Reaches_info["DA_Chn_L"].values[0], 0) + cur_Reach_info["RivLength"].values[0]

                    catinfo.loc[curcat_idx, "DA_Chn_Slp"] = (cur_Reach_info["RivSlope"].values[0] * np.maximum(cur_Reach_info["RivLength"].values[0], 0)
                                                             + Up_Reaches_info["DA_Chn_L"].values[0] * np.maximum(
                                                                 Up_Reaches_info["DA_Chn_Slp"].values[0], 0)
                                                             ) / catinfo.loc[curcat_idx, "DA_Chn_L"]

                    catinfo.loc[curcat_idx, "DA_Slope"] = (cur_Reach_info["BasSlope"].values[0]*cur_Reach_info["BasArea"].values[0]
                                                           + DA_ncl * slp_ncl
                                                           + Up_Reaches_info["DA_Slope"].values[0] * Up_Reaches_info["DrainArea"].values[0]
                                                           )/catinfo.loc[curcat_idx, "DrainArea"]

                # TC_1,TC_2,TC_3,TC_4,TC_5 = calculate_Tc(catinfo.loc[curcat_idx, "DrainArea"],catinfo.loc[curcat_idx, "DA_Chn_L"],catinfo.loc[curcat_idx, "DA_Chn_Slp"])
                #
                # catinfo.loc[curcat_idx, "Tc_1"] = TC_1
                # catinfo.loc[curcat_idx, "Tc_2"] = TC_2
                # catinfo.loc[curcat_idx, "Tc_3"] = TC_3
                # catinfo.loc[curcat_idx, "Tc_4"] = TC_4
                # catinfo.loc[curcat_idx, "Tc_5"] = TC_5

                #                print('1',catid,catinfo.loc[curcat_idx,'DA'].values,catinfo.loc[curcat_idx,'Strahler'].values,catinfo.loc[curcat_idx,'Sub_order'].values)
                catid = int(cur_Reach_info["DowSubId"].values[0])
            else:  # has mutiple upstram
                # all upstream has been processed
                if (np.min(Up_Reaches_info["Strahler"].values) > 0):

                    if 'DA_Chn_L' in Up_Reaches_info.columns:
                        Up_Reaches_info = Up_Reaches_info.sort_values(
                            by='DA_Chn_L', ascending=False)

                    catinfo.loc[curcat_idx, "DrainArea"] = (
                        cur_Reach_info["BasArea"].values[0]
                        + np.sum(Up_Reaches_info["DrainArea"].values)
                        + DA_ncl
                    )

                    if "DA_Chn_L" in catinfo.columns:
                        catinfo.loc[curcat_idx, "DA_Chn_L"] = np.maximum(
                            Up_Reaches_info["DA_Chn_L"].values[0], 0) + cur_Reach_info["RivLength"].values[0]

                        catinfo.loc[curcat_idx, "DA_Chn_Slp"] = (cur_Reach_info["RivSlope"].values[0] * np.maximum(cur_Reach_info["RivLength"].values[0], 0)
                                                                 + Up_Reaches_info["DA_Chn_L"].values[0] * np.maximum(
                                                                     Up_Reaches_info["DA_Chn_Slp"].values[0], 0)
                                                                 ) / catinfo.loc[curcat_idx, "DA_Chn_L"]

                        catinfo.loc[curcat_idx, "DA_Slope"] = (cur_Reach_info["BasSlope"].values[0]*cur_Reach_info["BasArea"].values[0]
                                                               + DA_ncl * slp_ncl
                                                               + Up_Reaches_info["DA_Slope"].values[0] * Up_Reaches_info["DrainArea"].values[0]
                                                               )/catinfo.loc[curcat_idx, "DrainArea"]

                    # TC_1,TC_2,TC_3,TC_4,TC_5 = calculate_Tc(catinfo.loc[curcat_idx, "DrainArea"],catinfo.loc[curcat_idx, "DA_Chn_L"],catinfo.loc[curcat_idx, "DA_Chn_Slp"])
                    #
                    # catinfo.loc[curcat_idx, "Tc_1"] = TC_1
                    # catinfo.loc[curcat_idx, "Tc_2"] = TC_2
                    # catinfo.loc[curcat_idx, "Tc_3"] = TC_3
                    # catinfo.loc[curcat_idx, "Tc_4"] = TC_4
                    # catinfo.loc[curcat_idx, "Tc_5"] = TC_5

                    if np.min(Up_Reaches_info["Strahler"].values) == np.max(
                        Up_Reaches_info["Strahler"].values
                    ):  # two reach has the same order
                        catinfo.loc[curcat_idx, "Strahler"] = (
                            Up_Reaches_info["Strahler"].values[0] + 1
                        )
                        catinfo.loc[curcat_idx, "Seg_order"] = 1
                        catinfo.loc[curcat_idx, "Seg_ID"] = iseg + 1
                        iseg = iseg + 1
                    #                        print('2',catid,catinfo.loc[catinfo['SubId'] == catid,'DA'].values,catinfo.loc[catinfo['SubId'] == catid,'Strahler'].values,catinfo.loc[catinfo['SubId'] == catid,'Sub_order'].values)
                    else:
                        max_straorder = np.max(
                            Up_Reaches_info["Strahler"].values)
                        catinfo.loc[curcat_idx, "Strahler"] = max_straorder
                        catinfo.loc[curcat_idx, "Seg_order"] = 1
                        catinfo.loc[curcat_idx, "Seg_ID"] = iseg + 1
                        iseg = iseg + 1
                    #                        print('3',catid,catinfo.loc[catinfo['SubId'] == catid,'DA'].values,catinfo.loc[catinfo['SubId'] == catid,'Strahler'].values,catinfo.loc[catinfo['SubId'] == catid,'Sub_order'].values)
                    catid = int(cur_Reach_info["DowSubId"].values[0])
                else:  # there are some reach has not been processed, save id to the list and wait for another loob
                    newcatlist[inewstart] = int(catid)
                    inewstart = inewstart + 1
                    F_intersect = 0

    mask = catinfoall["SubId"].isin(catinfo["SubId"].values)
    catinfoall.loc[mask, "Seg_ID"] = catinfo["Seg_ID"].values
    catinfoall.loc[mask, "Seg_order"] = catinfo["Seg_order"].values
    catinfoall.loc[mask, "Strahler"] = catinfo["Strahler"].values
    catinfoall.loc[mask, "Seg_ID"] = catinfo["Seg_ID"].values
    catinfoall.loc[mask, "DrainArea"] = catinfo["DrainArea"].values
    if "DA_Chn_L" in catinfoall.columns:
        catinfoall.loc[mask, "DA_Chn_L"] = catinfo["DA_Chn_L"].values
        catinfoall.loc[mask, "DA_Slope"] = catinfo["DA_Slope"].values
        catinfoall.loc[mask, "DA_Chn_Slp"] = catinfo["DA_Chn_Slp"].values
    catinfoall.loc[mask, "RivLength"] = catinfo["RivLength"].values
    catinfoall.loc[mask, "RivSlope"] = catinfo["RivSlope"].values

    # catinfoall.loc[mask, "Tc_1"] = catinfo["Tc_1"].values
    # catinfoall.loc[mask, "Tc_2"] = catinfo["Tc_2"].values
    # catinfoall.loc[mask, "Tc_3"] = catinfo["Tc_3"].values
    # catinfoall.loc[mask, "Tc_4"] = catinfo["Tc_4"].values
    # catinfoall.loc[mask, "Tc_5"] = catinfo["Tc_5"].values

    # calcuate channel manning's coefficient
    for i in range(0, len(catinfoall)):
        idx = catinfoall.index[i]
        # if catinfo['BkfWidth'].values[i] > 0 and catinfo['RivSlope'].values[i] > 0 :
        #     catinfo.loc[idx,'Ch_n'] = calculateChannaln(catinfo['BkfWidth'].values[i],catinfo['BkfDepth'].values[i],
        #                       catinfo['Q_Mean'].values[i],catinfo['RivSlope'].values[i])
        if catinfoall[Gauge_col_Name].values[i] > 0:
            if catinfoall["DA_Obs"].values[i] > 0:
                catinfoall.loc[idx, "DA_error"] = (
                    catinfoall["DrainArea"].values[i] / 1000.0 / 1000.0
                ) / catinfoall["DA_Obs"].values[i]

        # calcuate drainage area of each non connected lake

        if catinfoall["Lake_Cat"].values[i] == 2:
            catid = catinfoall["SubId"].values[i]

            # ncl subs to this ncl sub
            Upstreamcats = defcat(routing_ncl, catid)
            Up_cat_info = catinfo_ncl.loc[catinfo_ncl["SubId"].isin(
                Upstreamcats)].copy(deep=True)
            DA = sum(Up_cat_info["BasArea"].values)
            catinfoall.loc[idx, "DrainArea"] = DA

            catinfoall.loc[idx, "DA_Slope"] = np.average(
                Up_cat_info["BasSlope"].values, weights=Up_cat_info["BasArea"].values)

            # catinfoall.loc[idx, "RivLength"] = -1.2345
            # catinfoall.loc[idx, "RivSlope"] = -1.2345
            # catinfoall.loc[idx, "Min_DEM"] = -1.2345
            # catinfoall.loc[idx, "Max_DEM"] = -1.2345
            # catinfoall.loc[idx, "FloodP_n"] = -1.2345
            # catinfoall.loc[idx, "Ch_n"] = -1.2345

            # if "DA_Chn_L" in catinfoall.columns:
            #     catinfoall.loc[idx, "DA_Chn_Slp"] = -1.2345
            #     catinfoall.loc[idx, "DA_Chn_L"] = -1.2345

            # TC_1,TC_2,TC_3,TC_4,TC_5 = calculate_Tc(catinfoall.loc[idx, "DrainArea"],catinfoall.loc[idx, "DA_Chn_L"],catinfoall.loc[idx, "DA_Chn_Slp"])
            #
            # catinfoall.loc[idx, "Tc_1"] = TC_1
            # catinfoall.loc[idx, "Tc_2"] = TC_2
            # catinfoall.loc[idx, "Tc_3"] = TC_3
            # catinfoall.loc[idx, "Tc_4"] = TC_4
            # catinfoall.loc[idx, "Tc_5"] = TC_5

    return catinfoall


def defcat(out, outletid):
    """Functions will return upstream ids in out taht drainage
        to outletid
    ----------

    Notes
    -------

    Returns:
    -------
        Shedid
    """
    otsheds = np.full((1, 1), outletid)
    Shedid = np.full((len(out) + 1, 1), -99999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((len(out) + 1, 1), -99999999999)
        poshdid = 0
        #        print("################################################a")
        for i in range(0, len(otsheds)):
            #            print(otsheds)
            #            print(psid,outletid)
            Shedid[psid] = otsheds[i]
            #            print(Shedid[psid],otsheds[i])
            #            print("##################################################b")
            psid = psid + 1
            irow = np.argwhere(rout[:, 1] == otsheds[i]).astype(int)
            #            print(len(irow))
            for j in range(0, len(irow)):
                # if the catchment id already processed skip
                if rout[irow[j], 0] in Shedid:
                    continue
                noutshd[poshdid] = rout[irow[j], 0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd >= 0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid >= 0]
    return Shedid


def Connect_SubRegion_Update_DownSubId(AllCatinfo, DownCatinfo, Sub_Region_info):
    """Modify outlet subbasin's downstream subbasin ID for each subregion
    Parameters
    ----------
    AllCatinfo                        : dataframe
        it is a dataframe of the attribute table readed from finalcat_info
        .shp or finalrivply_info.shp
    DownCatinfo                       : dataframe
        It is a dataframe includes two columns:
        'Sub_Reg_ID': the subregion id
        'Dow_Sub_Reg_Id': downstream subbasin id of the outlet subbasin in
        this subregion
    Sub_Region_info                   : dataframe
        It is a dataframe includes subregion informations, with following
        columns:
        'Sub_Reg_ID' : the subregion id
        'Dow_Sub_Reg_Id': the downstream subregion id

    Notes
    -------

    Returns:
    -------
        Sub_Region_info      : dataframe
            An new columns indicate the outlet subbasin id of the sub region
            will be added
        AllCatinfo           : dataframe
            Downstream subbasin ID of each subregion's outlet subbasin will
            be updated.
    """

    # update downsteam subbasin id for each subregion outlet subbasins
    # current value is -1, updated to connected downstream subbasin ID
    # loop for each subregion
    routing_info = (
        Sub_Region_info[["Sub_Reg_ID", "Dow_Sub_Reg_Id"]].astype(
            "float").values
    )
    for i in range(0, len(Sub_Region_info)):

        # obtain all subbasin data for i isubregion
        isubregion = Sub_Region_info["Sub_Reg_ID"].values[i]
        # get link id for corresponding down stream sub region inlet point
        ILpt_ID = Sub_Region_info["ILpt_ID"].values[i]

        Sub_Region_cat_info = AllCatinfo.loc[
            AllCatinfo["Region_ID"] == isubregion
        ].copy()

        Sub_Region_info.loc[
            Sub_Region_info["Sub_Reg_ID"] == isubregion, "N_Up_SubRegion"
        ] = len(defcat(routing_info, isubregion))

        # check if this subregion exist in merged data
        if len(Sub_Region_cat_info) <= 0:
            continue

        # Subregion outlet subbasin ID
        outlet_mask = Sub_Region_cat_info["DowSubId"] == -1
        iReg_Outlet_Subid = Sub_Region_cat_info.loc[outlet_mask,
                                                    "SubId"].values[0]
        # (isubregion - 80000) * 200000 - 1
        Sub_Region_info.loc[
            Sub_Region_info["Sub_Reg_ID"] == isubregion, "Outlet_SubId"
        ] = iReg_Outlet_Subid
        # find downstream subregion id of current subregion
        Dow_Sub_Region_id = Sub_Region_info["Dow_Sub_Reg_Id"].values[i]

        # if this region has no down subregions. do not need to
        # modify the dowsubid of the outlet subbasin of this subregion
        if Dow_Sub_Region_id < 0:
            continue

        # find downstrem subbasin id of outlet subbasin
        Down_Sub_info = DownCatinfo.loc[DownCatinfo["ILpt_ID"] == ILpt_ID].copy(
        )

        if len(Down_Sub_info) == 1 and Dow_Sub_Region_id != 79999:
            DownSubid = Down_Sub_info["SubId"].values[0]
            AllCatinfo.loc[
                AllCatinfo["SubId"] == iReg_Outlet_Subid, "DowSubId"
            ] = DownSubid
        # two subregion drainage to the same downstream subregion,
        # the Down_Sub_info only contains one upper subregion
        # the rest do not exist in Down_Sub_info
        elif Dow_Sub_Region_id == 79999:
            AllCatinfo.loc[AllCatinfo["SubId"] ==
                           iReg_Outlet_Subid, "DowSubId"] = -1
        else:
            print("##################################################")
            print("Subregion : ", isubregion, "   To  ", Dow_Sub_Region_id)
            print(Dow_Sub_Region_id, ILpt_ID)
            print(Down_Sub_info)
            print(iReg_Outlet_Subid)
            print(DownCatinfo[["ILpt_ID", "SubId"]])
            print("Need be manually connected")
            print("##################################################")
    return AllCatinfo, Sub_Region_info


def Update_DA_Strahler_For_Combined_Result(AllCatinfo, Sub_Region_info, k, c):
    """Update Drainage area, strahler order of subbains

    Update drainage area and strahler order for combined routing product
    of each subregions
    Parameters
    ----------
    AllCatinfo                        : dataframe
        it is a dataframe of the attribute table readed from finalcat_info
        .shp or finalrivply_info.shp
    Sub_Region_info                   : dataframe
        It is a dataframe includes subregion informations, with following
        columns:
        'Sub_Reg_ID' : the subregion id
        'Dow_Sub_Reg_Id': the downstream subregion id

    Notes
    -------

    Returns:
    -------
    AllCatinfo           : dataframe
        Downstream DA and strahler order of each subregion along the flow
        pathway between subregions will be updated.
    """
    # start from head subregions with no upstream subregion
    Subregion_list = Sub_Region_info[Sub_Region_info["N_Up_SubRegion"] == 1][
        "Sub_Reg_ID"
    ].values
    Subregion_list = np.unique(Subregion_list)
    Subregion_list = Subregion_list.tolist()

    # loop until Subregion_list has no subregions
    # Subregion_list will be updated with downstream subregions of
    # current subregion in Subregion_list
    if len(AllCatinfo.loc[AllCatinfo["DowSubId"] == -1]) > 1:
        print("Wathersed has multiple outlet  ")
        print(AllCatinfo.loc[AllCatinfo["DowSubId"] == -1])
        return AllCatinfo
    elif len(AllCatinfo.loc[AllCatinfo["DowSubId"] == -1]) == 0:
        print("Watershed has no outlet")
        return AllCatinfo
    else:
        Watershedoutletsubid = (
            AllCatinfo.loc[AllCatinfo["DowSubId"]
                           == -1]["SubId"].values[0].astype(int)
        )

    # Area and DA check
    Total_DA_Subregion = 0.0
    for i in range(0, len(Sub_Region_info)):
        Outletsubid_csubr = Sub_Region_info["Outlet_SubId"].values[i]
        Total_DA_Subregion = (
            Total_DA_Subregion
            + AllCatinfo.loc[AllCatinfo["SubId"] ==
                             Outletsubid_csubr]["DrainArea"].values[0]
        )
        print(
            "######Area and DA check for subregion ",
            Sub_Region_info["Sub_Reg_ID"].values[i],
        )
        print(
            "DA at the subregion outlet is    ",
            AllCatinfo.loc[AllCatinfo["SubId"] ==
                           Outletsubid_csubr]["DrainArea"].values[0],
        )
        print(
            "Total subregion area is          ",
            sum(
                AllCatinfo.loc[
                    AllCatinfo["Region_ID"] == Sub_Region_info["Sub_Reg_ID"].values[i]
                ]["BasArea"].values
            ),
        )

    iloop = 1
    while len(Subregion_list) > 0:
        print("loop  ", iloop)
        print(Subregion_list)
        current_loop_list = copy.copy(Subregion_list)
        Subregion_list = []
        # loop for current subregion list
        for i in range(0, len(current_loop_list)):
            # current subregion id
            c_subr_id = current_loop_list[i]
            # down subregion id of current subregion
            if c_subr_id == 79999:
                continue

            d_subr_id = Sub_Region_info[Sub_Region_info["Sub_Reg_ID"] == c_subr_id][
                "Dow_Sub_Reg_Id"
            ].values[0]
            # add down subregon id to the list for next while loop
            Subregion_list.append(d_subr_id)

            # obtain outlet subid of this region
            Outletsubid_csubr = Sub_Region_info[
                Sub_Region_info["Sub_Reg_ID"] == c_subr_id
            ]["Outlet_SubId"].values[0]
            # obtain outlet sub info
            Outletsub_info = AllCatinfo.loc[
                AllCatinfo["SubId"] == Outletsubid_csubr
            ].copy()
            # obtain down subid of the outlet subbasin
            downsubid = Outletsub_info["DowSubId"].values[0]

            # downsubid did not exist.....
            if len(AllCatinfo.loc[AllCatinfo["SubId"] == downsubid]) <= 0:
                if int(c_subr_id) != int(Watershedoutletsubid):
                    print("Subregion:   ", c_subr_id)
                    print(
                        "SubId is ",
                        Outletsubid_csubr,
                        " DownSubId is  ",
                        downsubid,
                        Watershedoutletsubid,
                        int(c_subr_id) != int(Watershedoutletsubid),
                    )
                continue

            downsub_reg_id = AllCatinfo.loc[AllCatinfo["SubId"] == downsubid][
                "Region_ID"
            ].values[0]

            print("Subregion:   ", c_subr_id,
                  downsub_reg_id, d_subr_id, downsubid)

            if downsub_reg_id != d_subr_id:
                print(
                    "Subregion:   ",
                    c_subr_id,
                    "  did not connected with    ",
                    d_subr_id,
                )
                continue

            while downsub_reg_id == d_subr_id:
                csubid = downsubid  # update DA and Strahler for this subbasin
                print("Subregion:   ", c_subr_id, downsub_reg_id,
                      d_subr_id, downsubid, csubid)
                # current subid info
                C_sub_info = AllCatinfo.loc[AllCatinfo["SubId"] == csubid].copy(
                )
                # find all subbasin drainge to this csubid
                Upper_sub_info = AllCatinfo.loc[AllCatinfo["DowSubId"] == csubid].copy(
                )

                # new DA = basin area + DA of upper subregions

                NewDA = C_sub_info["BasArea"].values[0] + sum(
                    Upper_sub_info["DrainArea"].values
                )

                # calculate new Strahler orders
                # up stream Strahler orders
                Strahlers = Upper_sub_info["Strahler"].values
                maxStrahler = max(Strahlers)
                if np.sum(Strahlers == maxStrahler) >= 2:
                    newStrahler = maxStrahler + 1
                else:
                    newStrahler = maxStrahler
                # updateAllCatinfo
                AllCatinfo.loc[AllCatinfo["SubId"] ==
                               csubid, "Strahler"] = newStrahler
                AllCatinfo.loc[AllCatinfo["SubId"]
                               == csubid, "DrainArea"] = NewDA

                if AllCatinfo[AllCatinfo["SubId"] == csubid]["DA_Obs"].values[0] > 0:
                    da_obs = AllCatinfo[AllCatinfo["SubId"]
                                        == csubid]["DA_Obs"].values[0]
                    da_sim = AllCatinfo[AllCatinfo["SubId"] ==
                                        csubid]["DrainArea"].values[0]/1000.0/1000.0
                    if da_obs > 0:
                        AllCatinfo.loc[AllCatinfo["SubId"] ==
                                       csubid, "DA_error"] = da_sim/da_obs
                        AllCatinfo.loc[AllCatinfo["SubId"]
                                       == csubid, "Use_region"] = 1

                # da = AllCatinfo.loc[AllCatinfo["SubId"] == csubid, "DrainArea"].values[0] / 1000 / 1000  # m2 to km2
                # q = func_Q_DA(da, k, c)
                # bk_w= 7.2 * q ** 0.5
                # bk_d = 0.27 * q ** 0.3
                # bk_q = q
                # bk_slope = AllCatinfo.loc[AllCatinfo["SubId"] == csubid, "RivSlope"].values[0]
                # n_rch = calculateChannaln(bk_w, bk_d, bk_q, bk_slope)
                # if n_rch < min_manning_n:
                #     n_rch = 0.01
                # if n_rch > max_manning_n:
                #     nrch = max_manning_n
                # AllCatinfo.loc[AllCatinfo["SubId"] == csubid, "BkfWidth"] = 7.2 * q ** 0.5
                # AllCatinfo.loc[AllCatinfo["SubId"] == csubid, "BkfDepth"] = 0.27 * q ** 0.3
                # AllCatinfo.loc[AllCatinfo["SubId"] == csubid, "Q_Mean"] = q
                # AllCatinfo.loc[AllCatinfo["SubId"] == csubid, "Ch_n"] = n_rch

                # find next downsubbasin id
                downsubid = C_sub_info["DowSubId"].values[0]

                # downsubid did not exist.....
                if len(AllCatinfo.loc[AllCatinfo["SubId"] == downsubid]) <= 0:
                    if int(csubid) != int(Watershedoutletsubid):
                        print("Subregion:   ", d_subr_id)
                        print(
                            "SubId is ",
                            csubid,
                            " DownSubId is  ",
                            downsubid,
                            Watershedoutletsubid,
                            int(csubid) != int(Watershedoutletsubid),
                        )
                    break

                downsub_reg_id = AllCatinfo.loc[AllCatinfo["SubId"] == downsubid][
                    "Region_ID"
                ].values[0]

            # update list for next loop
        Subregion_list = list(set(Subregion_list))
        iloop = iloop + 1

    print("Check drainage area:")
    print("Total basin area is              ",
          sum(AllCatinfo["BasArea"].values))
    print(
        "DA of the watersehd outlet is    ",
        AllCatinfo.loc[AllCatinfo["SubId"] == int(Watershedoutletsubid)]["DrainArea"].values[
            0
        ],
    )
    print("Total DA of each subregion       ", Total_DA_Subregion)
    return AllCatinfo


def Determine_Lake_HRU_Id(Attribute_Table):
    """Function to determin hruid after combine lake polygon
    and subbasin polygon
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """

    Sub_ID = "SubId"
    Sub_Lake_ID = "HyLakeId"
    Lake_Id = "Hylak_id"

    # get current maximum subbasin Id
    max_subbasin_id = max(Attribute_Table[Sub_ID].values)
    # total number of new feature by overaly two polygon
    N_new_features = len(Attribute_Table)
    new_hruid = 1
    new_hruid_list = np.full(N_new_features + 100, np.nan)
    old_newhruid_list = np.full(N_new_features + 100, np.nan)

    # The new hru id is determined by following rules
    # 1 if subbasin hylake id < 0, means this feature do not have lake
    #  then, the new hru id  = old subid
    # 2 if subbasin hylke >0 and subbasin hylake = overlayied hyalke id
    #  then this feature is a lake hru
    #  the hru id is the hru id  = old subid + max_subbasin_id + 10
    # 3 if the subbasin hylakeid >0, but subbasin hylake id != overlaied lake id
    #  then this feature is not covered by the lake the new hru id  = old subid

    for i in range(0, len(Attribute_Table)):
        subid_sf_obj = Attribute_Table[Sub_ID].values[i]
        lakelakeid_sf_obj = Attribute_Table[Lake_Id].values[i]
        Sub_Lakeid_sf_obj = Attribute_Table[Sub_Lake_ID].values[i]

        # the null value in the attribute table is not convertable try and set
        # to -1
        try:
            subid_sf = float(subid_sf_obj)
        except:
            subid_sf = -1
            Attribute_Table.loc[i, Sub_ID] = -1
            pass

        try:
            lakelakeid_sf = float(lakelakeid_sf_obj)
        except:
            lakelakeid_sf = -1
            Attribute_Table.loc[i, Lake_Id] = -1
            pass

        try:
            Sub_Lakeid_sf = float(Sub_Lakeid_sf_obj)
        except:
            Sub_Lakeid_sf = -1
            Attribute_Table.loc[i, Sub_Lake_ID] = -1
            pass

        # first skip feature with subbasin id < 0
        # deal with this later
        if subid_sf < 0:
            continue

        # feature is not lake
        if Sub_Lakeid_sf <= 0:
            old_new_hruid = subid_sf
            Attribute_Table.loc[i, "HRU_IsLake"] = -1

        if Sub_Lakeid_sf > 0:
            if (
                lakelakeid_sf == Sub_Lakeid_sf
            ):  # the lakeid from lake polygon = lake id in subbasin polygon
                old_new_hruid = subid_sf + max_subbasin_id + 10
                Attribute_Table.loc[i, "HRU_IsLake"] = 1
            else:  # the lakeid from lake polygon != lake id in subbasin polygon, this lake do not belong to this subbasin, this part of subbasin treat as non lake hru
                old_new_hruid = float(subid_sf)
                Attribute_Table.loc[i, "HRU_IsLake"] = -1

        # if it is a new hru id
        if old_new_hruid not in old_newhruid_list:
            Attribute_Table.loc[i, "HRULake_ID"] = new_hruid
            old_newhruid_list[new_hruid] = old_new_hruid
            new_hruid = new_hruid + 1
        # if it is not new hru id
        else:
            Attribute_Table.loc[i, "HRULake_ID"] = int(
                np.argwhere(old_newhruid_list == old_new_hruid)[0]
            )

    # deal with feature with negative subbasin id
    # use the Sub_Lake_ID to find subbasin id,
    # if Sub_Lake_ID from lake polygon is also sammller than zero
    #    report an error
    for i in range(0, len(Attribute_Table)):
        subid_sf = Attribute_Table[Sub_ID].values[i]
        lakelakeid_sf = Attribute_Table[Lake_Id].values[i]
        Sub_Lakeid_sf = Attribute_Table[Sub_Lake_ID].values[i]
        if subid_sf > 0:
            continue
        if lakelakeid_sf <= 0:
            print("lake has unexpected holes ")
            print(Attribute_Table.loc[i, :])
            continue
        # find the subid of lakelakeid_sf
        tar_info = Attribute_Table.loc[
            (Attribute_Table[Lake_Id] == lakelakeid_sf)
            & (Attribute_Table["HRU_IsLake"] > 0)
        ]
        if len(tar_info) <= 0:
            print('Lakeid do not exist   ', lakelakeid_sf)
            continue
        Attribute_Table.loc[i, Sub_ID] = tar_info[Sub_ID].values[0]
        Attribute_Table.loc[i, "HRU_IsLake"] = tar_info["HRU_IsLake"].values[0]
        Attribute_Table.loc[i, "HRULake_ID"] = tar_info["HRULake_ID"].values[0]
        Attribute_Table.loc[i, Sub_Lake_ID] = tar_info[Sub_Lake_ID].values[0]

    return Attribute_Table


def Retrun_Validate_Attribute_Value(Attri_table_Lake_HRU_i, SubInfo, Col_NM, info_data):

    info_lake_hru = Attri_table_Lake_HRU_i.loc[Attri_table_Lake_HRU_i[Col_NM] != 0]
    if (
        len(info_lake_hru) > 0
    ):  # other part of this lake hru has validate soilid use the one with maximum area
        Updatevalue = info_lake_hru[Col_NM].values[0]
    else:  # check if the subbasin has a valid soil id
        SubInfo = SubInfo.loc[SubInfo[Col_NM] != 0]
        SubInfo = SubInfo.sort_values(by="HRU_Area", ascending=False)
        if len(SubInfo) > 0:
            Updatevalue = SubInfo[Col_NM].values[0]
        else:
            Updatevalue = info_data.loc[info_data[Col_NM]
                                        != 0][Col_NM].values[0]
    #            print("asdfasdfsadf2",Updatevalue)
    #    print("asdfasdfsadf1",Updatevalue)
    return Updatevalue


def Determine_HRU_Attributes(
    Attri_table,
    Sub_ID,
    Landuse_ID,
    Soil_ID,
    Veg_ID,
    Other_Ply_ID_1,
    Other_Ply_ID_2,
    Landuse_info_data,
    Soil_info_data,
    Veg_info_data,
):
    """Function to determine landuse,soil,and veg and other properties after union all polygons
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """

    # find invald rows
    inval_sub = Attri_table[Sub_ID] <= 0
    inval_landuse = (Attri_table[Landuse_ID] <= 0) & (
        Attri_table[Landuse_ID] != -1)
    inval_soil = (Attri_table[Soil_ID] <= 0) & (Attri_table[Soil_ID] != -1)
    inval_veg = (Attri_table[Veg_ID] <= 0) & (Attri_table[Veg_ID] != -1)
    inval_o2 = (Attri_table[Other_Ply_ID_2] <= 0) & (
        Attri_table[Other_Ply_ID_2] != -1)

    inval_rows = inval_sub | inval_landuse | inval_soil | inval_veg | inval_o2
    Lake_HRU_IDS = np.unique(Attri_table.loc[inval_rows, "HRULake_ID"].values)
    Lake_HRU_IDS = Lake_HRU_IDS[Lake_HRU_IDS > 0]
    # landuse,soil,and veg and other properties for lake hrus

    for i in range(0, len(Lake_HRU_IDS)):
        ilake_hru_id = Lake_HRU_IDS[i]
        if ilake_hru_id == 0:
            continue
        # obtain i lake hru attributes
        Attri_table_Lake_HRU_i = Attri_table.loc[
            Attri_table["HRULake_ID"] == ilake_hru_id
        ]
        Attri_table_Lake_HRU_i = Attri_table_Lake_HRU_i.sort_values(
            by="HRU_Area", ascending=False
        )
        Subid = Attri_table_Lake_HRU_i["SubId"].values[0]
        SubInfo = Attri_table.loc[Attri_table["SubId"] == Subid]
        for j in range(0, len(Attri_table_Lake_HRU_i)):
            ihru_id = Attri_table_Lake_HRU_i["HRU_ID"].values[j]
            is_Lake = Attri_table_Lake_HRU_i["HRU_IsLake"].values[j]
            if is_Lake > 0:

                if (
                    Attri_table_Lake_HRU_i[Soil_ID].values[0] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Soil_ID, Soil_info_data
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Soil_ID
                    ] = Vali_Value
                else:
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Soil_ID
                    ] = Attri_table_Lake_HRU_i[Soil_ID].values[0]

                if (
                    Attri_table_Lake_HRU_i[Other_Ply_ID_1].values[0] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Other_Ply_ID_1, Attri_table
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Other_Ply_ID_1
                    ] = Vali_Value
                else:
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Other_Ply_ID_1
                    ] = Attri_table_Lake_HRU_i[Other_Ply_ID_1].values[0]

                if (
                    Attri_table_Lake_HRU_i[Other_Ply_ID_2].values[0] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Other_Ply_ID_2, Attri_table
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Other_Ply_ID_2
                    ] = Vali_Value
                else:
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Other_Ply_ID_2
                    ] = Attri_table_Lake_HRU_i[Other_Ply_ID_2].values[0]

                Attri_table.loc[Attri_table["HRU_ID"]
                                == ihru_id, Landuse_ID] = int(-1)
                Attri_table.loc[Attri_table["HRU_ID"]
                                == ihru_id, Veg_ID] = int(-1)

            else:
                if (
                    Attri_table_Lake_HRU_i[Soil_ID].values[j] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Soil_ID, Soil_info_data
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Soil_ID
                    ] = Vali_Value
                if (
                    Attri_table_Lake_HRU_i[Landuse_ID].values[j] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Landuse_ID, Landuse_info_data
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Landuse_ID
                    ] = Vali_Value
                if (
                    Attri_table_Lake_HRU_i[Veg_ID].values[j] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Veg_ID, Veg_info_data
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Veg_ID
                    ] = Vali_Value
                if (
                    Attri_table_Lake_HRU_i[Other_Ply_ID_1].values[j] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Other_Ply_ID_1, Attri_table
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Other_Ply_ID_1
                    ] = Vali_Value
                if (
                    Attri_table_Lake_HRU_i[Other_Ply_ID_2].values[j] == 0
                ):  # the current hru has invalidate soil id
                    Vali_Value = Retrun_Validate_Attribute_Value(
                        Attri_table_Lake_HRU_i, SubInfo, Other_Ply_ID_2, Attri_table
                    )
                    Attri_table.loc[
                        Attri_table["HRU_ID"] == ihru_id, Other_Ply_ID_2
                    ] = Vali_Value

    Attri_table = Attri_table.drop(
        columns=['LAND_USE_C', 'VEG_C', 'SOIL_PROF'])
    Attri_table = pd.merge(Attri_table, Landuse_info_data,
                           how='inner', on=Landuse_ID).copy(deep=True)
    Attri_table = pd.merge(Attri_table, Soil_info_data,
                           how='inner', on=Soil_ID).copy(deep=True)
    Attri_table = pd.merge(Attri_table, Veg_info_data,
                           how='inner', on=Veg_ID).copy(deep=True)

    # for i in range(0, len(Attri_table)):
    #     if Attri_table["HRULake_ID"].values[i] == 0:
    #         continue
    #     Is_lake_hru = Attri_table["HRU_IsLake"].values[i]
    #     lake_hru_ID = Attri_table["HRULake_ID"].values[i]
    #     hruid = Attri_table["HRU_ID"].values[i]
    #     Landuse_ID_num = Attri_table[Landuse_ID].values[i]
    #     Soil_ID_num = Attri_table[Soil_ID].values[i]
    #     Veg_ID_num = Attri_table[Veg_ID].values[i]
    #
    #     if hruid == 0 or Landuse_ID_num == 0 or Soil_ID_num == 0 or Veg_ID_num == 0:
    #         continue
    #
    #     if not isinstance(hruid, numbers.Number) and hruid != hruid:
    #         continue
    #
    #     val = Attri_table.loc[Attri_table["HRU_ID"] == hruid][Landuse_ID].values[0]
    #     if isinstance(val,numbers.Number) and val == val:
    #         Attri_table.loc[i, Landuse_ID] = int(
    #             Attri_table.loc[Attri_table["HRU_ID"] == hruid][Landuse_ID].values[0]
    #         )
    #     else:
    #         Attri_table.loc[i, Landuse_ID] = 0
    #         Attri_table.loc[i, "SubId"] = int(0)
    #
    #     val = Attri_table.loc[Attri_table["HRU_ID"] == hruid][Veg_ID].values[0]
    #     if isinstance(val,numbers.Number) and val == val:
    #         Attri_table.loc[i, Veg_ID] = int(
    #             Attri_table.loc[Attri_table["HRU_ID"] == hruid][Veg_ID].values[0]
    #         )
    #     else:
    #         Attri_table.loc[i, Veg_ID] = 0
    #         Attri_table.loc[i, "SubId"] = int(0)
    #     val = Attri_table.loc[Attri_table["HRU_ID"] == hruid][Soil_ID].values[0]
    #     if isinstance(val,numbers.Number) and val == val:
    #         Attri_table.loc[i, Soil_ID] = int(
    #             Attri_table.loc[Attri_table["HRU_ID"] == hruid][Soil_ID].values[0]
    #         )
    #     else:
    #         Attri_table.loc[i, Soil_ID] = 0
    #         Attri_table.loc[i, "SubId"] = int(0)
    #     val = Attri_table.loc[Attri_table["HRU_ID"] == hruid][Other_Ply_ID_1].values[0]
    #     if isinstance(val,numbers.Number) and val == val:
    #         Attri_table.loc[i, Other_Ply_ID_1] = int(
    #             Attri_table.loc[Attri_table["HRU_ID"] == hruid][Other_Ply_ID_1].values[0]
    #         )
    #     else:
    #         Attri_table.loc[i, Other_Ply_ID_1] = 0
    #         Attri_table.loc[i, "SubId"] = int(0)
    #     val = Attri_table.loc[Attri_table["HRU_ID"] == hruid][Other_Ply_ID_2].values[0]
    #     if isinstance(val,numbers.Number) and val == val:
    #
    #         Attri_table.loc[i, Other_Ply_ID_2] = int(
    #             Attri_table.loc[Attri_table["HRU_ID"] == hruid][Other_Ply_ID_2].values[0]
    #         )
    #     else:
    #         Attri_table.loc[i, Other_Ply_ID_2] = 0
    #         Attri_table.loc[i, "SubId"] = int(0)
    #
    #     if isinstance(Landuse_ID_num, numbers.Number) and Landuse_ID_num ==Landuse_ID_num:
    #         Attri_table.loc[i, "LAND_USE_C"] = Landuse_info_data.loc[
    #             Landuse_info_data[Landuse_ID] == int(Landuse_ID_num), "LAND_USE_C"
    #         ].values[0]
    #     else:
    #         Attri_table.loc[i, "LAND_USE_C"] = int(0)
    #         Attri_table.loc[i, "SubId"] = int(0)
    #
    #     if isinstance(Veg_ID_num, numbers.Number) and Veg_ID_num == Veg_ID_num:
    #         Attri_table.loc[i, "VEG_C"] = Veg_info_data.loc[
    #             Veg_info_data[Veg_ID] == int(Veg_ID_num), "VEG_C"
    #         ].values[0]
    #     else:
    #         Attri_table.loc[i, "VEG_C"] = int(0)
    #         Attri_table.loc[i, "SubId"] = int(0)
    #
    #     if isinstance(Soil_ID_num, numbers.Number) and Soil_ID_num == Soil_ID_num:
    #         Attri_table.loc[i, "SOIL_PROF"] = (
    #             Soil_info_data.loc[
    #                 Soil_info_data[Soil_ID] == int(Soil_ID_num), "SOIL_PROF"
    #             ].values[0]
    #         )
    #     else:
    #         Attri_table.loc[i, "SOIL_PROF"] = int(0)
    #         Attri_table.loc[i, "SubId"] = int(0)

    Attri_table["facters"] = (
        Attri_table["HRULake_ID"].astype(str)
        + Attri_table[Landuse_ID].astype(str)
        + Attri_table[Soil_ID].astype(str)
        #        + Attri_table[Veg_ID].astype(str)
        + Attri_table[Other_Ply_ID_1].astype(str)
        #        + Attri_table[Other_Ply_ID_2].astype(str)
    )
    Attri_table["HRU_ID_New"] = pd.factorize(Attri_table["facters"])[0] + 1
    return Attri_table


def Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Increase_DA(
    finalriv_info, Conn_Lakes_ply, Area_Min
):
    """Modify attribute table mapoldnew_info, create new sub id for subbasin needs to be merged.
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    sub_colnm = "SubId"
    down_colnm = "DowSubId"
    DA_colnm = "DrainArea"
    SegID_colnm = "Seg_ID"

    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in finalriv_info.columns:
        Gauge_col_Name = "Has_Gauge"

    # Modify attribute table mapoldnew_info, create new sub id for
    # 1. catchment needs to be merged due to meet the drainage area thresthold
    # 2. connected lake catchment are transfered into non-connected lake catchment

    # copy attribute tabke and create two column to store new subid
    mapoldnew_info = finalriv_info.copy(deep=True)
    mapoldnew_info["nsubid"] = -1
    mapoldnew_info["ndownsubid"] = -1
    mapoldnew_info.reset_index(drop=True, inplace=True)
    max_da = np.max(mapoldnew_info['DrainArea'].values)

    # stupid larege drainage area are provied
    if Area_Min * 1000 * 1000 > max_da:
        Area_Min = max_da/1000.0/1000.0 - 1
        mapoldnew_info[Gauge_col_Name] = -1
        finalriv_info[Gauge_col_Name] = -1

    # select catchment segment that meet the drainage area
    Selected_riv = finalriv_info.loc[
        finalriv_info[DA_colnm] >= Area_Min * 1000 * 1000
    ].copy()  # find river with drainage area larger than area thresthold

    Selected_riv = update_the_selected_river_to_connect_upsub_with_largest_da(
        Selected_riv, mapoldnew_info)

    # calcuate topology for the selected catchment
    Selected_riv = UpdateTopology(Selected_riv, UpdateSubId=-1)
    Selected_riv = Selected_riv.sort_values(
        ["Strahler"], ascending=(True)
    )  # sort selected river by Strahler stream order
    Selected_riv = Selected_riv.loc[Selected_riv['Lake_Cat'] != 2].copy()
    Subid_main = Selected_riv[sub_colnm].values

    # Obtain connected lakes based on current river segment
    Connected_Lake_Mainriv = Selected_riv.loc[Selected_riv["Lake_Cat"] == 1][
        "HyLakeId"
    ].values
    Connected_Lake_Mainriv = np.unique(
        Connected_Lake_Mainriv[Connected_Lake_Mainriv > 0]
    )
    Lakecover_riv = finalriv_info.loc[
        finalriv_info["HyLakeId"].isin(Connected_Lake_Mainriv)
    ].copy()
    Subid_lakes = Lakecover_riv[sub_colnm].values
    #####

    ###
    # add removed gauges and remove gauge if gauge located within a lake on the main river
    unselected_gauges_subids_info = finalriv_info.loc[
        (~finalriv_info["SubId"].isin(Subid_main)) &
        (finalriv_info[Gauge_col_Name] > 0)
    ].copy(deep=True)
    unselected_gauges_subids = unselected_gauges_subids_info.loc[
        unselected_gauges_subids_info["HyLakeId"] <= 0]["SubId"].values

    finalriv_info_ncl = finalriv_info.copy(deep=True)
    # make unselected gauge to be a false lake
    mask1 = finalriv_info_ncl['SubId'].isin(unselected_gauges_subids)
    mask2 = finalriv_info_ncl['HyLakeId'] < 1
    mask = np.logical_and(mask1, mask2)
    finalriv_info_ncl.loc[mask, 'HyLakeId'] = - \
        finalriv_info_ncl.loc[finalriv_info_ncl['SubId'].isin(
            unselected_gauges_subids), 'SubId']
    fake_obs_hyalkeids = finalriv_info_ncl.loc[finalriv_info_ncl['SubId'].isin(
        unselected_gauges_subids), 'HyLakeId'].values
    ##
    ###
    # identify which connected lake will be moved to non connected lake due to remove of river network
    All_Conn_Lakeids = Conn_Lakes_ply["Hylak_id"].values
    mask = np.in1d(All_Conn_Lakeids, Connected_Lake_Mainriv)
    Conn_To_NonConlakeids = All_Conn_Lakeids[np.logical_not(mask)].copy()
    Conn_To_NonConlake_info = finalriv_info_ncl.loc[
        (finalriv_info_ncl["HyLakeId"].isin(Conn_To_NonConlakeids)) |
        (finalriv_info_ncl["HyLakeId"].isin(fake_obs_hyalkeids))
    ].copy()

    # Get origional non connected lake subid and lake id
    Old_Non_Connect_SubIds = finalriv_info.loc[finalriv_info["Lake_Cat"] == 2][
        "SubId"
    ].values
    Old_Non_Connect_SubIds = np.unique(
        Old_Non_Connect_SubIds[Old_Non_Connect_SubIds > 0]
    )
    Old_Non_Connect_LakeIds = finalriv_info.loc[finalriv_info["Lake_Cat"] == 2][
        "HyLakeId"
    ].values
    Old_Non_Connect_LakeIds = np.unique(
        Old_Non_Connect_LakeIds[Old_Non_Connect_LakeIds > 0]
    )

    # obtain rivsegments that covered by remaining lakes
    Selected_riv_ids = np.unique(
        Subid_main
    )  # np.unique(np.concatenate([Subid_main,Subid_lakes]))
    routing_info = finalriv_info[["SubId", "DowSubId"]].astype("float").values
    Seg_IDS = Selected_riv["Seg_ID"].values
    Seg_IDS = np.unique(Seg_IDS)

    # for Non connected lakes, catchment polygon do not need to change
    mapoldnew_info.loc[mapoldnew_info["Lake_Cat"] == 2, "nsubid"] = mapoldnew_info.loc[
        mapoldnew_info["Lake_Cat"] == 2
    ]["SubId"].values

    # for catchment polygon flow to lake with is changed from connected lake to non connected lakes
    idx = (
        Conn_To_NonConlake_info.groupby(
            ["HyLakeId"])["DrainArea"].transform(max)
        == Conn_To_NonConlake_info["DrainArea"]
    )
    Conn_To_NonConlake_info_outlet = Conn_To_NonConlake_info[idx]

    Conn_To_NonConlake_info_outlet = Conn_To_NonConlake_info_outlet.sort_values(
        ["Strahler", "DrainArea"], ascending=[True, True]
    )

    ###
    all_outlet_sub_info = finalriv_info[finalriv_info["DowSubId"] < 0].copy(
        deep=True)

    for i in range(0, len(all_outlet_sub_info)):
        if all_outlet_sub_info['DrainArea'].values[i] <= Area_Min * 1000 * 1000:
            tsubid = all_outlet_sub_info['SubId'].values[i]
            All_up_subids = defcat(routing_info, tsubid)

            mapoldnew_info = New_SubId_To_Dissolve(
                subid=tsubid,
                catchmentinfo=finalriv_info,
                mapoldnew_info=mapoldnew_info,
                ismodifids=1,
                mainriv=finalriv_info,
                modifiidin=All_up_subids,
                Lake_Cat=-1,
            )

    # process fron upstream lake to down stream lake
    for i in range(0, len(Conn_To_NonConlake_info_outlet)):
        processed_subid = np.unique(
            mapoldnew_info.loc[mapoldnew_info["nsubid"] > 0][sub_colnm].values
        )

        C_T_N_Lakeid = Conn_To_NonConlake_info_outlet["HyLakeId"].values[i]
        Lake_Cat_info = finalriv_info_ncl[finalriv_info_ncl["HyLakeId"]
                                          == C_T_N_Lakeid]
        Riv_Seg_IDS = np.unique(Lake_Cat_info["Seg_ID"].values)
        modifysubids = []
        for j in range(0, len(Riv_Seg_IDS)):

            iriv_seg = Riv_Seg_IDS[j]
            Lake_Cat_seg_info = Lake_Cat_info.loc[
                Lake_Cat_info["Seg_ID"] == iriv_seg
            ].copy()
            Lake_Cat_seg_info = Lake_Cat_seg_info.sort_values(
                ["DrainArea"], ascending=(False))
            tsubid = Lake_Cat_seg_info["SubId"].values[0]

            All_up_subids = defcat(routing_info, tsubid)
            All_up_subids = All_up_subids[All_up_subids > 0]

            mask = np.in1d(All_up_subids, processed_subid)
            seg_sub_ids = All_up_subids[np.logical_not(mask)]

            modifysubids = [
                *modifysubids,
                *seg_sub_ids.tolist(),
            ]  # combine two list not sum

        modifysubids = np.asarray(modifysubids)

        mapoldnew_info = New_SubId_To_Dissolve(
            subid=tsubid,
            catchmentinfo=finalriv_info_ncl,
            mapoldnew_info=mapoldnew_info,
            ismodifids=1,
            mainriv=finalriv_info_ncl,
            modifiidin=modifysubids,
            Lake_Cat=2,
        )
        # 3 for rest of the polygons dissolve to main river

    for iseg in range(0, len(Seg_IDS)):
        #            print('#########################################################################################33333')
        i_seg_id = Seg_IDS[iseg]
        i_seg_info = Selected_riv.loc[Selected_riv["Seg_ID"] == i_seg_id].copy(
        )
        i_seg_info = i_seg_info.sort_values(["Seg_order"], ascending=(True))
        sum_area = 0
        modifysubids = []
        seg_order = 1
        for iorder in range(0, len(i_seg_info)):
            tsubid = i_seg_info["SubId"].values[iorder]
            iorder_Lakeid = i_seg_info["HyLakeId"].values[iorder]
            modifysubids.append(tsubid)
            processed_subid = np.unique(
                mapoldnew_info.loc[mapoldnew_info["nsubid"]
                                   > 0][sub_colnm].values
            )
            sum_area = sum_area + i_seg_info["BasArea"].values[iorder]
            # two seg has the same HyLakeId id, can be merged

            # if has the same lake attribute with down stream subasins
            if iorder != len(i_seg_info) - 1:
                con_lake = i_seg_info["HyLakeId"].values[iorder] == i_seg_info["HyLakeId"].values[iorder + 1]
            else:
                con_lake = False
            # if do not have the gauge
            con_gauge = i_seg_info[Gauge_col_Name].values[iorder] <= 0

            # if do not meet the drainage area thresthold
            if iorder <= len(i_seg_info) - 2:
                if i_seg_info["BasArea"].values[iorder + 1] < 10 * 1000 * 1000:
                    con_area = True
                else:
                    con_area = sum_area < Area_Min * 1000 * 1000
            else:
                con_area = sum_area < Area_Min * 1000 * 1000

            con_area_lake = sum_area < 10 * 1000 * 1000
            # conditions check if subbasin is between two lakes
            if iorder >= 1 and iorder <= len(i_seg_info) - 2:
                # if not havve the same attribute with upstream lake
                con_lake_up = i_seg_info["HyLakeId"].values[iorder -
                                                            1] != i_seg_info["HyLakeId"].values[iorder]
                # if not havve the same attribute with down lake
                con_lake_down = i_seg_info["HyLakeId"].values[iorder] != i_seg_info["HyLakeId"].values[iorder + 1]
                # if this is not a lake subbasin
                is_lake = i_seg_info["HyLakeId"].values[iorder] <= 0

                # i_seg_info["BasArea"].values[iorder] < 10 * 1000 * 1000
                con_area_lake = False

            if iorder == 0 and iorder < len(i_seg_info) - 1:
                # if not havve the same attribute with upstream lake
                con_lake_up = True
                # if not havve the same attribute with down lake
                con_lake_down = i_seg_info["HyLakeId"].values[iorder] != i_seg_info["HyLakeId"].values[iorder + 1]
                # if this is not a lake subbasin
                is_lake = i_seg_info["HyLakeId"].values[iorder] <= 0

                # i_seg_info["BasArea"].values[iorder] < 10 * 1000 * 1000
                con_area_lake = False

            if iorder == len(i_seg_info) - 2 and i_seg_info["HyLakeId"].values[iorder] > 0:
                # if not havve the same attribute with upstream lake
                con_lake_up = True
                # if not havve the same attribute with down lake
                con_lake_down = i_seg_info["HyLakeId"].values[iorder] != i_seg_info["HyLakeId"].values[iorder + 1]
                # if this is not a lake subbasin
                is_lake = i_seg_info["HyLakeId"].values[iorder] > 0

                # if tsubid == 9023607:
                #     print(con_lake_up,con_lake_down,is_lake,iorder,len(i_seg_info) - 2,i_seg_info["HyLakeId"].values[iorder] > 0)
                # change add lake to downstream info
                downsubid = i_seg_info["DowSubId"].values[iorder]

                # i_seg_info["BasArea"].values[iorder +1] < 10 * 1000 * 1000
                con_area_lake = False

                if con_area_lake and con_lake_down:
                    finalriv_info.loc[
                        finalriv_info["SubId"] == downsubid, 'HyLakeId'] = i_seg_info["HyLakeId"].values[iorder]
                    finalriv_info.loc[
                        finalriv_info["SubId"] == downsubid, 'Lake_Cat'] = i_seg_info["Lake_Cat"].values[iorder]
                    finalriv_info.loc[
                        finalriv_info["SubId"] == downsubid, 'LakeVol'] = i_seg_info["LakeVol"].values[iorder]
                    finalriv_info.loc[
                        finalriv_info["SubId"] == downsubid, 'LakeDepth'] = i_seg_info["LakeDepth"].values[iorder]
                    finalriv_info.loc[
                        finalriv_info["SubId"] == downsubid, 'LakeArea'] = i_seg_info["LakeArea"].values[iorder]
                    finalriv_info.loc[
                        finalriv_info["SubId"] == downsubid, 'Laketype'] = i_seg_info["Laketype"].values[iorder]

            if iorder == len(i_seg_info) - 1:
                con_lake_up = True
                con_lake_down = False
                is_lake = True

            # if tsubid == 4703912:
            #     print(con_lake_up,con_lake_down,is_lake,iorder,len(i_seg_info) - 2,i_seg_info["HyLakeId"].values[iorder] > 0)
            #     asdf

            if iorder == len(i_seg_info) - 1:
                sum_area = 0
                seg_sub_ids = np.asarray(modifysubids)
                # if needs to add lake sub around the main stream
                if iorder_Lakeid > 0:
                    subid_cur_lake_info = finalriv_info.loc[
                        finalriv_info["HyLakeId"] == iorder_Lakeid
                    ].copy()
                    routing_info_lake = (
                        subid_cur_lake_info[["SubId", "DowSubId"]]
                        .astype("float")
                        .values
                    )
                    UpstreamLakeids = defcat(routing_info_lake, tsubid)
                    seg_sub_ids = np.unique(
                        np.concatenate([seg_sub_ids, UpstreamLakeids])
                    )
                    seg_sub_ids = seg_sub_ids[seg_sub_ids > 0]

                    # merge all subbasin not connected to the main river but drainarge to this tsubid
                All_up_subids = defcat(routing_info, tsubid)
                All_up_subids = All_up_subids[All_up_subids > 0]
                mask1 = np.in1d(
                    All_up_subids, Subid_main
                )  # exluced ids that belongs to main river stream
                All_up_subids_no_main = All_up_subids[np.logical_not(mask1)]

                seg_sub_ids = np.unique(
                    np.concatenate([seg_sub_ids, All_up_subids_no_main])
                )
                seg_sub_ids = seg_sub_ids[seg_sub_ids > 0]
                mask = np.in1d(seg_sub_ids, processed_subid)
                seg_sub_ids = seg_sub_ids[np.logical_not(mask)]

                mask_old_nonLake = np.in1d(seg_sub_ids, Old_Non_Connect_SubIds)
                seg_sub_ids = seg_sub_ids[np.logical_not(mask_old_nonLake)]

                mapoldnew_info = New_SubId_To_Dissolve(
                    subid=tsubid,
                    catchmentinfo=finalriv_info,
                    mapoldnew_info=mapoldnew_info,
                    ismodifids=1,
                    modifiidin=seg_sub_ids,
                    mainriv=Selected_riv,
                    Lake_Cat=3,
                    seg_order=seg_order,
                )
                #                    New_NonConn_Lakes = ConnectLake_to_NonConnectLake_Updateinfo(NonC_Lakeinfo = New_NonConn_Lakes,finalriv_info = finalriv_info ,Merged_subids = seg_sub_ids,Connect_Lake_ply_info = Conn_Lakes_ply,ConLakeId = iorder_Lakeid)
                modifysubids = []
                seg_order = seg_order + 1

            elif con_lake_up and con_lake_down and is_lake and con_gauge and con_area_lake:
                continue

            elif con_lake and con_gauge and con_area:
                continue

            else:
                sum_area = 0
                seg_sub_ids = np.asarray(modifysubids)
                # if needs to add lake sub around the main stream
                if iorder_Lakeid > 0:
                    subid_cur_lake_info = finalriv_info.loc[
                        finalriv_info["HyLakeId"] == iorder_Lakeid
                    ].copy()
                    routing_info_lake = (
                        subid_cur_lake_info[["SubId", "DowSubId"]]
                        .astype("float")
                        .values
                    )
                    UpstreamLakeids = defcat(routing_info_lake, tsubid)
                    seg_sub_ids = np.unique(
                        np.concatenate([seg_sub_ids, UpstreamLakeids])
                    )
                    seg_sub_ids = seg_sub_ids[seg_sub_ids > 0]

                All_up_subids = defcat(routing_info, tsubid)
                All_up_subids = All_up_subids[All_up_subids > 0]
                mask1 = np.in1d(
                    All_up_subids, Subid_main
                )  # exluced ids that belongs to main river stream
                All_up_subids_no_main = All_up_subids[np.logical_not(mask1)]

                seg_sub_ids = np.unique(
                    np.concatenate([seg_sub_ids, All_up_subids_no_main])
                )
                seg_sub_ids = seg_sub_ids[seg_sub_ids > 0]
                mask = np.in1d(seg_sub_ids, processed_subid)
                seg_sub_ids = seg_sub_ids[np.logical_not(mask)]

                mask_old_nonLake = np.in1d(seg_sub_ids, Old_Non_Connect_SubIds)
                seg_sub_ids = seg_sub_ids[np.logical_not(mask_old_nonLake)]

                mapoldnew_info = New_SubId_To_Dissolve(
                    subid=tsubid,
                    catchmentinfo=finalriv_info,
                    mapoldnew_info=mapoldnew_info,
                    ismodifids=1,
                    modifiidin=seg_sub_ids,
                    mainriv=Selected_riv,
                    Lake_Cat=3,
                    seg_order=seg_order,
                )
                modifysubids = []
                seg_order = seg_order + 1

    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0, 'HyLakeId'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0, 'Lake_Cat'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0, 'LakeVol'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0, 'LakeDepth'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0, 'LakeArea'] = 0
    mapoldnew_info.loc[mapoldnew_info['HyLakeId'] <= 0, 'Laketype'] = 0

    return (
        mapoldnew_info,
        Selected_riv_ids,
        Connected_Lake_Mainriv,
        Old_Non_Connect_LakeIds,
        Conn_To_NonConlakeids,
    )


def update_the_selected_river_to_connect_upsub_with_largest_da(Selected_riv, mapoldnew_info):
    Selected_riv_ids = np.unique(Selected_riv['SubId'].values)

    # target to find all subid has do not have subbasins
    # in the selected river network.
    lakeid_in_new_network = np.unique(Selected_riv['HyLakeId'].values)
    lakeid_in_new_network = lakeid_in_new_network[lakeid_in_new_network > 0]

    mask_lakes = np.logical_and(
        mapoldnew_info['Lake_Cat'] == 1, ~mapoldnew_info['HyLakeId'].isin(lakeid_in_new_network))

    mask_pois = mapoldnew_info['Has_POI'] > 0

    potential_sub_to_extend = mapoldnew_info[np.logical_or(
        mask_lakes, mask_pois)].copy(deep=True)

    potential_sub_to_extend = potential_sub_to_extend[~potential_sub_to_extend['SubId'].isin(
        Selected_riv_ids)]
    # extend each headwater stream
    subid_with_upstream = np.unique(potential_sub_to_extend['SubId'].values)
    for tsubid in subid_with_upstream:

        # find all upstream subbasins
        upstream_sub = mapoldnew_info[mapoldnew_info['SubId'] == tsubid].copy(
            deep=True)

        # check if one of the upstream subbasin already have river network
        upstream_sub_withriver = upstream_sub[upstream_sub['SubId'].isin(
            Selected_riv_ids)]

        # if one of the upstream subbasin already have river network, skip
        if len(upstream_sub_withriver) > 0 or len(upstream_sub) == 0:
            continue

        # remove channels all channels in this sub
        Selected_riv_ids = np.append(Selected_riv_ids, tsubid)
        # create new channel for this sub
        # sort upstream subbasins by DrainArea
        # the new channel start from a upstream subbasin with largest drainage area
        upstream_sub = upstream_sub.sort_values(
            by='DrainArea', ascending=False)

        cur_subid = upstream_sub['SubId'].values[0]
        # if tsubid == 17313:
        #     print(mapoldnew_info.columns)
        #     print(Selected_riv_ids)
        #     print(mapoldnew_info[mapoldnew_info['SubId'] == tsubid][['Old_SubId','Old_DowSubId','SubId','DowSubId','nsubid', 'ndownsubid']])
        # get the first channel which is the downsubid if the
        cur_downsubid = mapoldnew_info[mapoldnew_info['SubId']
                                       == cur_subid]['DowSubId'].values[0]
        # print(tsubid,cur_downsubid,cur_downsubid not in Selected_riv_ids)
        while cur_downsubid not in Selected_riv_ids and cur_downsubid != -1:
            # print(tsubid,cur_downsubid)
            Selected_riv_ids = np.append(Selected_riv_ids, cur_downsubid)
            cur_subinfo = mapoldnew_info[mapoldnew_info['SubId'] == cur_downsubid].copy(
                deep=True)
            if len(cur_subinfo) <= 0:
                print("check this subid ", cur_downsubid)
                break
            cur_downsubid = cur_subinfo['DowSubId'].values[0]
    Selected_riv_out = mapoldnew_info[mapoldnew_info['SubId'].isin(
        Selected_riv_ids)]
    return Selected_riv_out


def Add_River_Segment_Between_Lakes_And_Observations(mapoldnew_info, Selected_riv_ids, finalriv_infoply):
    # function to add river network between lake or poi subbasin to non lake/poi subbasins
    # used attributes in mapoldnew_info
    # nsubid the new sub id for each subbasin
    # ndownsubid the new dowsubid for each subbasin
    # Selected_riv_ids is the river segment already included in the network
    # mapoldnew_info["Old_SubId"] = mapoldnew_info["SubId"].values
    # mapoldnew_info["Old_DowSubId"] = mapoldnew_info["DowSubId"].values
    # mapoldnew_info["SubId"] = mapoldnew_info["nsubid"].values

    # mapoldnew_info["DowSubId"] = mapoldnew_info["ndownsubid"].values

    # target to find all subid has upstream subbasins
    # and none upstream subbains have river network
    # get all subbasins have a upstream subbasin
    mask = mapoldnew_info['SubId'].isin(mapoldnew_info['DowSubId'].values)
    sub_with_upstream = mapoldnew_info[mask].copy(deep=True)

    subid_with_upstream = np.unique(sub_with_upstream['SubId'].values)
    for tsubid in subid_with_upstream:

        # find all upstream subbasins
        upstream_sub = mapoldnew_info[mapoldnew_info['DowSubId'] == tsubid].copy(
            deep=True)

        # check if one of the upstream subbasin already have river network
        upstream_sub_withriver = upstream_sub[upstream_sub['SubId'].isin(
            Selected_riv_ids)]

        # if one of the upstream subbasin already have river network, skip
        if len(upstream_sub_withriver) > 0:
            continue

        # remove channels all channels in this sub
        mask = np.isin(
            Selected_riv_ids, mapoldnew_info[mapoldnew_info['SubId'] == tsubid]['Old_SubId'].values)
        Selected_riv_ids = Selected_riv_ids[~mask]

        # create new channel for this sub
        # sort upstream subbasins by DrainArea
        # the new channel start from a upstream subbasin with largest drainage area
        upstream_sub = upstream_sub.sort_values(
            by='DrainArea', ascending=False)
        cur_subid = upstream_sub['SubId'].values[0]
        # if tsubid == 17313:
        #     print(mapoldnew_info.columns)
        #     print(Selected_riv_ids)
        #     print(mapoldnew_info[mapoldnew_info['SubId'] == tsubid][['Old_SubId','Old_DowSubId','SubId','DowSubId','nsubid', 'ndownsubid']])
        # get the first channel which is the downsubid if the
        cur_downsubid = finalriv_infoply[finalriv_infoply['SubId']
                                         == cur_subid]['DowSubId'].values[0]
        # print(tsubid,cur_downsubid,cur_downsubid not in Selected_riv_ids)
        while cur_downsubid not in Selected_riv_ids and cur_downsubid != -1:
            # print(tsubid,cur_downsubid)
            Selected_riv_ids = np.append(Selected_riv_ids, cur_downsubid)
            cur_subinfo = finalriv_infoply[finalriv_infoply['SubId'] == cur_downsubid].copy(
                deep=True)
            if len(cur_subinfo) <= 0:
                print("check this subid ", cur_downsubid)
                break
            cur_downsubid = cur_subinfo['DowSubId'].values[0]
    # print(17313 in Selected_riv_ids)
    # remove all subid that are already included in the river network
    mapoldnew_info = update_channel_attributes(
        Selected_riv_ids, mapoldnew_info, finalriv_infoply)
    return Selected_riv_ids, mapoldnew_info


def update_channel_attributes(Selected_riv_ids, mapoldnew_info, finalriv_infoply):

    # obtain subid in the new network that needs to be updated

    subid_update_riv = np.unique(
        mapoldnew_info[mapoldnew_info['Old_SubId'].isin(Selected_riv_ids)]['SubId'].values)
    all_river_segments = finalriv_infoply[finalriv_infoply['SubId'].isin(
        Selected_riv_ids)]
    for tsubid in subid_update_riv:
        mask_new = mapoldnew_info['SubId'] == tsubid
        old_subids = np.unique(mapoldnew_info[mask_new]['Old_SubId'].values)
        riv_in_new_sub = all_river_segments[all_river_segments['SubId'].isin(
            old_subids)]

        if len(riv_in_new_sub) > 0:
            mapoldnew_info.loc[mask_new, "RivLength"] = np.sum(
                riv_in_new_sub["RivLength"].values)
            mapoldnew_info.loc[mask_new, "RivSlope"] = np.average(
                riv_in_new_sub["RivSlope"].values,
                weights=riv_in_new_sub["RivLength"].values,
            )
            mapoldnew_info.loc[mask_new, "FloodP_n"] = np.average(
                riv_in_new_sub["FloodP_n"].values,
                weights=riv_in_new_sub["RivLength"].values,
            )
            mapoldnew_info.loc[mask_new, "Q_Mean"] = np.average(
                riv_in_new_sub["Q_Mean"].values,
                weights=riv_in_new_sub["RivLength"].values,
            )
            mapoldnew_info.loc[mask_new, "Ch_n"] = np.average(
                riv_in_new_sub["Ch_n"].values,
                weights=riv_in_new_sub["RivLength"].values,
            )
            mapoldnew_info.loc[mask_new, "BkfWidth"] = np.max(
                riv_in_new_sub["BkfWidth"].values)
            mapoldnew_info.loc[mask_new, "BkfDepth"] = np.max(
                riv_in_new_sub["BkfDepth"].values)

            mapoldnew_info.loc[mask_new, "Max_DEM"] = np.max(
                riv_in_new_sub["Max_DEM"].values)
            mapoldnew_info.loc[mask_new, "Min_DEM"] = np.min(
                riv_in_new_sub["Min_DEM"].values)
    return mapoldnew_info


def Return_Selected_Lakes_Attribute_Table_And_Id(
    finalcat_info,
    Thres_Area_Conn_Lakes,
    Thres_Area_Non_Conn_Lakes,
    Selected_Lake_List_in,
):
    """Retrun selected lake's attribute table and ID
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    finalcat_info["LakeArea"] = finalcat_info["LakeArea"].astype(float)
    finalcat_info["HyLakeId"] = finalcat_info["HyLakeId"].astype(int)
    finalcat_info["Seg_ID"] = finalcat_info["Seg_ID"].astype(int)

    Non_ConnL_info = finalcat_info.loc[finalcat_info["Lake_Cat"] == 2].copy(
        deep=True)
    ConnL_info = finalcat_info.loc[finalcat_info["Lake_Cat"] == 1].copy(
        deep=True)

    if "Has_POI" in finalcat_info.columns:
        Gauge_col_Name = "Has_POI"
    else:
        Gauge_col_Name = "Has_Gauge"
    # first obtain selected lakes based on lake area
    # process connected lakes first
    mask1 = ConnL_info["LakeArea"] >= Thres_Area_Conn_Lakes
    mask2 = ConnL_info[Gauge_col_Name] == 1
    mask = np.logical_or(mask1, mask2)
    Selected_ConnLakes = ConnL_info.loc[
        mask
    ]["HyLakeId"].values

    Selected_ConnLakes = Selected_ConnLakes[Selected_ConnLakes > 0]
    Selected_ConnLakes = np.unique(Selected_ConnLakes)

    mask1 = Non_ConnL_info["LakeArea"] >= Thres_Area_Non_Conn_Lakes
    mask2 = Non_ConnL_info[Gauge_col_Name] == 1
    mask = np.logical_or(mask1, mask2)
    # process non connected selected lakes
    Selected_Non_ConnLakes = Non_ConnL_info[
        mask
    ]["HyLakeId"].values
    Selected_Non_ConnLakes = Selected_Non_ConnLakes[Selected_Non_ConnLakes > 0]
    Selected_Non_ConnLakes = np.unique(Selected_Non_ConnLakes)

    # selecte lake by list
    if len(Selected_Lake_List_in) > 0:
        Selected_Lake_List = Selected_Lake_List_in
    else:
        # invalide lake id, no lake will be selected by list
        Selected_Lake_List = [-9999999]
    All_ConnL = ConnL_info["HyLakeId"].values
    All_Non_ConnL = Non_ConnL_info["HyLakeId"].values
    Selected_Lake_List_in_array = np.array(Selected_Lake_List)

    mask_CL = np.in1d(All_ConnL, Selected_Lake_List_in_array)
    mask_NCL = np.in1d(All_Non_ConnL, Selected_Lake_List_in_array)

    Selected_ConnLakes_List = All_ConnL[mask_CL]
    Selected_ConnLakes_List = np.unique(Selected_ConnLakes_List)
    Selected_Non_ConnLakes_List = All_Non_ConnL[mask_NCL]
    Selected_Non_ConnLakes_List = np.unique(Selected_Non_ConnLakes_List)

    # combine lake from list and area
    Selected_Non_ConnLakes = np.concatenate(
        (Selected_Non_ConnLakes, Selected_Non_ConnLakes_List), axis=0)
    Selected_ConnLakes = np.concatenate(
        (Selected_ConnLakes, Selected_ConnLakes_List), axis=0)

    Un_Selected_ConnLakes_info = finalcat_info.loc[
        (finalcat_info["Lake_Cat"] == 1)
        & (np.logical_not(finalcat_info["HyLakeId"].isin(Selected_ConnLakes)))
    ]
    Un_Selected_Non_ConnL_info = finalcat_info.loc[
        (finalcat_info["Lake_Cat"] == 2)
        & (np.logical_not(finalcat_info["HyLakeId"].isin(Selected_Non_ConnLakes)))
    ]

    return (
        Selected_Non_ConnLakes,
        Selected_ConnLakes,
        Un_Selected_ConnLakes_info,
        Un_Selected_Non_ConnL_info,
    )


def Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Remove_CL(
    finalcat_info_temp, Un_Selected_ConnLakes_info
):
    """Change attributes for catchments that needs to be merged due to remove
        connected lakes
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    mapoldnew_info = finalcat_info_temp.copy(deep=True)
    mapoldnew_info["nsubid"] = -1
    mapoldnew_info["ndownsubid"] = -1
    mapoldnew_info["nsubid2"] = mapoldnew_info["SubId"]
    mapoldnew_info.reset_index(drop=True, inplace=True)
    # Loop each unselected lake stream seg

    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in mapoldnew_info.columns:
        Gauge_col_Name = "Has_Gauge"

    Seg_IDS = Un_Selected_ConnLakes_info["Seg_ID"].values
    Seg_IDS = np.unique(Seg_IDS)
    for iseg in range(0, len(Seg_IDS)):
        #            print('#########################################################################################33333')
        i_seg_id = Seg_IDS[iseg]
        i_seg_info = finalcat_info_temp.loc[
            (finalcat_info_temp["Seg_ID"] == i_seg_id)
            & (finalcat_info_temp["Lake_Cat"] != 2)
        ].copy()
        i_seg_info = i_seg_info.sort_values(["Seg_order"], ascending=(True))

        # each part of the segment are not avaiable to be merged
        N_Hylakeid = np.unique(i_seg_info["HyLakeId"].values)
        if len(i_seg_info) == len(N_Hylakeid):
            continue

        # All lakes in this segment are removed
        if np.max(N_Hylakeid) and np.max(i_seg_info[Gauge_col_Name].values) <= 0:
            tsubid = i_seg_info["SubId"].values[len(i_seg_info) - 1]
            seg_sub_ids = i_seg_info["SubId"].values
            mapoldnew_info = New_SubId_To_Dissolve(
                subid=tsubid,
                catchmentinfo=finalcat_info_temp,
                mapoldnew_info=mapoldnew_info,
                ismodifids=1,
                modifiidin=seg_sub_ids,
                mainriv=finalcat_info_temp,
                Lake_Cat=-1,
                seg_order=1,
            )

        # loop from the first order of the current segment
        modifysubids = []
        seg_order = 1
        for iorder in range(0, len(i_seg_info)):
            tsubid = i_seg_info["SubId"].values[iorder]
            modifysubids.append(tsubid)

            # two seg has the same HyLakeId id, can be merged
            if iorder == len(i_seg_info) - 1:
                seg_sub_ids = np.asarray(modifysubids)
                mapoldnew_info = New_SubId_To_Dissolve(
                    subid=tsubid,
                    catchmentinfo=finalcat_info_temp,
                    mapoldnew_info=mapoldnew_info,
                    ismodifids=1,
                    modifiidin=seg_sub_ids,
                    mainriv=finalcat_info_temp,
                    Lake_Cat=i_seg_info["HyLakeId"].values[iorder],
                    seg_order=seg_order,
                )
                modifysubids = []
                seg_order = seg_order + 1

            elif (
                i_seg_info["HyLakeId"].values[iorder]
                == i_seg_info["HyLakeId"].values[iorder + 1]
                and i_seg_info[Gauge_col_Name].values[iorder] <= 0
            ):
                continue
            else:
                seg_sub_ids = np.asarray(modifysubids)
                mapoldnew_info = New_SubId_To_Dissolve(
                    subid=tsubid,
                    catchmentinfo=finalcat_info_temp,
                    mapoldnew_info=mapoldnew_info,
                    ismodifids=1,
                    modifiidin=seg_sub_ids,
                    mainriv=finalcat_info_temp,
                    Lake_Cat=i_seg_info["HyLakeId"].values[iorder],
                    seg_order=seg_order,
                )
                modifysubids = []
                seg_order = seg_order + 1
    return mapoldnew_info


def Change_Attribute_Values_For_Catchments_Need_To_Be_Merged_By_Remove_NCL(
    mapoldnew_info, finalcat_info_temp, Un_Selected_Non_ConnL_info
):
    """Change attributes for catchments that needs to be merged due to remove
        non connected lakes
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """

    Un_Selected_Non_ConnL_info = Un_Selected_Non_ConnL_info.sort_values(
        ["DrainArea"], ascending=[False]
    )
    for i in range(0, len(Un_Selected_Non_ConnL_info)):
        Remove_Non_ConnL_Lakeid = Un_Selected_Non_ConnL_info["HyLakeId"].values[
            i
        ]  # select one non connected lake
        Remove_Non_ConnL_Lake_Sub_info = finalcat_info_temp.loc[
            finalcat_info_temp["HyLakeId"] == Remove_Non_ConnL_Lakeid
        ].copy()  # obtain cat info of this non connected lake catchment
        if len(Remove_Non_ConnL_Lake_Sub_info) != 1:
            print("It is not a non connected lake catchment")
            print(Remove_Non_ConnL_Lake_Sub_info)
            continue
        modifysubids = []  # array store all catchment id needs to be merged
        # intial subid
        csubid = Remove_Non_ConnL_Lake_Sub_info["SubId"].values[0]
        downsubid = Remove_Non_ConnL_Lake_Sub_info["DowSubId"].values[
            0
        ]  # downstream subid
        modifysubids.append(csubid)  # add inital subid into the list
        tsubid = -1
        is_pre_modified = 0

        Down_Sub_info = finalcat_info_temp.loc[
            finalcat_info_temp["SubId"] == downsubid
        ].copy()  # obtain downstream infomation
        if len(Down_Sub_info) == 0 and downsubid < 0:
            continue
        if Down_Sub_info["Lake_Cat"].values[0] != 2:
            # check if this downsubid has a new subid
            nsubid = mapoldnew_info.loc[mapoldnew_info["SubId"] == downsubid][
                "nsubid"
            ].values[
                0
            ]  # check if this catchment has modified: either merged to other catchment

            if nsubid > 0:  # if it been already processed
                tsubid = nsubid  # the target catchment id will be the newsunid
                is_pre_modified = 1  # set is modifed to 1
                modifysubids.append(tsubid)  # add this subid to the list
            else:
                tsubid = downsubid
                is_pre_modified = -1
                modifysubids.append(tsubid)

        else:  # down stream cat is a non connected lake cat

            nsubid = mapoldnew_info.loc[mapoldnew_info["SubId"] == downsubid][
                "nsubid"
            ].values[0]
            if nsubid > 0:
                tsubid = nsubid
                is_pre_modified = 1
                modifysubids.append(tsubid)
            else:
                tsubid = downsubid
                #                    print("This value should always be zero:    ", np.sum(Un_Selected_Non_ConnL_info['SubId'] == downsubid))
                is_pre_modified = -1
                modifysubids.append(tsubid)

        Tar_Lake_Id = mapoldnew_info[mapoldnew_info["SubId"] == tsubid][
            "HyLakeId"
        ].values[0]

        if is_pre_modified > 0:
            mapoldnew_info = New_SubId_To_Dissolve(
                subid=tsubid,
                catchmentinfo=mapoldnew_info,
                mapoldnew_info=mapoldnew_info,
                ismodifids=1,
                modifiidin=modifysubids,
                mainriv=mapoldnew_info,
                Lake_Cat=Tar_Lake_Id,
            )
        else:
            mapoldnew_info = New_SubId_To_Dissolve(
                subid=tsubid,
                catchmentinfo=mapoldnew_info,
                mapoldnew_info=mapoldnew_info,
                ismodifids=1,
                modifiidin=modifysubids,
                mainriv=finalcat_info_temp,
                Lake_Cat=Tar_Lake_Id,
            )

    unprocessedmask = mapoldnew_info["nsubid"] <= 0

    mapoldnew_info.loc[unprocessedmask, "nsubid"] = mapoldnew_info.loc[
        unprocessedmask, "nsubid2"
    ].values

    mapoldnew_info.drop(columns=["nsubid2"])

    return mapoldnew_info


def Change_Attribute_Values_For_Catchments_Covered_By_Same_Lake(finalrivply_info):
    """Change attributes for catchments that covered by the same lake.
    ----------

    Notes
    -------
        For example, lake 'la' covering catchment a,b,c. the lake outlet catchment
        is a. then this function will change attribute of b and c to a.
    Returns:
    -------
        None,
    """

    sub_colnm = "SubId"
    mapoldnew_info = finalrivply_info.copy(deep=True)
    mapoldnew_info["nsubid"] = mapoldnew_info["SubId"].values
    mapoldnew_info['nsubid'] = mapoldnew_info['nsubid'].astype('int32')
    AllConnectLakeIDS = finalrivply_info["HyLakeId"].values
    AllConnectLakeIDS = AllConnectLakeIDS[AllConnectLakeIDS > 0]
    AllConnectLakeIDS = np.unique(AllConnectLakeIDS)

    # process connected lakes  merge polygons
    for i in range(0, len(AllConnectLakeIDS)):
        lakeid = AllConnectLakeIDS[i]
        Lakesub_info = finalrivply_info.loc[finalrivply_info["HyLakeId"] == lakeid]
        Lakesub_info = Lakesub_info.sort_values(
            ["DrainArea"], ascending=(False))
        tsubid = Lakesub_info[sub_colnm].values[
            0
        ]  # outlet subbasin id with highest acc
        lakesubids = Lakesub_info[sub_colnm].values
        if len(lakesubids) > 1:  # only for connected lakes
            mapoldnew_info = New_SubId_To_Dissolve(
                subid=tsubid,
                catchmentinfo=finalrivply_info,
                mapoldnew_info=mapoldnew_info,
                ismodifids=1,
                modifiidin=lakesubids,
                mainriv=finalrivply_info,
                Lake_Cat=1,
            )
    return mapoldnew_info


def Return_SubIds_Between_Two_Subbasins_In_Rouing_Network(
    routing_info, subid_downstream, subid_upstream
):
    """Return subid of subbasins between two subid in the routing network
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """

    HydroBasins1 = defcat(
        routing_info, subid_downstream
    )  # return fid of polygons that needs to be select

    if subid_upstream > 0:
        HydroBasins2 = defcat(routing_info, subid_upstream)
        for i in range(len(HydroBasins2)):
            if HydroBasins2[i] == subid_upstream:
                continue
            rows = np.argwhere(HydroBasins1 == HydroBasins2[i])
            HydroBasins1 = np.delete(HydroBasins1, rows)
        HydroBasins = HydroBasins1
    else:
        HydroBasins = HydroBasins1

    return HydroBasins


def return_non_lake_inflow_segs_and_segs_within_lakes(
    riv_lake_id, str_id, riv_lake_id2, cl_id, routing_info, str_start_pt_lakeid
):

    routing_info_rout = routing_info[[
        "SubId", "DowSubId"]].astype("int").values
    # combine river lake id and str id and lake id
    riv_lake = np.column_stack((riv_lake_id, str_id))
    riv_lake = riv_lake[riv_lake[:, 0].argsort()]
    riv_lake_cl = np.column_stack((riv_lake_id2, cl_id))
    riv_lake_cl = riv_lake_cl[riv_lake_cl[:, 0].argsort()]
    riv_lake = np.append(riv_lake, riv_lake_cl, axis=1)

    # assign acc for each row
    str_id_unique = np.unique(np.array(riv_lake[:, 1]))
    cl_id_unique = np.unique(np.array(riv_lake[:, 3]))

    acc_str_cl = np.full((len(riv_lake), 1), np.nan)
    for i in range(0, len(str_id_unique)):
        strid = str_id_unique[i]
        # print(strid)
        # print(routing_info.loc[routing_info["SubId"] == strid])
        maxacc = routing_info.loc[routing_info["SubId"] == strid]["MaxAcc_cat"].values[
            0
        ]
        acc_str_cl[riv_lake[:, 1] == strid] = maxacc
    riv_lake = np.append(riv_lake, acc_str_cl, axis=1)

    # loop for each lake, identify, outlet seg, and inlet segs
    # and segs between outlet and inlet segs

    non_lake_inflow_segs = []
    str_id_lake_inlfow = []
    str_id_within_lakes = []
    for i in range(0, len(cl_id_unique)):
        # obtain current lakeid
        cl_id = cl_id_unique[i]
        # get all new riv seg id and str id covered by this lake
        riv_lake_il = riv_lake[riv_lake[:, 3] == cl_id]
        riv_lake_il = riv_lake_il[riv_lake_il[:, 4].argsort()]
        outlet_str = riv_lake_il[len(riv_lake_il) - 1, 1]
        sub_routing_info = routing_info.loc[
            routing_info["SubId"].isin(riv_lake_il[:, 1])
        ]
        # if lake only overlay with one str, skip
        # this river is lake infow str and lake outflow str
        if len(riv_lake_il) == 1:
            str_id_cl_j = riv_lake_il[0, 1]
            # the river can not start within the lake
            if (
                str_start_pt_lakeid.loc[str_start_pt_lakeid["IL_SubId"] == str_id_cl_j][
                    "sl_cl_id"
                ].values[0]
                != cl_id
            ):
                str_id_lake_inlfow.append(int(riv_lake_il[0, 1]))
            else:
                non_lake_inflow_segs.append(int(riv_lake_il[0, 0]))
            continue

        strs_to_lake_outlet = defcat(routing_info_rout, outlet_str)

        # excpet str of the lake outlet, all str covered by the lake
        # will be removed
        # at the same time, all dwon str if it of each str coverd by the lakes
        # when they are flow to the lake outlet
        # it will aslo be removed
        unique_str_id_cv_by_lake = np.unique(riv_lake_il[:, 1])

        for j in range(0, len(riv_lake_il)):
            str_id_cl_j = riv_lake_il[j, 1]
            # check if there is any upstream river id in the str covered by the lake and there is some upstream
            # and the begining of the str is not located in current lake
            if str_id_cl_j != outlet_str:
                str_id_within_lakes.append(int(riv_lake_il[j, 1]))
                # check down stream id of the currunt str covered by the lake
                if len(routing_info.loc[routing_info['SubId'] == str_id_cl_j]) > 0:
                    # get down stream of the corrent str covered by the lake
                    down_sun_of_lake_cover_str = routing_info.loc[routing_info['SubId']
                                                                  == str_id_cl_j]['DowSubId'].values[0]
                    # if the it is also flow to the lake outlet , remove it
                    # and it is not already in the list
                    if down_sun_of_lake_cover_str in strs_to_lake_outlet and down_sun_of_lake_cover_str not in riv_lake_il[:, 1]:
                        # remove all dnow str unitl it reach to a lake coverd str_r
                        # or a str already be removed
                        str_id_within_lakes.append(
                            int(down_sun_of_lake_cover_str))
                        #
                        cstrid = down_sun_of_lake_cover_str
                        k_d = 0
                        while k_d < 100:
                            k_d = k_d + 1
                            if len(routing_info.loc[routing_info['SubId'] == cstrid]) > 0:
                                dstr_id = routing_info.loc[routing_info['SubId']
                                                           == cstrid]['DowSubId'].values[0]
                                # check if it reach a lake covered str id
                                if dstr_id not in riv_lake_il[:, 1] and dstr_id in strs_to_lake_outlet:
                                    str_id_within_lakes.append(int(dstr_id))
                                    cstrid = dstr_id
                                else:
                                    break
                            else:
                                break

            # obtain all upstream of current lake str
            up_str_cl_j = defcat(routing_info_rout, str_id_cl_j)

            if len(up_str_cl_j) > 0:
                # check if there is more than 1 up stream is
                # covreed by the lake
                #
                mask = np.isin(up_str_cl_j, unique_str_id_cv_by_lake)
                if sum(mask) >= 2:
                    has_up_str_inlake = True
                else:
                    has_up_str_inlake = False
            else:
                has_up_str_inlake = False

            if (
                len(sub_routing_info.loc[sub_routing_info["DowSubId"] == str_id_cl_j])
                > 0
                or str_start_pt_lakeid.loc[
                    str_start_pt_lakeid["IL_SubId"] == str_id_cl_j
                ]["sl_cl_id"].values[0]
                == cl_id
                or has_up_str_inlake
            ):
                # there is some str drainge to this str,
                # so this lake-riv seg is not the lake inflow segment
                non_lake_inflow_segs.append(int(riv_lake_il[j, 0]))
            else:
                str_id_lake_inlfow.append(int(riv_lake_il[j, 1]))

    return list(set(str_id_within_lakes)), list(set(non_lake_inflow_segs)), list(set(str_id_lake_inlfow))


def Check_If_Lake_Have_Multi_OutLet(CL_Id, Str_Id, Routing_info):

    # create a emppty array to store lake id with multi outlet
    Lakes_WIth_Multi_Outlet = []
    Remove_Str = []
    # combine lake id and coresponding str id into two dim array
    CL_Str = np.column_stack((CL_Id, Str_Id))
    # obtain unique str id
    Str_Id_unique = np.unique(np.array(Str_Id))
    # obtain unique lake ids
    CL_Id_unique = np.unique(np.array(CL_Id))

    # add maxacc of each str to CL_str
    Acc_Str_CL = np.full((len(CL_Str), 1), np.nan)
    for i in range(0, len(Str_Id_unique)):
        strid = Str_Id_unique[i]
        maxacc = Routing_info.loc[Routing_info["SubId"] == strid]["MaxAcc_cat"].values[
            0
        ]
        Acc_Str_CL[CL_Str[:, 1] == strid] = maxacc
    CL_Str = np.append(CL_Str, Acc_Str_CL, axis=1)

    # sort CL_str based on max acc of
    CL_Str = CL_Str[CL_Str[:, 2].argsort()]

    # routing info
    routing_info_only = Routing_info[[
        "SubId", "DowSubId"]].astype("int").values
    #    print(routing_info_only)
    # check if lakes have multi outlet
    for i in range(0, len(CL_Id_unique)):
        # got lake id
        lake_id = CL_Id_unique[i]
        # got all str overlaied by this lake id
        i_CL_Str = CL_Str[CL_Str[:, 0] == lake_id]

        # lake is overlaied by one str, so must only have one lake outlet
        if len(i_CL_Str) <= 1:
            continue

        # len(i_CL_Str)>1
        # check if all str covered by this lake will drainage to the str with
        # maximum acc. if it is ture than lake has only one outelt,
        # if not the lake have multi outlet

        # obtain str id with max acc among strs covered by the lake
        str_max_acc = i_CL_Str[len(i_CL_Str) - 1, 1]
        #        print(str_max_acc,len(routing_info_only),"#################")
        # obtain all upstream str id of max acc str
        str_to_str_max_acc = defcat(routing_info_only, str_max_acc)

        #        if len(str_to_str_max_acc) <= 1:
        #            str_to_str_max_acc = defcat(
        #                routing_info_only, i_CL_Str[len(i_CL_Str) - 2, 1]
        #            )
        # create a mask for i_CL_Str[:,1], it will be true for in positon
        # where it's value in str_to_str_max_acc
        mask = np.isin(i_CL_Str[:, 1], str_to_str_max_acc)
        # not all element in i_CL_Str[:,1] exist in str_to_str_max_acc
        # the lake have muli outlet

        Remove_Str_i = []

        # if not all str drainage to outlet str
        if len(mask[mask == True]) < len(i_CL_Str[:, 1]):
            # obtain strids of str not drainage to outlet str
            str_notflowto_lakeoutlet = i_CL_Str[np.logical_not(mask), 1]
            # loop for strids of str not drainage to outlet str
            for istr in range(0, len(str_notflowto_lakeoutlet)):
                # get i str id
                strid = str_notflowto_lakeoutlet[istr]
                # check  streams drainage to current str
                #                upstrs  = Defcat(routing_info_only,strid)
                # no upstream str, ### remove str instead remove lake
                #                if len(upstrs) == 1:
                if strid != str_max_acc:
                    Remove_Str_i.append(strid)

            # if len(Remove_Str_i) == len(str_notflowto_lakeoutlet):
            #     Remove_Str = Remove_Str + Remove_Str_i
            # else:
            Remove_Str = Remove_Str + Remove_Str_i
            Lakes_WIth_Multi_Outlet.append(lake_id)

    return Lakes_WIth_Multi_Outlet, Remove_Str


def change_attribute_values_for_catchments_covered_by_same_lake(finalrivply_info):
    """Change attributes for catchments that covered by the same lake.
    ----------

    Notes
    -------
        For example, lake 'la' covering catchment a,b,c. the lake outlet catchment
        is a. then this function will change attribute of b and c to a.
    Returns:
    -------
        None,
    """

    sub_colnm = "SubId"
    mapoldnew_info = finalrivply_info.copy(deep=True)
    mapoldnew_info["nsubid"] = mapoldnew_info["SubId"].values
    AllConnectLakeIDS = finalrivply_info["HyLakeId"].values
    AllConnectLakeIDS = AllConnectLakeIDS[AllConnectLakeIDS > 0]
    AllConnectLakeIDS = np.unique(AllConnectLakeIDS)

    # process connected lakes  merge polygons
    for i in range(0, len(AllConnectLakeIDS)):
        lakeid = AllConnectLakeIDS[i]
        Lakesub_info = finalrivply_info.loc[finalrivply_info["HyLakeId"] == lakeid]
        Lakesub_info = Lakesub_info.sort_values(
            ["DrainArea"], ascending=(False))
        tsubid = Lakesub_info[sub_colnm].values[
            0
        ]  # outlet subbasin id with highest acc
        lakesubids = Lakesub_info[sub_colnm].values
        if len(lakesubids) > 1:  # only for connected lakes
            mapoldnew_info = New_SubId_To_Dissolve(
                subid=tsubid,
                catchmentinfo=finalrivply_info,
                mapoldnew_info=mapoldnew_info,
                ismodifids=1,
                modifiidin=lakesubids,
                mainriv=finalrivply_info,
                Lake_Cat=1,
            )
    return mapoldnew_info


def update_topology(mapoldnew_info, UpdateStreamorder=1, UpdateSubId=1):
    """Functions will update subid,downsubid, calcuate stream order and
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
        for i in range(0, len(idx)):
            nsubid = mapoldnew_info.loc[idx[i], "nsubid"]
            subid = mapoldnew_info.loc[idx[i], "SubId"]
            odownsubid = mapoldnew_info.loc[idx[i], "DowSubId"]

            donsubidinfo = mapoldnew_info.loc[
                mapoldnew_info["SubId"] == odownsubid
            ].copy()

            if len(donsubidinfo) > 0:
                mapoldnew_info.loc[idx[i], "ndownsubid"] = donsubidinfo[
                    "nsubid"
                ].values[0]
            else:
                mapoldnew_info.loc[idx[i], "ndownsubid"] = -1

            if nsubid == mapoldnew_info.loc[idx[i], "ndownsubid"]:
                mapoldnew_info.loc[idx[i], "ndownsubid"] = -1

            if nsubid == mapoldnew_info.loc[idx[i], "ndownsubid"]:
                mapoldnew_info.loc[idx[i], "ndownsubid"] = -1

        mapoldnew_info["Old_SubId"] = mapoldnew_info["SubId"].values
        mapoldnew_info["Old_DowSubId"] = mapoldnew_info["DowSubId"].values
        mapoldnew_info["SubId"] = mapoldnew_info["nsubid"].values

        mapoldnew_info["DowSubId"] = mapoldnew_info["ndownsubid"].values

        riv_pd_nncls_routing_info = mapoldnew_info[mapoldnew_info['Lake_Cat'] != 2].copy(
            deep=True)
        mask = ~mapoldnew_info['SubId'].isin(
            riv_pd_nncls_routing_info['DowSubId'])

        # if "DA_Chn_L" in mapoldnew_info.columns:
        #     mapoldnew_info.loc[mask, "DA_Chn_L"] = -1.2345
        #     mapoldnew_info.loc[mask, "DA_Chn_Slp"] = -1.2345
        # mapoldnew_info.loc[mask, "RivLength"] = -1.2345
        # mapoldnew_info.loc[mask, "RivSlope"] = -1.2345
        # mapoldnew_info.loc[mask, "FloodP_n"] = -1.2345
        # mapoldnew_info.loc[mask, "Ch_n"] = -1.2345
        # mapoldnew_info.loc[mask, "Max_DEM"] = -1.2345
        # mapoldnew_info.loc[mask, "Min_DEM"] = -1.2345

    if UpdateStreamorder < 0:
        return mapoldnew_info

    mapoldnew_info_unique = mapoldnew_info.drop_duplicates(
        "SubId", keep="first")

    mapoldnew_info_unique = streamorderanddrainagearea(mapoldnew_info_unique)

    for i in range(0, len(mapoldnew_info_unique)):
        isubid = mapoldnew_info_unique["SubId"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "Strahler"
        ] = mapoldnew_info_unique["Strahler"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "Seg_ID"
        ] = mapoldnew_info_unique["Seg_ID"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "Seg_order"
        ] = mapoldnew_info_unique["Seg_order"].values[i]
        mapoldnew_info.loc[
            mapoldnew_info["SubId"] == isubid, "DrainArea"
        ] = mapoldnew_info_unique["DrainArea"].values[i]

    return mapoldnew_info


def func_Q_DA(A, k, c):
    return k * A ** c


def return_k_and_c_in_q_da_relationship(da_q):

    try:
        popt2, pcov2 = curve_fit(func_Q_DA, da_q[:, 0], da_q[:, 1])
    except RuntimeError:
        print("#######################################################")
        popt2 = np.full(2, -1)

    return popt2[0], popt2[1]


def calculateChannaln(width, depth, Q, slope):
    zch = 2
    sidwd = zch * depth  # river side width
    tab = "          "
    botwd = width - 2 * sidwd  # river
    if botwd < 0:
        botwd = 0.5 * width
        sidwd = 0.5 * 0.5 * width
        zch = (width - botwd) / 2 / depth
    Ach = botwd * depth + 2 * zch * depth * depth / 2
    #    arcpy.AddMessage(depth)
    #    arcpy.AddMessage(zch)
    #    arcpy.AddMessage(botwd)
    #    arcpy.AddMessage(width)
    #    arcpy.AddMessage(slope)

    Pch = botwd + 2 * depth * (1 + zch ** 2) ** 0.5
    Rch = float(Ach) / float(Pch)  # in meter
    V = float(Q) / float(Ach)
    if V > 0:
        n = (Rch ** (2.0 / 3.0)) * (slope ** (1.0 / 2.0)) / V
    else:
        n = -1.2345
    return n


def return_interest_catchments_info(catinfo, outlet_obs_id, path_sub_reg_outlets_v="#"):

    if outlet_obs_id < 0:
        return catinfo

    routing_info = catinfo[["SubId", "DowSubId"]].astype("float").values

    Gauge_col_Name = "Has_POI"
    if "Has_POI" not in catinfo.columns:
        Gauge_col_Name = "Has_Gauge"

    if path_sub_reg_outlets_v != "#":

        Sub_reg_outlets = Dbf_To_Dataframe(path_sub_reg_outlets_v)[
            "reg_subid"].values
        Sub_reg_outlets_ids = np.unique(Sub_reg_outlets)
        Sub_reg_outlets_ids = Sub_reg_outlets_ids[Sub_reg_outlets_ids > 0]

        # Find all obervation id that is subregion outlet
        reg_outlet_info = catinfo.loc[catinfo[Gauge_col_Name].isin(
            Sub_reg_outlets_ids)]

        # Define outlet ID
        outletid = -1

        if outlet_obs_id < 0:
            print("To use subregion, the Subregion Id MUST provided as Outlet_Obs_ID")
            return catinfo

        outletID_info = catinfo.loc[catinfo[Gauge_col_Name] == outlet_obs_id]
        if len(outletID_info) > 0:
            outletid = outletID_info["SubId"].values[0]
        else:
            print("No Outlet id is founded for subregion   ", outletID_info)
            return catinfo

        # find all subregion drainge to this outlet id
        HydroBasins1 = defcat(routing_info, outletid)
        # if there is other subregion outlet included in current sturcture
        # remove subbasins drainge to them

        if len(reg_outlet_info) >= 2:  # has upstream regin outlet s
            # remove subbains drainage to upstream regin outlet s
            for i in range(0, len(reg_outlet_info)):
                upregid = reg_outlet_info["SubId"].values[i]
                # the subregion ouetlet not within the target domain neglect
                if upregid == outletid or np.sum(np.in1d(HydroBasins1, upregid)) < 1:
                    continue
                HydroBasins_remove = defcat(routing_info, upregid)
                mask = np.in1d(
                    HydroBasins1, HydroBasins_remove
                )  # exluced ids that belongs to main river stream
                HydroBasins1 = HydroBasins1[np.logical_not(mask)]
        HydroBasins = HydroBasins1

        catinfo = catinfo.loc[catinfo["SubId"].isin(HydroBasins)]
        return catinfo
    # selected based on observation guage obs id
    else:
        outletid = -1
        outletID_info = catinfo.loc[catinfo[Gauge_col_Name] == outlet_obs_id]
        if len(outletID_info) > 0:
            outletid = outletID_info["SubId"].values[0]
            # find upsteam catchment id
            HydroBasins = defcat(routing_info, outletid)
        else:
            HydroBasins = catinfo["SubId"].values

        catinfo = catinfo.loc[catinfo["SubId"].isin(HydroBasins)]
        return catinfo


def create_grid_weight_main(Mapforcing, Forcinfo):

    grid_weight_string_list = []
    hruids = Mapforcing["HRU_ID"].values
    hruids = np.unique(hruids)
    #    Lakeids = np.unique(Lakeids)
    grid_weight_string_list = []

    os.environ["JOBLIB_TEMP_FOLDER"] = tempfile.gettempdir()

    grid_weight_string_list.append(":GridWeights")
    grid_weight_string_list.append("   #      ")
    grid_weight_string_list.append("   # [# HRUs]")

    sNhru = len(hruids)
    grid_weight_string_list.append("   :NumberHRUs       " + str(sNhru))

    sNcell = (max(Forcinfo["Row"].values) + 1) * \
        (max(Forcinfo["Col"].values) + 1)
    grid_weight_string_list.append("   :NumberGridCells  " + str(sNcell))
    grid_weight_string_list.append("   #            ")
    grid_weight_string_list.append("   # [HRU ID] [Cell #] [w_kl]")
    max_col = max(Forcinfo["Col"].values)
    Avafgid = Forcinfo["FGID"].values
    # for id in hruids:
    #     grid_weight_string_ihru = create_grid_weight_hru(i,Mapforcing,Forcinfo)
    #     grid_weight_string_list.append(grid_weight_string_ihru)
    n_hru_group = int(len(hruids)/10)

    if len(hruids) > 1000:
        hru_ids_groups = np.array_split(hruids, n_hru_group)
    else:
        hru_ids_groups = [hruids]

    for i in range(0, len(hru_ids_groups)):
        i_hru_group = hru_ids_groups[i]
        grid_weight_string_hrusi = Parallel(n_jobs=4)(delayed(create_grid_weight_hru)(
            i, Mapforcing.copy(deep=True), Avafgid, max_col) for i in i_hru_group)
        grid_weight_string_list = grid_weight_string_list + grid_weight_string_hrusi

    grid_weight_string_list.append(":EndGridWeights")
    grid_weight_string = "\n".join(grid_weight_string_list)
    return grid_weight_string


def create_grid_weight_hru(hruid, Mapforcing, Avafgid, max_col):

    grid_weight_string_hru = ' '

    cats = Mapforcing.loc[Mapforcing["HRU_ID"] == hruid].copy(deep=True)

#    cats = cats_hru[cats_hru["Map_FGID"].isin(Avafgid)].copy(deep=True)

    if len(cats) <= 0:
        cats = Mapforcing.loc[Mapforcing["HRU_ID"] == hruid].copy(deep=True)
        print("Following Grid has to be inluded:.......")
        print(cats["Map_FGID"])
        return grid_weight_string_hru

    tarea = sum(cats["s_area"].values)
    fids = cats["Map_FGID"].values
    fids = np.unique(fids)
    sumwt = 0.0
    for j in range(0, len(fids)):
        scat = cats[cats["Map_FGID"] == fids[j]]
        if j < len(fids) - 1:
            sarea = sum(scat["s_area"].values)
            wt = float(sarea) / float(tarea)
            sumwt = sumwt + wt
        else:
            wt = 1 - sumwt

        if len(scat["Map_Row"].values) > 1:  # should be 1
            print(
                str(catid)
                + "error: 1 hru, 1 grid, produce muti polygon need to be merged "
            )
            Strcellid = (
                str(
                    int(
                        scat["Map_Row"].values[0]
                        * (max_col + 1)
                        + scat["Map_Col"].values[0]
                    )
                )
                + "      "
            )
        else:
            Strcellid = (
                str(
                    int(
                        scat["Map_Row"].values[0]
                        * (max_col + 1)
                        + scat["Map_Col"].values[0]
                    )
                )
                + "      "
            )
        if len(fids) == 1:
            grid_weight_string_hru = grid_weight_string_hru + "    " + \
                str(int(hruid)) + "     " + Strcellid + "      " + str(wt)
        else:
            if j == len(fids) - 1:
                grid_weight_string_hru = grid_weight_string_hru + "    " + \
                    str(int(hruid)) + "     " + Strcellid + "      " + str(wt)
            else:
                grid_weight_string_hru = grid_weight_string_hru + "    " + \
                    str(int(hruid)) + "     " + Strcellid + \
                    "      " + str(wt) + "\n"

    return grid_weight_string_hru
