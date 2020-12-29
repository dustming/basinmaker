import numpy as np
import pandas as pd
import copy


def defcat(out, outletid):
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
                #### if the catchment id already processed skip
                if rout[irow[j], 0] in Shedid:
                    continue
                noutshd[poshdid] = rout[irow[j], 0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd >= 0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid >= 0]
    return Shedid


def streamorderanddrainagearea(catinfoall):
    catinfo = catinfoall.loc[
        catinfoall["IsLake"] != 2
    ].copy()  ### remove none connected lake catchments, which do not connected to the river system
    catinfo_ncl = catinfoall.loc[catinfoall["IsLake"] == 2].copy()
    routing_ncl = catinfo_ncl[["SubId", "DowSubId"]].astype("float").values

    catlist = np.full((len(catinfo)), -9)
    icat = 0
    iseg = 1
    ### find first segments of all reaches, no upstream reaches
    for i in range(0, len(catinfo)):
        idx = catinfo.index[i]
        if catinfo["SubId"].values[i] == catinfo["DowSubId"].values[i]:
            catinfo.loc[idx, "DowSubId"] = -1
        catid = catinfo["SubId"].values[i]
        if (
            len(catinfo[catinfo["DowSubId"] == catid]) == 0
        ):  ### the river seg has no upstream segment
            catlist[icat] = int(
                catinfo["DowSubId"].values[i]
            )  #### store next reach segment

            #### calculate DA of head watershed include None connected lakes
            if len(routing_ncl) == 0:
                DA_ncl = 0.0
            else:
                Upstreamcats = defcat(routing_ncl, catid)  ### alll subuds
                Up_cat_info = catinfo_ncl.loc[
                    catinfo_ncl["SubId"].isin(Upstreamcats)
                ].copy()

                if len(Up_cat_info) > 0:
                    DA_ncl = sum(Up_cat_info["BasArea"].values)
                else:
                    DA_ncl = 0.0

            catinfo.loc[idx, "DA"] = DA_ncl + catinfo["BasArea"].values[i]
            catinfo.loc[idx, "Strahler"] = 1
            catinfo.loc[idx, "Seg_order"] = 1
            catinfo.loc[idx, "Seg_ID"] = iseg
            icat = icat + 1
            iseg = iseg + 1

    catlist = np.unique(catlist)
    catlist = catlist[catlist > 0]
    #    print(catlist)
    ### Loop for each first reach, until go to reach intersection
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

            #### calculate DA of None connected lakes
            if len(routing_ncl) == 0:
                DA_ncl = 0.0
            else:
                #                print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
                #                print(catid)
                #                print(routing_ncl)
                Upstreamcats = defcat(routing_ncl, catid)  ### alll subuds
                Up_cat_info = catinfo_ncl.loc[
                    catinfo_ncl["SubId"].isin(Upstreamcats)
                ].copy()
                if len(Up_cat_info) > 0:
                    DA_ncl = sum(Up_cat_info["BasArea"].values)
                else:
                    DA_ncl = 0.0

            if (
                len(cur_Reach_info) <= 0
            ):  ### reach the most downstream of the watersheds
                break

            if len(Up_Reaches_info) == 1:  ### only have one upstream
                catinfo.loc[curcat_idx, "DA"] = (
                    cur_Reach_info["BasArea"].values[0]
                    + Up_Reaches_info["DA"].values[0]
                    + DA_ncl
                )
                catinfo.loc[curcat_idx, "Strahler"] = Up_Reaches_info[
                    "Strahler"
                ].values[0]
                catinfo.loc[curcat_idx, "Seg_order"] = (
                    Up_Reaches_info["Seg_order"].values[0] + 1
                )
                catinfo.loc[curcat_idx, "Seg_ID"] = Up_Reaches_info["Seg_ID"].values[0]
                #                print('1',catid,catinfo.loc[curcat_idx,'DA'].values,catinfo.loc[curcat_idx,'Strahler'].values,catinfo.loc[curcat_idx,'Sub_order'].values)
                catid = int(cur_Reach_info["DowSubId"].values[0])
            else:  ### has mutiple upstram
                if (
                    np.min(Up_Reaches_info["Strahler"].values) > 0
                ):  ### all upstream has been processed
                    catinfo.loc[catinfo["SubId"] == catid, "DA"] = (
                        cur_Reach_info["BasArea"].values[0]
                        + np.sum(Up_Reaches_info["DA"].values)
                        + DA_ncl
                    )
                    if np.min(Up_Reaches_info["Strahler"].values) == np.max(
                        Up_Reaches_info["Strahler"].values
                    ):  ### two reach has the same order
                        catinfo.loc[curcat_idx, "Strahler"] = (
                            Up_Reaches_info["Strahler"].values[0] + 1
                        )
                        catinfo.loc[curcat_idx, "Seg_order"] = 1
                        catinfo.loc[curcat_idx, "Seg_ID"] = iseg + 1
                        iseg = iseg + 1
                    #                        print('2',catid,catinfo.loc[catinfo['SubId'] == catid,'DA'].values,catinfo.loc[catinfo['SubId'] == catid,'Strahler'].values,catinfo.loc[catinfo['SubId'] == catid,'Sub_order'].values)
                    else:
                        max_straorder = np.max(Up_Reaches_info["Strahler"].values)
                        catinfo.loc[curcat_idx, "Strahler"] = max_straorder
                        catinfo.loc[curcat_idx, "Seg_order"] = 1
                        catinfo.loc[curcat_idx, "Seg_ID"] = iseg + 1
                        iseg = iseg + 1
                    #                        print('3',catid,catinfo.loc[catinfo['SubId'] == catid,'DA'].values,catinfo.loc[catinfo['SubId'] == catid,'Strahler'].values,catinfo.loc[catinfo['SubId'] == catid,'Sub_order'].values)
                    catid = int(cur_Reach_info["DowSubId"].values[0])
                else:  ## there are some reach has not been processed, save id to the list and wait for another loob
                    newcatlist[inewstart] = int(catid)
                    inewstart = inewstart + 1
                    F_intersect = 0

    mask = catinfoall["SubId"].isin(catinfo["SubId"].values)
    catinfoall.loc[mask, "Seg_ID"] = catinfo["Seg_ID"].values
    catinfoall.loc[mask, "Seg_order"] = catinfo["Seg_order"].values
    catinfoall.loc[mask, "Strahler"] = catinfo["Strahler"].values
    catinfoall.loc[mask, "Seg_ID"] = catinfo["Seg_ID"].values
    catinfoall.loc[mask, "DA"] = catinfo["DA"].values

    return catinfoall
