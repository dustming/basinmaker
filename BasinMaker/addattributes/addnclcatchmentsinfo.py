import numpy as np
import pandas as pd
import copy
from addattributes.adddaandstreamorder import defcat


def update_non_connected_catchment_info(catinfo):
    routing_info = catinfo[["SubId", "DowSubId"]].astype("float").values
    catinfo_non_connected = catinfo.loc[catinfo["IsLake"] == 2].copy()

    catids_nc = catinfo_non_connected["SubId"].copy()

    catinfo.loc[
        catinfo["SubId"].isin(catids_nc), "RivLength"
    ] = 0.0  ## no reiver length since not connected.

    for i in range(0, len(catinfo_non_connected)):
        c_subid = catinfo_non_connected["SubId"].values[i]
        d_subid = catinfo_non_connected["DowSubId"].values[i]
        d_sub_info = catinfo.loc[catinfo["SubId"] == d_subid].copy()

        lc_subid = d_subid

        Upstreamcats = defcat(routing_info, c_subid)  ### alll subuds

        Up_cat_info = catinfo.loc[catinfo["SubId"].isin(Upstreamcats)].copy()

        DA = sum(Up_cat_info["BasArea"].values)

        catinfo.loc[catinfo["SubId"] == c_subid, "DA"] = DA

        if len(d_sub_info) < 1:
            continue

        ## add nonconnected lake catchment area to downsubbasin drinage area
        #        if d_sub_info['IsLake'].values[0]  != 2:
        #            catinfo.loc[catinfo['SubId'] == d_subid,'DA'] = d_sub_info['DA'].values[0] + DA

        while d_sub_info["IsLake"].values[0] == 2:

            lc_subid_info = catinfo.loc[catinfo["SubId"] == lc_subid].copy()
            d_subid = lc_subid_info["DowSubId"].values[0]
            d_sub_info = catinfo.loc[catinfo["SubId"] == d_subid].copy()
            if len(d_sub_info) < 1:
                lc_subid = -1
                break
            lc_subid = d_subid

        if lc_subid == -1:
            continue

        catinfo.loc[catinfo["SubId"] == c_subid, "RivSlope"] = d_sub_info[
            "RivSlope"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Ch_n"] = d_sub_info["Ch_n"].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Q_Mean"] = d_sub_info[
            "Q_Mean"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "BkfWidth"] = d_sub_info[
            "BkfWidth"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "BkfDepth"] = d_sub_info[
            "BkfDepth"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Strahler"] = d_sub_info[
            "Strahler"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Seg_ID"] = d_sub_info[
            "Seg_ID"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Seg_order"] = d_sub_info[
            "Seg_order"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Max_DEM"] = d_sub_info[
            "Max_DEM"
        ].values[0]
        catinfo.loc[catinfo["SubId"] == c_subid, "Min_DEM"] = d_sub_info[
            "Min_DEM"
        ].values[0]

    return catinfo
