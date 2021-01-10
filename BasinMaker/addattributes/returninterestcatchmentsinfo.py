import numpy as np
import pandas as pd
import copy
from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
from addattributes.adddaandstreamorder import defcat


def return_interest_catchments_info(catinfo, outlet_obs_id, path_sub_reg_outlets_v="#"):

    if outlet_obs_id < 0:
        return catinfo

    routing_info = catinfo[["SubId", "DowSubId"]].astype("float").values

    if path_sub_reg_outlets_v != "#":

        Sub_reg_outlets = Dbf_To_Dataframe(path_sub_reg_outlets_v)["value"].values
        Sub_reg_outlets = np.unique(Sub_reg_outlets)
        Sub_reg_outlets_ids = Sub_reg_outlets_ids[Sub_reg_outlets_ids > 0]

        #### Find all obervation id that is subregion outlet
        reg_outlet_info = catinfo.loc[catinfo["IsObs"].isin(Sub_reg_outlets_ids)]

        #### Define outlet ID
        outletid = -1

        if outlet_obs_id < 0:
            print("To use subregion, the Subregion Id MUST provided as Outlet_Obs_ID")
            return catinfo

        outletID_info = catinfo.loc[catinfo["IsObs"] == outlet_obs_id]

        if len(outletID_info) > 0:
            outletid = outletID_info["SubId"].values[0]
        else:
            print("No Outlet id is founded for subregion   ", Outlet_Obs_ID)
            return catinfo

        ### find all subregion drainge to this outlet id
        HydroBasins1 = defcat(routing_info, outletid)

        ### if there is other subregion outlet included in current sturcture
        ### remove subbasins drainge to them

        if len(reg_outlet_info) >= 2:  ### has upstream regin outlet s
            ### remove subbains drainage to upstream regin outlet s
            for i in range(0, len(reg_outlet_info)):
                upregid = reg_outlet_info["SubId"].values[i]
                ### the subregion ouetlet not within the target domain neglect
                if upregid == outletid or np.sum(np.in1d(HydroBasins1, upregid)) < 1:
                    continue
                HydroBasins_remove = Defcat(routing_info, upregid)
                mask = np.in1d(
                    HydroBasins1, HydroBasins_remove
                )  ### exluced ids that belongs to main river stream
                HydroBasins1 = HydroBasins1[np.logical_not(mask)]

        HydroBasins = HydroBasins1

        catinfo = catinfo.loc[catinfo["SubId"].isin(HydroBasins)]
        return catinfo
    ### selected based on observation guage obs id
    else:
        outletid = -1
        outletID_info = catinfo.loc[catinfo["IsObs"] == outlet_obs_id]
        if len(outletID_info) > 0:
            outletid = outletID_info["SubId"].values[0]
            ##find upsteam catchment id
            HydroBasins = defcat(routing_info, outletid)
        else:
            HydroBasins = catinfo["SubId"].values

        catinfo = catinfo.loc[catinfo["SubId"].isin(HydroBasins)]
        return catinfo
