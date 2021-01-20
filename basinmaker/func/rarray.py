import copy

import numpy as np


def Is_Point_Close_To_Id_In_Raster(prow, pcol, nrows, ncols, id, raster_array):
    """Check if the point is around grids with value equal to Id in raster_array

    Parameters
    ----------

    Returns:
    -------
    Isclose               : logical
       True  : close to grids with value euqal to id
       False : not close to grids with value to id
    """
    Isclose = False
    n_grids_eq_id = 0
    if prow != 0 and prow != nrows - 1 and pcol != 0 and pcol != ncols - 1:
        if raster_array[prow - 1, pcol + 1] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow - 1, pcol - 1] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow - 1, pcol] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow, pcol + 1] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow, pcol - 1] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow + 1, pcol - 1] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow + 1, pcol + 1] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1
        if raster_array[prow + 1, pcol] == id:
            Isclose = True
            n_grids_eq_id = n_grids_eq_id + 1

    return Isclose, n_grids_eq_id


def CE_mcat4lake2(
    cat1, lake, fac, fdir, bsid, nrows, ncols, Pourpoints, noncnlake, str_array
):
    cat = copy.copy(cat1)
    Non_con_lake_cat = copy.copy(cat1)
    Non_con_lake_cat[:, :] = -9999
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid > 0]
    noncnlake = np.unique(noncnlake)
    noncnlake = noncnlake[noncnlake > 0]
    for i in range(0, len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake == lakeid).astype(int)
        lakacc = np.full((len(lrowcol), 3), -9999)
        lakacc[:, 0] = lrowcol[:, 0]
        lakacc[:, 1] = lrowcol[:, 1]
        lakacc[:, 2] = fac[lrowcol[:, 0], lrowcol[:, 1]]
        lakacc = lakacc[lakacc[:, 2].argsort()]
        lorow = lakacc[len(lakacc) - 1, 0]
        locol = lakacc[len(lakacc) - 1, 1]  ###### lake outlet row and col
        arclakeid = cat1[lorow, locol]
        pp = Pourpoints[lorow, locol]
        pp = np.unique(pp)
        pp = pp[pp > 0]
        #        print(lakeid,len(np.argwhere(noncnlake==lakeid)),pp,Pourpoints[lorow,locol],arclakeid)
        if len(pp) == 1:
            if arclakeid < 0:
                cat[lrowcol[:, 0], lrowcol[:, 1]] = pp
            else:
                cat[lrowcol[:, 0], lrowcol[:, 1]] = arclakeid
            if (
                len(np.argwhere(noncnlake == lakeid)) > 0
            ):  ##if lake belong to non connected lakes
                if arclakeid < 0:
                    nonlrowcol = np.argwhere(cat == pp).astype(int)
                    strids = np.unique(str_array[nonlrowcol[:, 0], nonlrowcol[:, 1]])
                    if (
                        len(strids) <= 1
                    ):  ## if none connected lake catchment overlay with stream, remove this lake
                        Non_con_lake_cat[nonlrowcol[:, 0], nonlrowcol[:, 1]] = lakeid
                    else:
                        cat[lrowcol[:, 0], lrowcol[:, 1]] = -99
                else:
                    nonlrowcol = np.argwhere(cat == arclakeid).astype(int)
                    strids = np.unique(str_array[nonlrowcol[:, 0], nonlrowcol[:, 1]])

                    if (
                        len(strids) <= 1
                    ):  ## if none connected lake catchment overlay with stream, remove this lake
                        Non_con_lake_cat[nonlrowcol[:, 0], nonlrowcol[:, 1]] = lakeid
                    else:
                        cat[lrowcol[:, 0], lrowcol[:, 1]] = -99
    #### For some reason pour point was missing in non-contribute catchment
    Pours = np.unique(Pourpoints)
    Pours = Pours[Pours > 0]
    for i in range(0, len(Pours)):
        pourid = Pours[i]
        rowcol = Pourpoints == pourid
        if cat[rowcol] < 0:
            temp_notused, nout = Is_Point_Close_To_Id_In_Raster(
                rowcol[0, 0], rowcol[0, 1], nrows, ncols, pourid, cat
            )
            if len(cat[cat == pourid]) > 0 and nout < 8:
                cat[rowcol] = pourid
    rowcol1 = fac > 0
    rowcol2 = cat < 0
    noncontribuite = np.logical_and(rowcol1, rowcol2)
    cat[noncontribuite] = 2 * max(np.unique(cat)) + 1
    return cat, Non_con_lake_cat


def GenerateFinalPourpoints(
    fac, fdir, lake, cat3, bsid, blid, boid, nrows, ncols, cat, obs, Is_divid_region
):
    Poups = copy.copy(cat3)
    Poups[:, :] = -9999
    GWat = copy.copy(cat3)
    GWatids = np.unique(cat3)
    GWatids = GWatids[GWatids >= 0]
    ncatid = 1
    for i in range(0, len(GWatids)):
        trow, tcol = Getbasinoutlet(GWatids[i], GWat, fac, fdir, nrows, ncols)
        #        arcpy.AddMessage("result     " +str(trow) +"    " +str(tcol))
        Poups[trow, tcol] = ncatid
        ncatid = ncatid + 1
    OWat = copy.copy(cat)
    OWatids = np.unique(cat)
    OWatids = OWatids[OWatids >= 0]
    for i in range(0, len(OWatids)):
        trow, tcol = Getbasinoutlet(OWatids[i], OWat, fac, fdir, nrows, ncols)
        if not GWat[trow, tcol] >= blid:
            if Poups[trow, tcol] < 0:
                Poups[trow, tcol] = ncatid
                ncatid = ncatid + 1
    if Is_divid_region > 0:
        return Poups
    obsids = np.unique(obs)
    obsids = obsids[obsids > 0]
    for i in range(0, len(obsids)):
        rowcol = np.argwhere(obs == obsids[i]).astype(int)
        if (
            Poups[rowcol[0, 0], rowcol[0, 1]] < 0
            and lake[rowcol[0, 0], rowcol[0, 1]] < 0
        ):
            Poups[rowcol[0, 0], rowcol[0, 1]] = ncatid
            ncatid = ncatid + 1
    return Poups


def Check_If_Str_Is_Head_Stream(prow, pcol, nrows, ncols, str, strid):

    noout = True
    ### the point prow and pcol is not surrounded by oter strems, except
    ### stream with stream id = strid
    if prow != 0 and prow != nrows - 1 and pcol != 0 and pcol != ncols - 1:

        if str[prow - 1, pcol + 1] > 0 and str[prow - 1, pcol + 1] != strid:
            noout = False
        if str[prow - 1, pcol - 1] > 0 and str[prow - 1, pcol - 1] != strid:
            noout = False
        if str[prow - 1, pcol] > 0 and str[prow - 1, pcol] != strid:
            noout = False
        if str[prow, pcol + 1] > 0 and str[prow, pcol + 1] != strid:
            noout = False
        if str[prow, pcol - 1] > 0 and str[prow, pcol - 1] != strid:
            noout = False
        if str[prow + 1, pcol - 1] > 0 and str[prow + 1, pcol - 1] != strid:
            noout = False
        if str[prow + 1, pcol + 1] > 0 and str[prow + 1, pcol + 1] != strid:
            noout = False
        if str[prow + 1, pcol] > 0 and str[prow + 1, pcol] != strid:
            noout = False
    return noout


###
def GenerPourpoint(
    cat,
    lake,
    Str,
    nrows,
    ncols,
    blid,
    bsid,
    bcid,
    fac,
    hydir,
    Is_divid_region,
    Remove_Str,
    Min_Grid_Number=50,
):
    GP_cat = copy.copy(cat)
    sblid = copy.copy(blid)
    Str_new = copy.copy(Str)
    ############### Part 1 Get all pourpoints of hydroshed catchment
    arcatid = np.unique(cat)  #### cat all catchment idd
    arcatid = arcatid[arcatid >= 0]
    catoutloc = np.full((len(arcatid), 3), -9999)
    for i in range(0, len(arcatid)):
        catid = arcatid[i]
        catrowcol = np.argwhere(cat == catid).astype(int)
        trow, tcol = Getbasinoutlet(catid, cat, fac, hydir, nrows, ncols)
        GP_cat[
            catrowcol[:, 0], catrowcol[:, 1]
        ] = -9999  ### set the catment cells into null
        GP_cat[trow, tcol] = bcid  #### change the outlet of catchment into wid
        bcid = bcid + 1
        catoutloc[i, 0] = catid  ## catchment id
        catoutloc[i, 1] = trow  #### catchment pourpont row
        catoutloc[i, 2] = tcol  #### catchment pourpont col
    #    writeraster(outFolder+subid+"_Pourpoints_1.asc",GP_cat)

    #### make the end point of stream in the head water shed to be a pourpoints

    Strids = np.unique(Str)
    Strids = Strids[Strids > 0]
    newstr_id = max(Strids) + 10
    for i in range(0, len(Strids)):
        strid = Strids[i]
        strrowcol = np.argwhere(Str == strid).astype(int)
        nstrrow = strrowcol.shape[0]

        ### if number of the stream grid smaller than 50

        if nstrrow < Min_Grid_Number:
            continue
        iStr = np.full((nstrrow, 4), -9999)  ##### 0 row, 1 col, 2 fac,
        iStr[:, 0] = strrowcol[:, 0]
        iStr[:, 1] = strrowcol[:, 1]
        iStr[:, 2] = fac[strrowcol[:, 0], strrowcol[:, 1]]
        iStr = iStr[iStr[:, 2].argsort()].astype(int)
        ### starting point of the stream
        Start_Pt_row = iStr[0, 0]
        Start_Pt_col = iStr[0, 1]

        Is_Head_Stream = Check_If_Str_Is_Head_Stream(
            Start_Pt_row, Start_Pt_col, nrows, ncols, Str, strid
        )
        if Is_Head_Stream == True:
            Str_new[
                iStr[0, 0], iStr[0, 1]
            ] = newstr_id  #### change the outlet of catchment into wid
            Str_new[iStr[1, 0], iStr[1, 1]] = newstr_id
            Str_new[iStr[2, 0], iStr[2, 1]] = newstr_id
            Str_new[iStr[3, 0], iStr[3, 1]] = newstr_id
            Str_new[iStr[4, 0], iStr[4, 1]] = newstr_id
            newstr_id = newstr_id + 1
    ######
    ##################Part 2 Get pourpoints of Lake inflow streams
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid >= 0]
    Lakemorestream = copy.copy(arlakeid)
    Lakemorestream[:] = -9999
    idxlakemorestream = 0
    for i in range(0, len(arlakeid)):  #### loop for each lake
        lid = arlakeid[i]  ##### obtain lake id
        rowcol = np.argwhere(lake == lid.astype(int))  ### got row anc col of lake grids
        nrow = rowcol.shape[0]  ### number of grids of the lake
        Stridinlake = np.full(nrow, -9999)  ####
        Stridinlake[:] = Str[
            rowcol[:, 0], rowcol[:, 1]
        ]  ### get all stream with in the lake domain
        Strid_L = np.unique(
            Stridinlake[np.argwhere(Stridinlake > 0).astype(int)]
        )  ### Get unique stream id of sach stream
        ##### find the intercept point of stream and lake
        if Is_divid_region > 0:
            if len(Strid_L) <= 1:
                Lakemorestream[idxlakemorestream] = lid
                idxlakemorestream = idxlakemorestream + 1
                continue

        for j in range(0, len(Strid_L)):  #### loop for each stream intercept with lake
            strid = Strid_L[j]
            ### get [row,col] of strid in str grids
            strrowcol = np.argwhere(Str == strid).astype(int)
            ### get number of grids of this str
            nstrrow = strrowcol.shape[0]
            ### create an array 0 will store row 1 will be col and 2 will be acc
            ### 3 will be the lake id
            Strchek = np.full((nstrrow, 4), -9999)  ##### 0 row, 1 col, 2 fac,
            Strchek[:, 0] = strrowcol[:, 0]
            Strchek[:, 1] = strrowcol[:, 1]
            Strchek[:, 2] = fac[strrowcol[:, 0], strrowcol[:, 1]]
            Strchek[:, 3] = lake[strrowcol[:, 0], strrowcol[:, 1]]
            ### sort the grids by flow accumulation
            Strchek = Strchek[Strchek[:, 2].argsort()].astype(int)

            ### find the first intersection point between river and lake
            Grids_Lake_river = np.where(Strchek[:, 3] == lid)
            irowst = np.nanmin(Grids_Lake_river)
            Intersection_Grid_row = Strchek[irowst, 0]
            Intersection_Grid_col = Strchek[irowst, 1]

            noout = 0
            ### double check if the head stream cell is nearby the lake
            if (
                Strchek[0, 0] != 0
                and Strchek[0, 0] != nrows - 1
                and Strchek[0, 1] != 0
                and Strchek[0, 1] != ncols - 1
            ):
                noout, temp_notused = Is_Point_Close_To_Id_In_Raster(
                    Strchek[0, 0], Strchek[0, 1], nrows, ncols, lid, lake
                )
            ####

            if irowst != 0:  ## is the stream is not start within the lake
                if (
                    lake[Strchek[irowst - 1, 0], Strchek[irowst - 1, 1]] == -9999
                ):  ### double check
                    ##### this means the stream connect two lakes, so must assign an pourpoints
                    if (
                        len(np.unique(lake[Strchek[0:irowst, 0], Strchek[0:irowst, 1]]))
                        >= 3
                    ):
                        GP_cat[Strchek[irowst - 1, 0], Strchek[irowst - 1, 1]] = bsid
                        bsid = bsid + 1

                    ##### the head stream celll is not nearby the lake
                    if noout == False:
                        GP_cat[Strchek[irowst - 1, 0], Strchek[irowst - 1, 1]] = bsid
                        bsid = bsid + 1
            #### it is possible that two steam combine together near the lake, double check if the stream conncet to
            # anotehr stream and this steam is not witin the lake
            if irowst == 0 or noout == True:
                nostr = Str[Strchek[0, 0], Strchek[0, 1]]
                a = 0
                orowcol = np.full((8, 3), -9999)
                if (
                    Strchek[0, 0] != 0
                    and Strchek[0, 0] != nrows - 1
                    and Strchek[0, 1] != 0
                    and Strchek[0, 1] != ncols - 1
                ):
                    if (
                        Str[Strchek[0, 0] - 1, Strchek[0, 1] + 1] != -9999
                        and lake[Strchek[0, 0] - 1, Strchek[0, 1] + 1] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] - 1
                        orowcol[a, 1] = Strchek[0, 1] + 1
                        orowcol[a, 2] = Str[Strchek[0, 0] - 1, Strchek[0, 1] + 1]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0] - 1, Strchek[0, 1] - 1] != -9999
                        and lake[Strchek[0, 0] - 1, Strchek[0, 1] - 1] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] - 1
                        orowcol[a, 1] = Strchek[0, 1] - 1
                        orowcol[a, 2] = Str[Strchek[0, 0] - 1, Strchek[0, 1] - 1]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0] - 1, Strchek[0, 1]] != -9999
                        and lake[Strchek[0, 0] - 1, Strchek[0, 1]] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] - 1
                        orowcol[a, 1] = Strchek[0, 1] - 0
                        orowcol[a, 2] = Str[Strchek[0, 0] - 1, Strchek[0, 1] - 0]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0], Strchek[0, 1] + 1] != -9999
                        and lake[Strchek[0, 0], Strchek[0, 1] + 1] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] - 0
                        orowcol[a, 1] = Strchek[0, 1] + 1
                        orowcol[a, 2] = Str[Strchek[0, 0] - 0, Strchek[0, 1] + 1]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0], Strchek[0, 1] - 1] != -9999
                        and lake[Strchek[0, 0], Strchek[0, 1] - 1] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] - 0
                        orowcol[a, 1] = Strchek[0, 1] - 1
                        orowcol[a, 2] = Str[Strchek[0, 0] - 0, Strchek[0, 1] - 1]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0] + 1, Strchek[0, 1] - 1] != -9999
                        and lake[Strchek[0, 0] + 1, Strchek[0, 1] - 1] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] + 1
                        orowcol[a, 1] = Strchek[0, 1] - 1
                        orowcol[a, 2] = Str[Strchek[0, 0] + 1, Strchek[0, 1] - 1]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0] + 1, Strchek[0, 1] + 1] != -9999
                        and lake[Strchek[0, 0] + 1, Strchek[0, 1] + 1] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] + 1
                        orowcol[a, 1] = Strchek[0, 1] + 1
                        orowcol[a, 2] = Str[Strchek[0, 0] + 1, Strchek[0, 1] + 1]
                        a = a + 1
                    if (
                        Str[Strchek[0, 0] + 1, Strchek[0, 1]] != -9999
                        and lake[Strchek[0, 0] + 1, Strchek[0, 1]] == -9999
                    ):
                        orowcol[a, 0] = Strchek[0, 0] + 1
                        orowcol[a, 1] = Strchek[0, 1] - 0
                        orowcol[a, 2] = Str[Strchek[0, 0] + 1, Strchek[0, 1] - 0]
                        a = a + 1
                if a > 0:
                    for ka in range(0, a):
                        nostr = orowcol[ka, 2]  ## up stream stram id
                        srowcol = np.argwhere(Str == nostr).astype(int)
                        snrow = srowcol.shape[0]
                        iStrchek = np.full(
                            (snrow, 4), -9999
                        )  ##### 0 row, 1 col, 2 fac,
                        iStrchek[:, 0] = srowcol[:, 0]
                        iStrchek[:, 1] = srowcol[:, 1]
                        iStrchek[:, 2] = fac[srowcol[:, 0], srowcol[:, 1]]
                        iStrchek = iStrchek[iStrchek[:, 2].argsort()]
                        noout, temp_notused = Is_Point_Close_To_Id_In_Raster(
                            iStrchek[0, 0], iStrchek[0, 1], nrows, ncols, lid, lake
                        )
                        Lakinstr = np.full(snrow, -9999)
                        Lakinstr[:] = lake[srowcol[:, 0], srowcol[:, 1]]
                        d = np.argwhere(Lakinstr == lid).astype(
                            int
                        )  #### the connected stream should not within the lake
                        if len(d) < 1 and noout == False:
                            GP_cat[orowcol[ka, 0], orowcol[ka, 1]] = bsid
                            bsid = bsid + 1
    ################################################################################
    ################## Part 3Get Lake pourpoint id and remove cat pourpoint that contribute to lake
    #    writeraster(outFolder+"Pourpoints_2.asc",GP_cat)
    Lakemorestream = Lakemorestream[Lakemorestream > 0]
    for i in range(0, len(arlakeid)):
        lakeid = arlakeid[i]

        if Is_divid_region > 0:
            mask = Lakemorestream == lakeid
            if np.sum(mask) >= 1:
                continue

        lrowcol = np.argwhere(lake == lakeid).astype(int)
        arcatid = np.unique(
            cat[lrowcol[:, 0], lrowcol[:, 1]]
        )  ## Get all catchment that intercept with lake
        lakeacc = np.full((len(lrowcol), 3), -9999)
        lakeacc[:, 0] = lrowcol[:, 0]
        lakeacc[:, 1] = lrowcol[:, 1]
        lakeacc[:, 2] = fac[lrowcol[:, 0], lrowcol[:, 1]]
        lakeacc = lakeacc[lakeacc[:, 2].argsort()]
        maxcatid = cat[
            lakeacc[len(lakeacc) - 1, 0], lakeacc[len(lakeacc) - 1, 1]
        ]  ### Get the catchment id that will keeped

        ### remove the all pourpoints close to lake pourpoints
        if (
            lakeacc[len(lakeacc) - 1, 0] != 0
            and lakeacc[len(lakeacc) - 1, 0] != nrows - 1
            and lakeacc[len(lakeacc) - 1, 1] != 0
            and lakeacc[len(lakeacc) - 1, 1] != ncols - 1
        ):
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0] - 1, lakeacc[len(lakeacc) - 1, 1] - 1
            ] = -9999
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0] + 1, lakeacc[len(lakeacc) - 1, 1] - 1
            ] = -9999
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0], lakeacc[len(lakeacc) - 1, 1] - 1
            ] = -9999
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0] + 1, lakeacc[len(lakeacc) - 1, 1]
            ] = -9999
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0] + 1, lakeacc[len(lakeacc) - 1, 1] + 1
            ] = -9999
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0], lakeacc[len(lakeacc) - 1, 1] + 1
            ] = -9999
            GP_cat[
                lakeacc[len(lakeacc) - 1, 0] - 1, lakeacc[len(lakeacc) - 1, 1]
            ] = -9999

        ### remove pourpoints that not equal to maxcatid
        for j in range(0, len(arcatid)):
            if arcatid[j] != maxcatid:
                crowcol = np.argwhere(cat == arcatid[j]).astype(int)
                checkcat = np.full((len(crowcol), 4), -9999)
                checkcat[:, 0] = crowcol[:, 0]
                checkcat[:, 1] = crowcol[:, 1]
                checkcat[:, 2] = GP_cat[crowcol[:, 0], crowcol[:, 1]]
                checkcat[:, 3] = fac[crowcol[:, 0], crowcol[:, 1]]
                checkcat = checkcat[checkcat[:, 3].argsort()]
                dele = checkcat[checkcat[:, 2] < blid].astype(
                    int
                )  #### do not delete the lake and strem ids
                GP_cat[dele[:, 0], dele[:, 1]] = -9999

                ### add keep catchment in case the catchment outlet is at the boundary of the domain
                nrow, ncol = Nextcell(
                    hydir,
                    checkcat[len(checkcat) - 1, 0],
                    checkcat[len(checkcat) - 1, 1],
                )
                if nrow > 0 or ncol > 0:
                    if nrow >= nrows or ncols >= ncols:
                        continue
                    if cat[nrow, ncol] < 0:
                        GP_cat[
                            checkcat[len(checkcat) - 1, 0],
                            checkcat[len(checkcat) - 1, 1],
                        ] = bcid
                        bcid = bcid + 1

        #### add lake outlet
        GP_cat[lakeacc[len(lakeacc) - 1, 0], lakeacc[len(lakeacc) - 1, 1]] = sblid
        sblid = sblid + 1
    return GP_cat, Lakemorestream, Str_new


def Getbasinoutlet(ID, basin, fac, dir, nrows, ncols):
    catrowcol = np.argwhere(basin == ID).astype(int)
    catacc = np.full((len(catrowcol), 3), -9999)
    catacc[:, 0] = catrowcol[:, 0]
    catacc[:, 1] = catrowcol[:, 1]
    catacc[:, 2] = fac[catrowcol[:, 0], catrowcol[:, 1]]
    catacc = catacc[catacc[:, 2].argsort()]
    #    print(ID,len(catrowcol))
    ### check if it is a real basin outlet
    crow = catacc[len(catrowcol) - 1, 0]
    ccol = catacc[len(catrowcol) - 1, 1]

    nrow, ncol = Nextcell(dir, crow, ccol)
    #    print(ID,basin[nrow,ncol],basin[crow,ccol],fac[nrow,ncol],fac[crow,ccol],crow,ccol)
    if nrow < 0 or ncol < 0:
        return crow, ccol
    elif nrow >= nrows or ncol >= ncols:
        return crow, ccol
    elif basin[nrow, ncol] < 0:
        return crow, ccol
    elif basin[nrow, ncol] != ID:  #  all above means the outlet is the real loutlet
        return crow, ccol
    else:
        crow = nrow
        ccol = ncol
        ifound = 0
        for i in range(0, 1000):  #### find next 1000 grids, to find the basin outlet
            nrow, ncol = Nextcell(dir, crow, ccol)
            if nrow < 0 or ncol < 0:
                ifound = 1
                break
            elif nrow >= nrows or ncol >= ncols:
                ifound = 1
                break
            elif basin[nrow, ncol] < 0:
                ifound = 1
                break
            elif basin[nrow, ncol] != ID:
                ifound = 1  #     all above means the outlet is the real loutlet
                break
            else:
                crow = nrow
                ccol = ncol
                continue
        if ifound == 0:
            print(" true basin outlet not found for ID...." + str(ID))
        return crow, ccol


def Nextcell(N_dir, N_row, N_col):
    if N_dir[N_row, N_col] == 1:
        N_nrow = N_row + 0
        N_ncol = N_col + 1
    elif N_dir[N_row, N_col] == 2:
        N_nrow = N_row + 1
        N_ncol = N_col + 1
    elif N_dir[N_row, N_col] == 4:
        N_nrow = N_row + 1
        N_ncol = N_col + 0
    elif N_dir[N_row, N_col] == 8:
        N_nrow = N_row + 1
        N_ncol = N_col - 1
    elif N_dir[N_row, N_col] == 16:
        N_nrow = N_row + 0
        N_ncol = N_col - 1
    elif N_dir[N_row, N_col] == 32:
        N_nrow = N_row - 1
        N_ncol = N_col - 1
    elif N_dir[N_row, N_col] == 64:
        N_nrow = N_row - 1
        N_ncol = N_col + 0
    elif N_dir[N_row, N_col] == 128:
        N_nrow = N_row - 1
        N_ncol = N_col + 1
    else:
        N_nrow = -9999
        N_ncol = -9999
    return N_nrow, N_ncol


def return_subid_of_next_down_stream_grids(
    cat_array, catid, nfdr_arcgis_array, cat_outlet_array, ncols, nrows
):
    # get outlet row cols
    maxtry = 10
    outlet_row_col = np.argwhere(cat_outlet_array == catid)
    if len(outlet_row_col) <= 0:
        return -1
    outlet_row = outlet_row_col[0, 0]
    outlet_col = outlet_row_col[0, 1]
    downcatid = -1
    for i in range(0, maxtry):
        nrow, ncol = Nextcell(nfdr_arcgis_array, outlet_row, outlet_col)
        if nrow >= nrows or ncol >= ncols:
            downsubid = -1
            return downsubid
        elif nrow <= 0 or ncol <= 0:
            downsubid = -1
            return downsubid
        else:
            downcatid = cat_array[nrow, ncol]
            if downcatid != catid:
                return downcatid
            else:
                outlet_row = nrow
                outlet_col = ncol
    return downcatid
