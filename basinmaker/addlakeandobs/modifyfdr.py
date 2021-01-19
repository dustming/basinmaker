import copy
import os

import numpy as np


def modify_lakes_flow_direction(
    cat3,
    lake,
    fac,
    dir,
    str_array,
    nrows,
    ncols,
    LakeBD_array,
    Pec_Grid_outlier,
    MaximumLakegrids,
    Lakemorestream,
):
    cat = copy.copy(cat3)
    ndir = copy.copy(dir)
    changed_ndir = copy.copy(dir)
    changed_ndir[:, :] = -9999
    BD_problem = copy.copy(dir)
    BD_problem[:, :] = -9999
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid > 0]
    outlakeids = np.full((1000000, 2), -99999.999)
    stream_mask = str_array > 0
    for i in range(0, len(arlakeid)):
        lakeid = arlakeid[i]
        if lakeid in Lakemorestream:
            continue

        lrowcol = np.argwhere(lake == lakeid).astype(int)
        lakacc = np.full((len(lrowcol), 3), -9999)
        lakacc[:, 0] = lrowcol[:, 0]
        lakacc[:, 1] = lrowcol[:, 1]
        lakacc[:, 2] = fac[lrowcol[:, 0], lrowcol[:, 1]]
        lakacc = lakacc[lakacc[:, 2].argsort()]
        lorow = lakacc[len(lakacc) - 1, 0]
        locol = lakacc[len(lakacc) - 1, 1]  ###### lake outlet row and col
        arclakeid = cat[lorow, locol]  ####### lake catchment id
        lakecatrowcol = np.argwhere(cat == arclakeid).astype(
            int
        )  #### lake catchment cells
        Lakeincat1 = lake[lakecatrowcol[:, 0], lakecatrowcol[:, 1]]
        nlake = np.argwhere(Lakeincat1 == lakeid).astype(int)
        outlakeids[i, 0] = lakeid
        outlakeids[i, 1] = float(len(nlake) / len(lrowcol))
        #        print("########################################################################3")
        #        print(lakeid,arclakeid,len(nlake),len(lrowcol),float(len(nlake)/len(lrowcol)),Pec_Grid_outlier)
        if outlakeids[i, 1] > Pec_Grid_outlier:
            continue

        #        if outlakeids[i,1] > Pec_Grid_outlier: ### smaller than 0.97
        #            continue

        #        print(lakeid,arclakeid,len(nlake),len(lrowcol),float(len(nlake)/len(lrowcol)))
        BD_mask = LakeBD_array == lakeid
        Lake_mask = lake == lakeid
        cat_lake_mask = cat == arclakeid

        Lakeincat_mask = np.logical_and(Lake_mask, cat_lake_mask)
        nlake2 = np.sum(Lakeincat_mask)

        Lakeoutcat_mask = np.logical_and(Lake_mask, np.logical_not(cat_lake_mask))

        BD_Out_Lakecat_mask = np.logical_and(Lakeoutcat_mask, BD_mask)

        #        BD_Out_Lakecat_Nriv_mask = np.logical_and(BD_Out_Lakecat_mask,np.logical_not(stream_mask)) do not allow modify stream grids
        BD_Out_Lakecat_Nriv_mask = BD_Out_Lakecat_mask
        BD_problem[BD_Out_Lakecat_Nriv_mask] = 1

        if (
            np.sum(BD_Out_Lakecat_mask) > MaximumLakegrids
        ):  ### smaller than nlakegrids or smaller than 0.9
            continue

        print("Total # of lakes: ", len(arlakeid), " Processing    ", i, "th lake")
        print(
            "Lake ID : ",
            lakeid,
            "Lake Cat ID   ",
            arclakeid,
            "Total numer of Lake grids   ",
            len(lrowcol),
            "Numer of Lake grids in Lake Cat:  ",
            nlake2,
            len(nlake),
        )
        print(
            "# of Lake boundary grids:   ",
            np.sum(BD_mask),
            "# of grids do not flow to lake catchments",
            np.sum(Lakeoutcat_mask),
            "# of lake boundary grids not flow to lake catchment   ",
            np.sum(BD_Out_Lakecat_mask),
            "  # of lake boundary grids not flow to lake catchment not a river gird  ",
            np.sum(BD_Out_Lakecat_Nriv_mask),
        )

        #####  Locate the grids that at the target grids ege
        Grid_Nee_Mo = np.argwhere(BD_Out_Lakecat_Nriv_mask == 1)
        #        print(print(Grid_Nee_Mo[:,0], Grid_Nee_Mo[:,1]))
        
        # problem_points is the point needs to modfiy flow direction 
        problem_points = copy.copy(Grid_Nee_Mo)        
        iter = 0

        # loop until all len problem_points is zero or iteration times exceed        
        while (iter < np.sum(BD_Out_Lakecat_mask) + 10 and len(problem_points) >0):

            # loop for each problem points, and check if it can flow to an lake 
            # catchment, by chaning it's flow direction only.
            #
            iter = iter + 1
            idx = 0
            new_problem_points = np.full((len(problem_points)+10,2),-9999)
            for i in range(0, len(problem_points)):
                # get problem point row anc col
                trow = problem_points[i, 0]
                tcol = problem_points[i, 1]
                #                print(i,ipo,trow,tcol)
                if trow >= nrows - 1 or tcol == ncols - 1:
                    continue
                    
                # try to modify the flow direction of this row and col 
                # see if it can flow to the lake catchment 
                
                ndir, IS_Change, changed_ndir = changeflowdirectionofedgegrids(
                    ndir,
                    trow,
                    tcol,
                    Lakeincat_mask,
                    1,
                    ncols,
                    nrows,
                    BD_Out_Lakecat_Nriv_mask,
                    changed_ndir,
                )
                
                # if the point can flow to the lake catchment 
                # store the points in to goodpoint2 
                if IS_Change > 0:
                    Lakeincat_mask[trow,tcol] = 1
                else:
                    new_problem_points[idx,0] = trow
                    new_problem_points[idx,1] = tcol 
                    idx = idx + 1
            # keep only larger         
            new_problem_points=new_problem_points[new_problem_points[:,0] > 0]        
            problem_points = new_problem_points
        print(" # of uncorrect lake boundary grids     ", len(problem_points))
    outlakeids = outlakeids[outlakeids[:, 0] > 0]
    return outlakeids, changed_ndir, ndir, BD_problem


def changeflowdirectionofedgegrids(
    N_dir,
    p_row,
    p_col,
    lake1,
    lid,
    ncols,
    nrows,
    BD_Out_Lakecat_Nriv_mask,
    Changed_ndir,
):
    ndir = copy.copy(N_dir)
    changed_ndir = copy.copy(Changed_ndir)
    IS_Change = 0

    if (
        lake1[p_row + 0, p_col + 1] == lid
    ):  #### it is a broudary grids and did not flow to the lake outlet
        ndir[p_row, p_col] = 1  ### change the flow direction of new point to old points
        changed_ndir[p_row, p_col] = 1
        IS_Change = 1
    # dir 2
    if lake1[p_row + 1, p_col + 1] == lid:

        ndir[p_row, p_col] = 2
        changed_ndir[p_row, p_col] = 2
        IS_Change = 1
    # dir 3
    if lake1[p_row + 1, p_col + 0] == lid:

        ndir[p_row, p_col] = 4
        changed_ndir[p_row, p_col] = 4
        IS_Change = 1
    # dir 4
    if lake1[p_row + 1, p_col - 1] == lid:

        ndir[p_row, p_col] = 8
        changed_ndir[p_row, p_col] = 8
        IS_Change = 1

    # dir 5
    if lake1[p_row + 0, p_col - 1] == lid:

        ndir[p_row, p_col] = 16
        changed_ndir[p_row, p_col] = 16
        IS_Change = 1

    # dir 6
    if lake1[p_row - 1, p_col - 1] == lid:

        ndir[p_row, p_col] = 32
        changed_ndir[p_row, p_col] = 32
        IS_Change = 1
    # dir 7
    if lake1[p_row - 1, p_col + 0] == lid:

        ndir[p_row, p_col] = 64
        changed_ndir[p_row, p_col] = 64
        IS_Change = 1
    # dir 8
    if lake1[p_row - 1, p_col + 1] == lid:

        ndir[p_row, p_col] = 128
        changed_ndir[p_row, p_col] = 128
        IS_Change = 1
    #    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir, IS_Change, changed_ndir


def Dirpoints_v3(
    N_dir,
    p_row,
    p_col,
    lake1,
    lid,
    goodpoint,
    k,
    ncols,
    nrows,
    BD_Out_Lakecat_Nriv_mask,
    Changed_ndir,
):
    ### this function change flow direction of some grids around p_row, p_col, flow to p_row, p_col

    ndir = copy.copy(N_dir)
    changed_ndir = copy.copy(Changed_ndir)
    ip = copy.copy(k) + 1
    # dir 1
    if (
        BD_Out_Lakecat_Nriv_mask[p_row + 0, p_col + 1] == lid
    ):  #### it is a broudary grids and did not flow to the lake outlet
        tt = goodpoint[
            goodpoint[:, 0] == p_row + 0,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col + 1,
                ]
            )
            < 1
        ):  #### the point not exist in good points which store all boundary points
            ndir[
                p_row + 0, p_col + 1
            ] = 16  ### change the flow direction of new point to old points
            changed_ndir[p_row + 0, p_col + 1] = 16
            goodpoint[ip, 0] = p_row + 0
            goodpoint[ip, 1] = p_col + 1
            ip = ip + 1
    # dir 2
    if BD_Out_Lakecat_Nriv_mask[p_row + 1, p_col + 1] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row + 1,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col + 1,
                ]
            )
            < 1
        ):
            ndir[p_row + 1, p_col + 1] = 32
            changed_ndir[p_row + 1, p_col + 1] = 32
            goodpoint[ip, 0] = p_row + 1
            goodpoint[ip, 1] = p_col + 1
            ip = ip + 1
    # dir 3
    if BD_Out_Lakecat_Nriv_mask[p_row + 1, p_col + 0] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row + 1,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col + 0,
                ]
            )
            < 1
        ):
            ndir[p_row + 1, p_col + 0] = 64
            changed_ndir[p_row + 1, p_col + 0] = 64
            goodpoint[ip, 0] = p_row + 1
            goodpoint[ip, 1] = p_col + 0
            ip = ip + 1
    # dir 4
    if BD_Out_Lakecat_Nriv_mask[p_row + 1, p_col - 1] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row + 1,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col - 1,
                ]
            )
            < 1
        ):
            ndir[p_row + 1, p_col - 1] = 128
            changed_ndir[p_row + 1, p_col - 1] = 128
            goodpoint[ip, 0] = p_row + 1
            goodpoint[ip, 1] = p_col - 1
            ip = ip + 1
    # dir 5
    if BD_Out_Lakecat_Nriv_mask[p_row + 0, p_col - 1] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row + 0,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col - 1,
                ]
            )
            < 1
        ):
            ndir[p_row + 0, p_col - 1] = 1
            changed_ndir[p_row + 0, p_col - 1] = 1
            goodpoint[ip, 0] = p_row + 0
            goodpoint[ip, 1] = p_col - 1
            ip = ip + 1
    # dir 6
    if BD_Out_Lakecat_Nriv_mask[p_row - 1, p_col - 1] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row - 1,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col - 1,
                ]
            )
            < 1
        ):
            ndir[p_row - 1, p_col - 1] = 2
            changed_ndir[p_row - 1, p_col - 1] = 2
            goodpoint[ip, 0] = p_row - 1
            goodpoint[ip, 1] = p_col - 1
            ip = ip + 1
    # dir 7
    if BD_Out_Lakecat_Nriv_mask[p_row - 1, p_col + 0] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row - 1,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col + 0,
                ]
            )
            < 1
        ):
            ndir[p_row - 1, p_col + 0] = 4
            changed_ndir[p_row - 1, p_col + 0] = 4
            goodpoint[ip, 0] = p_row - 1
            goodpoint[ip, 1] = p_col + 0
            ip = ip + 1
    # dir 8
    if BD_Out_Lakecat_Nriv_mask[p_row - 1, p_col + 1] == lid:

        tt = goodpoint[
            goodpoint[:, 0] == p_row - 1,
        ]
        if (
            len(
                tt[
                    tt[:, 1] == p_col + 1,
                ]
            )
            < 1
        ):
            ndir[p_row - 1, p_col + 1] = 8
            changed_ndir[p_row - 1, p_col + 1] = 8
            goodpoint[ip, 0] = p_row - 1
            goodpoint[ip, 1] = p_col + 1
            ip = ip + 1
    #    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir, goodpoint, ip, changed_ndir
