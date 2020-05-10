
from qgis.core import *
import qgis
from qgis.analysis import QgsNativeAlgorithms
from simpledbf import Dbf5
import os
import sys
import numpy as np
import shutil
from shutil import copyfile
from distutils.dir_util import copy_tree
import tempfile
import copy
import pandas as pd
import sqlite3
from GetBasinoutlet import Getbasinoutlet,Nextcell
from Generatecatinfo import Generatecatinfo,Generatecatinfo_riv,calculateChannaln,Writecatinfotodbf,Streamorderanddrainagearea
from WriteRavenInputs import writelake,Writervhchanl
from WriteRavenInputs import WriteObsfiles
from RavenOutputFuctions import plotGuagelineobs







def checklakeboundary(lake1,lid,p_row,p_col):
    numnonlake = 0
    nonlakedir = np.full(10,-999)
    if lake1[p_row + 0,p_col + 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 0
    if lake1[p_row + 1,p_col + 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 1
    if lake1[p_row + 1,p_col + 0] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 2
    if lake1[p_row + 1,p_col - 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 3
    if lake1[p_row + 0,p_col - 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 4
    if lake1[p_row - 1,p_col - 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 5
    if lake1[p_row - 1,p_col + 0] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 6
    if lake1[p_row - 1,p_col + 1] != lid:
        numnonlake = numnonlake + 1
        nonlakedir[numnonlake] = 7
    return numnonlake,nonlakedir
    


########## Return row and column number of a specific catchment's outlet  
# hylake       : a dataframe from HydroLAKE database containing all infromation about lakes 
# noncnlake    : a two dimension array of lakes that are not connected by river network 
# NonConLThres : the flow accumulation 2D array  
# dir    : the flow direction 2D array
# nrows, ncols  : maximum number of rows and columns in the basin 2D array   
# return crow,ccol  : row and column id of the catchment outlet.

def Dirpoints2(N_dir,p_row,p_col,lake1,lid,goodpoint,k,ncols,nrows):
    ndir = copy.copy(N_dir)
    ip = copy.copy(k) + 1
#dir 1
    if lake1[p_row + 0,p_col + 1] == lid:  #### it is a lake cell
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 0,p_col + 1)
        if numnonlake != 0: # it is a boundary lake cell
            tt = goodpoint[goodpoint[:,0] == p_row + 0,]
            if len(tt[tt[:,1] == p_col + 1,]) < 1:#### the point not exist in good points which store all boundary points
                ndir[p_row + 0,p_col + 1] = 16   ### change the flow direction of new point to old points
                goodpoint[ip,0] = p_row + 0
                goodpoint[ip,1] = p_col + 1
                ip = ip + 1
#dir 2
    if lake1[p_row + 1,p_col + 1] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 1,p_col + 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 1,]
            if len(tt[tt[:,1] == p_col + 1,]) < 1:
                ndir[p_row + 1,p_col + 1] = 32
                goodpoint[ip,0] = p_row + 1
                goodpoint[ip,1] = p_col + 1
                ip = ip + 1
#dir 3
    if lake1[p_row + 1,p_col + 0] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 1,p_col + 0)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 1,]
            if len(tt[tt[:,1] == p_col + 0,]) < 1:
                ndir[p_row + 1,p_col + 0] = 64
                goodpoint[ip,0] = p_row + 1
                goodpoint[ip,1] = p_col + 0
                ip = ip + 1
#dir 4
    if lake1[p_row + 1,p_col - 1] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 1,p_col - 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 1,]
            if len(tt[tt[:,1] == p_col - 1,]) < 1:
                ndir[p_row + 1,p_col - 1] = 128
                goodpoint[ip,0] = p_row + 1
                goodpoint[ip,1] = p_col - 1
                ip = ip + 1
#dir 5
    if lake1[p_row + 0,p_col - 1] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:   ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row + 0,p_col - 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row + 0,]
            if len(tt[tt[:,1] == p_col - 1,]) < 1:
                ndir[p_row + 0,p_col - 1] = 1
                goodpoint[ip,0] = p_row + 0
                goodpoint[ip,1] = p_col - 1
                ip = ip + 1
#dir 6
    if lake1[p_row - 1,p_col - 1] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:  ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row - 1,p_col - 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row - 1,]
            if len(tt[tt[:,1] == p_col - 1,]) < 1:
                ndir[p_row - 1,p_col - 1] = 2
                goodpoint[ip,0] = p_row - 1
                goodpoint[ip,1] = p_col - 1
                ip = ip + 1
#dir 7
    if lake1[p_row - 1,p_col + 0] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2: ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row - 1,p_col + 0)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row - 1,]
            if len(tt[tt[:,1] == p_col + 0,]) < 1:
                ndir[p_row - 1,p_col + 0] = 4
                goodpoint[ip,0] = p_row - 1
                goodpoint[ip,1] = p_col + 0
                ip = ip + 1
#dir 8
    if lake1[p_row - 1,p_col + 1] == lid:
        if p_row + 1 >= nrows-3 or p_col + 1 >= ncols - 3 or p_row - 1 <= 2 or p_col - 1 <= 2:  ### lake not at the boundary of domain
            numnonlake = 99
        else:
            numnonlake,nonlakedir = checklakeboundary(lake1,lid,p_row - 1,p_col + 1)
        if numnonlake != 0:
            tt = goodpoint[goodpoint[:,0] == p_row - 1,]
            if len(tt[tt[:,1] == p_col + 1,]) < 1:
                ndir[p_row - 1,p_col + 1] = 8
                goodpoint[ip,0] = p_row - 1
                goodpoint[ip,1] = p_col + 1
                ip = ip + 1
#    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir,goodpoint,ip

def ChangeDIR(dir,lake1,acc,ncols,nrows,outlakeids,nlakegrids,cat3):
    ndir = copy.copy(dir)
    for i in range(0,len(outlakeids)):
        lid = outlakeids[i,0]
        goodpoint = np.full((20000,2),-99999)
        lrowcol = np.argwhere(lake1==lid).astype(int)
        
        if len(lrowcol) > nlakegrids and outlakeids[i,1] > 0.9:   ### smaller than nlakegrids or smaller than 0.9
            continue
        if outlakeids[i,1] > 0.97: ### smaller than 0.97
            continue
            
        prow,pcol = Getbasinoutlet(lid,lake1,acc,dir,nrows,ncols)

        goodpoint[0,0] = prow
        goodpoint[0,1] = pcol
        if prow >= nrows - 1 or pcol == ncols  - 1:
            continue
        ip = 0
        k = 0
        ipo = -1
        while ip > ipo:
            for i in range(0,len(goodpoint[goodpoint[:,0]>0,])):
                if i > ipo:
                    trow = goodpoint[i,0]
                    tcol = goodpoint[i,1]
                    if trow >= nrows - 1 or tcol == ncols  - 1:
                        continue
                    ndir,goodpoint,k1= Dirpoints2(ndir,trow,tcol,lake1,lid,goodpoint,k,ncols,nrows)
                    k = k1 - 1
            ipo = ip
            ip = len(goodpoint[goodpoint[:,0]>0,]) - 1
            k = ip
    return ndir

  
def Checklake(prow,pcol,nrows,ncols,lid,lake):
    noout=0
    ### double check if the head stream cell is nearby the lake, if near by the lake the stream was ignored
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if lake[prow-1,pcol+1] == lid:
            noout=1
        if lake[prow-1,pcol-1] == lid:
            noout=1
        if lake[prow-1,pcol] == lid:
            noout=1
        if lake[prow,pcol+1] == lid:
            noout=1
        if lake[prow,pcol-1] == lid:
            noout=1
        if lake[prow+1,pcol-1] == lid:
            noout=1
        if lake[prow+1,pcol+1] == lid:
            noout=1
        if lake[prow+1,pcol] == lid:
            noout=1
    return noout
############################################################################33

####################################################################3
def Checkcat(prow,pcol,nrows,ncols,lid,lake):
    noout=0
    ### double check if the head stream cell is nearby the lake, if near by the lake the stream was ignored
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if not lake[prow-1,pcol+1] == lid:
            noout=1 + noout
        if not lake[prow-1,pcol-1] == lid:
            noout=1 + noout
        if not lake[prow-1,pcol] == lid:
            noout=1 + noout
        if not lake[prow,pcol+1] == lid:
            noout=1 + noout
        if not lake[prow,pcol-1] == lid:
            noout=1 + noout
        if not lake[prow+1,pcol-1] == lid:
            noout=1 + noout
        if not lake[prow+1,pcol+1] == lid:
            noout=1 + noout
        if not lake[prow+1,pcol] == lid:
            noout=1 + noout
    return noout

###################################################################3

def selectlake(hylake,noncnlake,NonConLThres,hylakeinfo):
#    arcpy.AddMessage("123asdfasdfasdfasd")
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(noncnlake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(noncnlake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['Hylak_id'] == sl_lid]
        if len(slakeinfo) <=0:
            continue
        if slakeinfo.iloc[0]['Lake_area'] >= NonConLThres:
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = sl_lid
    return sl_lake

def selectlake2(hylake,Lakehres,hylakeinfo):
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(sl_lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(sl_lake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['Hylak_id'] == sl_lid]
        if len(slakeinfo)<=0:
#            print("Lake excluded     asdfasd " + str(sl_lid))
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999
            continue
        if slakeinfo.iloc[0]['Lake_area'] < Lakehres:
#            print("Lake excluded     due to area " + str(sl_lid))
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999
    return sl_lake

##################################################################3

def check_lakecatchment(cat3,lake,fac,fdir,bsid,nrows,ncols):
    cat = copy.copy(cat3)
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    outlakeids = np.full((1000000,2),-99999.999)
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        lakacc = np.full((len(lrowcol),3),-9999)
        lakacc[:,0] = lrowcol[:,0]
        lakacc[:,1] = lrowcol[:,1]
        lakacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakacc = lakacc[lakacc[:,2].argsort()]
        lorow = lakacc[len(lakacc)-1,0]
        locol = lakacc[len(lakacc)-1,1]  ###### lake outlet row and col
        arclakeid = cat[lorow,locol]  ####### lake catchment id 
        lakecatrowcol = np.argwhere(cat==arclakeid).astype(int)    #### lake catchment cells 
        Lakeincat1 = lake[lakecatrowcol[:,0],lakecatrowcol[:,1]]
        nlake = np.argwhere(Lakeincat1==lakeid).astype(int)
        outlakeids[i,0] = lakeid
        outlakeids[i,1] = float(len(nlake)/len(lrowcol)) 
#        print(lakeid,arclakeid,len(nlake),len(lrowcol),float(len(nlake)/len(lrowcol))) 
    outlakeids= outlakeids[outlakeids[:,0] > 0]
    return outlakeids
######################################################################3
def Addobspoints(obs,pourpoints,boid,cat):
    obsids = np.unique(obs)
    obsids = obsids[obsids>=0]
    for i in range(0,len(obsids)):
        rowcol = np.argwhere(obs==obsids[i]).astype(int)
        if cat[rowcol[0,0],rowcol[0,1]] > 0:
            if pourpoints[rowcol[0,0],rowcol[0,1]] < 0:
                pourpoints[rowcol[0,0],rowcol[0,1]] = boid + obsids[i]
    return pourpoints

####################################################################
def CE_mcat4lake(cat1,lake,fac,fdir,bsid,nrows,ncols,Pourpoints):
    #####adjust for some lakes was divided into two catchment beacuse of flow direction and in stream. double check each lake
    ##### and decide if need to merge these catchment into one lake catchment.
    cat = copy.copy(cat1)
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    outlakeids = np.full((1000000,2),-99999.999)
    outi = 0
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        lakacc = np.full((len(lrowcol),3),-9999)
        lakacc[:,0] = lrowcol[:,0]
        lakacc[:,1] = lrowcol[:,1]
        lakacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakacc = lakacc[lakacc[:,2].argsort()]
        lorow = lakacc[len(lakacc)-1,0]
        locol = lakacc[len(lakacc)-1,1]  ###### lake outlet row and col
        arclakeid = cat[lorow,locol]  ####### lake catchment id 
        lakecatrowcol = np.argwhere(cat==arclakeid).astype(int)    #### lake catchment cells 
        Lakeincat1 = lake[lakecatrowcol[:,0],lakecatrowcol[:,1]]
        nlake = np.argwhere(Lakeincat1==lakeid).astype(int) 

        lakenrow,lakencol = Nextcell(fdir,lorow,locol)
        if lakenrow < 0 or lakencol < 0:
            lakedowcatid = -999
        else:
            if lakenrow >= nrows or lakencol >= ncols:
                continue
            lakedowcatid = cat[lakenrow,lakencol]
        # if not arclakeid < bsid and arclakeid > blid:
        #     continue
        arcatid,catcounts = np.unique(cat[lrowcol[:,0],lrowcol[:,1]],return_counts=True) ###### all catchment id containing this lake
        tarid = 0
        ### if there are more than 1 catchment in cat1, determine if they need to be combined
        ### check if these catchment flow into the lake if it is true, change catchment id into lake catchment id
        if len(arcatid)>1:  #
#            if float(len(lakecatrowcol))/float(len(lrowcol)) < 0.9: #and len(lrowcol) < 10000: #: float(max(catcounts))/float(len(lrowcol)) < 0.8 and 
            outlakeids[outi,0] = lakeid
            outlakeids[outi,1] = float(len(nlake))/float(len(lrowcol))
#            print(outlakeids[outi,0],outlakeids[outi,1])
            outi = outi + 1
            for j in range(0,len(arcatid)):
                crowcol = np.argwhere(cat==arcatid[j]).astype(int)
                catacc = np.full((len(crowcol),3),-9999)
                catacc[:,0] = crowcol[:,0]
                catacc[:,1] = crowcol[:,1]
                catacc[:,2] = fac[crowcol[:,0],crowcol[:,1]]
                catacc = catacc[catacc[:,2].argsort()]
                catorow = catacc[len(catacc)-1,0]
                catocol = catacc[len(catacc)-1,1] ### catchment outlet
                Lakeincat = lake[crowcol[:,0],crowcol[:,1]]
                nlake = np.argwhere(Lakeincat==lakeid).astype(int)
                nrow,ncol = Nextcell(fdir,catorow,catocol) #####Get the next row and col of downstream catchment
                if nrow < 0 or ncol < 0:
                    continue
                if nrow < nrows and ncol < ncols:
               ### if downstream catchment is target lake,and this catchment is an lakeinflow catchment combine them
                    if cat[nrow,ncol] == arclakeid and float(len(nlake))/float(len(crowcol)) > 0.1 and cat[catorow,catocol] > bsid:
                        cat[crowcol[:,0],crowcol[:,1]] = arclakeid
                    if float(len(nlake))/float(len(lrowcol)) > 0.1 and cat[catorow,catocol] > bsid and cat[catorow,catocol] != lakedowcatid:
#                        arcpy.AddMessage("2")
                        cat[crowcol[:,0],crowcol[:,1]] = arclakeid
#                        if cat[catorow,catocol] != arclakeid and cat[nrow,ncol] != arclakeid:
#                            print lakeid
                    if cat[nrow,ncol] > bsid and arcatid[j] > bsid:  #### lake input cat route to another lake input catch
                        cat[crowcol[:,0],crowcol[:,1]] = cat[nrow,ncol]
        pp = Pourpoints[lrowcol[:,0],lrowcol[:,1]]
        pp = np.unique(pp)
        pp = pp[pp > 0]
        if len(pp) == 1:
            cat[lrowcol[:,0],lrowcol[:,1]] = arclakeid
    outlakeids= outlakeids[outlakeids[:,0] > 0]
    return cat,outlakeids
###################################################33
def CE_mcat4lake2(cat1,lake,fac,fdir,bsid,nrows,ncols,Pourpoints,noncnlake):
    cat = copy.copy(cat1)
    Non_con_lake_cat = copy.copy(cat1)
    Non_con_lake_cat[:,:] = -9999
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    noncnlake = np.unique(noncnlake)
    noncnlake = noncnlake[noncnlake>=0]
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        lakacc = np.full((len(lrowcol),3),-9999)
        lakacc[:,0] = lrowcol[:,0]
        lakacc[:,1] = lrowcol[:,1]
        lakacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakacc = lakacc[lakacc[:,2].argsort()]
        lorow = lakacc[len(lakacc)-1,0]
        locol = lakacc[len(lakacc)-1,1]  ###### lake outlet row and col
        arclakeid = cat1[lorow,locol]
        pp = Pourpoints[lorow,locol]
        pp = np.unique(pp)
        pp = pp[pp > 0]
#        print(lakeid,len(np.argwhere(noncnlake==lakeid)),pp,Pourpoints[lorow,locol],arclakeid)
        if len(pp) == 1:
            if arclakeid < 0:
                cat[lrowcol[:,0],lrowcol[:,1]] = pp
            else:
                cat[lrowcol[:,0],lrowcol[:,1]] = arclakeid
            if len(np.argwhere(noncnlake==lakeid)) > 0: ##if lake belong to non connected lakes
                if arclakeid < 0:
                    nonlrowcol = np.argwhere(cat==pp).astype(int)
                    Non_con_lake_cat[nonlrowcol[:,0],nonlrowcol[:,1]] = lakeid
                else:
                    nonlrowcol = np.argwhere(cat==arclakeid).astype(int)
                    Non_con_lake_cat[nonlrowcol[:,0],nonlrowcol[:,1]] = lakeid
#### For some reason pour point was missing in non-contribute catchment 
    Pours = np.unique(Pourpoints)
    Pours = Pours[Pours>0]
    for i in range(0,len(Pours)):
        pourid = Pours[i]
        rowcol = Pourpoints == pourid
        if cat[rowcol] < 0:
            nout = Checkcat(rowcol[0,0],rowcol[0,1],nrows,ncols,pourid,cat)
            if len(cat[cat == pourid]) > 0 and nout < 8:
                cat[rowcol] = pourid
    rowcol1 = fac > 0
    rowcol2 = cat < 0
    noncontribuite = np.logical_and(rowcol1, rowcol2)
    cat[noncontribuite] = 2*max(np.unique(cat)) + 1
    return cat,Non_con_lake_cat
######################################################
def CE_Lakeerror(fac,fdir,lake,cat2,bsid,blid,boid,nrows,ncols,cat):
    Watseds = copy.copy(cat2)
    Poups = np.unique(Watseds)
    Poups = Poups[Poups>=0]
    ##### Part 2, remove some small catchment which is not lake catchment
    out = np.full((len(Poups),4),-9999)
    for i in range(0,len(Poups)):
        catid = Poups[i]
        if catid > boid:
            continue #### do nothing for observation catchments
        rowcol = np.argwhere(Watseds==catid).astype(int)
        catacc = np.full((len(rowcol),3),-9999)
        catacc[:,0] = rowcol[:,0]
        catacc[:,1] = rowcol[:,1]
        catacc[:,2] = fac[rowcol[:,0],rowcol[:,1]]
        catacc = catacc[catacc[:,2].argsort()].astype(int)
        rowcol[0,0] = catacc[len(catacc)-1,0]
        rowcol[0,1] = catacc[len(catacc)-1,1]
        nrow,ncol = Nextcell(fdir,rowcol[0,0],rowcol[0,1])### get the downstream catchment id
        if nrow < 0 or ncol < 0:
            continue
        if nrow < nrows and ncol < ncols:
            if len(rowcol) < 10 and Watseds[rowcol[0,0],rowcol[0,1]] > bsid:
                Watseds[catacc[:,0],catacc[:,1]] = Watseds[nrow,ncol]
            if len(rowcol) < 10 and Watseds[rowcol[0,0],rowcol[0,1]] < blid:
                Watseds[catacc[:,0],catacc[:,1]] = Watseds[nrow,ncol]
    return Watseds

#########################################33
def GenerateFinalPourpoints(fac,fdir,lake,cat3,bsid,blid,boid,nrows,ncols,cat,obs):
    Poups = copy.copy(cat3)
    Poups[:,:]=-9999
    GWat = copy.copy(cat3)
    GWatids = np.unique(cat3)
    GWatids = GWatids[GWatids>=0]
    ncatid = 1
    for i in range(0,len(GWatids)):
        trow,tcol = Getbasinoutlet(GWatids[i],GWat,fac,fdir,nrows,ncols)
#        arcpy.AddMessage("result     " +str(trow) +"    " +str(tcol))
        Poups[trow,tcol] = ncatid
        ncatid = ncatid + 1
    OWat = copy.copy(cat)
    OWatids = np.unique(cat)
    OWatids = OWatids[OWatids>=0]
    for i in range(0,len(OWatids)):
        trow,tcol = Getbasinoutlet(OWatids[i],OWat,fac,fdir,nrows,ncols)
        if not GWat[trow,tcol] >= blid:
            if Poups[trow,tcol] < 0:
                Poups[trow,tcol] = ncatid
                ncatid = ncatid + 1
    obsids = np.unique(obs)
    obsids = obsids[obsids>0]
    for i in range(0,len(obsids)):
        rowcol = np.argwhere(obs==obsids[i]).astype(int)
        if Poups[rowcol[0,0],rowcol[0,1]] < 0 and lake[rowcol[0,0],rowcol[0,1]] < 0:
            Poups[rowcol[0,0],rowcol[0,1]] = ncatid
            ncatid = ncatid + 1
    return Poups
#######
####
# ####################################################33


def GenerPourpoint(cat,lake,Str,nrows,ncols,blid,bsid,bcid,fac,hydir):
    GP_cat = copy.copy(cat)
    sblid = copy.copy(blid)
    ############### Part 1 Get all pourpoints of hydroshed catchment
    arcatid = np.unique(cat)#### cat all catchment idd
    arcatid = arcatid[arcatid>=0]
    catoutloc = np.full((len(arcatid),3),-9999)
    for i in range(0,len(arcatid)):
        catid = arcatid[i]
        catrowcol = np.argwhere(cat==catid).astype(int)
        trow,tcol = Getbasinoutlet(catid,cat,fac,hydir,nrows,ncols)
        GP_cat[catrowcol[:,0],catrowcol[:,1]]=-9999   ### set the catment cells into null
        GP_cat[trow,tcol]=bcid #### change the outlet of catchment into wid
        bcid = bcid + 1
        catoutloc[i,0] = catid  ## catchment id
        catoutloc[i,1] = trow  #### catchment pourpont row
        catoutloc[i,2] = tcol  #### catchment pourpont col
#    writeraster(outFolder+subid+"_Pourpoints_1.asc",GP_cat)

    ##################Part 2 Get pourpoints of Lake inflow streams
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)): #### loop for each lake
        lid = arlakeid[i]       ##### obtain lake id 
        rowcol = np.argwhere(lake==lid.astype(int))  ### got row anc col of lake grids 
        nrow = rowcol.shape[0]  ### number of grids of the lake 
        Stridinlake = np.full(nrow,-9999)  #### 
        Stridinlake[:] = Str[rowcol[:,0],rowcol[:,1]]  ### get all stream with in the lake domain 
        Strid_L  = np.unique(Stridinlake[np.argwhere(Stridinlake > 0).astype(int)]) ### Get unique stream id of sach stream 
        ##### find the intercept point of stream and lake
        for j in range(0,len(Strid_L)):  #### loop for each stream intercept with lake
            strid = Strid_L[j]
            strrowcol = np.argwhere(Str == strid).astype(int)
            nstrrow = strrowcol.shape[0]
            Strchek = np.full((nstrrow,4),-9999)##### 0 row, 1 col, 2 fac,
            Strchek[:,0] = strrowcol[:,0]
            Strchek[:,1] = strrowcol[:,1]
            Strchek[:,2] = fac[strrowcol[:,0],strrowcol[:,1]]
            Strchek[:,3] = lake[strrowcol[:,0],strrowcol[:,1]]
            Strchek = Strchek[Strchek[:,2].argsort()].astype(int)
            
            ###3 find the first intersection point between river and lake 
            Grids_Lake_river = np.where(Strchek[:,3] == lid)
            irowst = np.nanmin(Grids_Lake_river)
            Intersection_Grid_row = Strchek[irowst,0]
            Intersection_Grid_col = Strchek[irowst,1]

            noout = 0
            ### double check if the head stream cell is nearby the lake
            if Strchek[0,0] != 0 and Strchek[0,0] != nrows -1 and Strchek[0,1] != 0 and Strchek[0,1] != ncols-1:
                noout = Checklake(Strchek[0,0],Strchek[0,1],nrows,ncols,lid,lake)            
            ####
            
            if irowst != 0: ## is the stream is not start within the lake 
                if lake[Strchek[irowst-1,0],Strchek[irowst-1,1]] == -9999:  ### double check 
                    ##### this means the stream connect two lakes, so must assign an pourpoints
                    if len(np.unique(lake[Strchek[0:irowst,0],Strchek[0:irowst,1]])) >= 3:
                        GP_cat[Strchek[irowst-1,0],Strchek[irowst-1,1]] = bsid
                        bsid = bsid + 1

                    ##### the head stream celll is not nearby the lake
                    if noout == 0:
                        GP_cat[Strchek[irowst-1,0],Strchek[irowst-1,1]] = bsid
                        bsid = bsid + 1
            #### it is possible that two steam combine together near the lake, double check if the stream conncet to
            # anotehr stream and this steam is not witin the lake
            if irowst == 0 or noout == 1:
                nostr = Str[Strchek[0,0],Strchek[0,1]]
                a = 0
                orowcol = np.full((8,3),-9999)
                if Strchek[0,0] != 0 and Strchek[0,0] != nrows -1 and Strchek[0,1] != 0 and Strchek[0,1] != ncols-1:
                    if Str[Strchek[0,0]-1,Strchek[0,1]+1] != -9999 and lake[Strchek[0,0]-1,Strchek[0,1]+1] == -9999:
                        orowcol[a,0] = Strchek[0,0]-1
                        orowcol[a,1] = Strchek[0,1]+1
                        orowcol[a,2] = Str[Strchek[0,0]-1,Strchek[0,1]+1]
                        a = a+1
                    if Str[Strchek[0,0]-1,Strchek[0,1]-1] != -9999 and lake[Strchek[0,0]-1,Strchek[0,1]-1] == -9999:
                        orowcol[a,0] = Strchek[0,0]-1
                        orowcol[a,1] = Strchek[0,1]-1
                        orowcol[a,2] = Str[Strchek[0,0]-1,Strchek[0,1]-1]
                        a = a+1
                    if Str[Strchek[0,0]-1,Strchek[0,1]] != -9999 and lake[Strchek[0,0]-1,Strchek[0,1]] == -9999:
                        orowcol[a,0] = Strchek[0,0]-1
                        orowcol[a,1] = Strchek[0,1]-0
                        orowcol[a,2] = Str[Strchek[0,0]-1,Strchek[0,1]-0]
                        a = a+1
                    if Str[Strchek[0,0],Strchek[0,1]+1] != -9999 and lake[Strchek[0,0],Strchek[0,1]+1] == -9999:
                        orowcol[a,0] = Strchek[0,0]-0
                        orowcol[a,1] = Strchek[0,1]+1
                        orowcol[a,2] = Str[Strchek[0,0]-0,Strchek[0,1]+1]
                        a = a+1
                    if Str[Strchek[0,0],Strchek[0,1]-1] != -9999 and lake[Strchek[0,0],Strchek[0,1]-1] == -9999:
                        orowcol[a,0] = Strchek[0,0]-0
                        orowcol[a,1] = Strchek[0,1]-1
                        orowcol[a,2] = Str[Strchek[0,0]-0,Strchek[0,1]-1]
                        a = a+1
                    if Str[Strchek[0,0]+1,Strchek[0,1]-1] != -9999 and lake[Strchek[0,0]+1,Strchek[0,1]-1] == -9999:
                        orowcol[a,0] = Strchek[0,0]+1
                        orowcol[a,1] = Strchek[0,1]-1
                        orowcol[a,2] = Str[Strchek[0,0]+1,Strchek[0,1]-1]
                        a = a+1
                    if Str[Strchek[0,0]+1,Strchek[0,1]+1] != -9999 and lake[Strchek[0,0]+1,Strchek[0,1]+1] == -9999:
                        orowcol[a,0] = Strchek[0,0]+1
                        orowcol[a,1] = Strchek[0,1]+1
                        orowcol[a,2] = Str[Strchek[0,0]+1,Strchek[0,1]+1]
                        a = a+1
                    if Str[Strchek[0,0]+1,Strchek[0,1]] != -9999 and lake[Strchek[0,0]+1,Strchek[0,1]] == -9999:
                        orowcol[a,0] = Strchek[0,0]+1
                        orowcol[a,1] = Strchek[0,1]-0
                        orowcol[a,2] = Str[Strchek[0,0]+1,Strchek[0,1]-0]
                        a = a+1
                if a > 0:
                    for ka in range(0,a):
                        nostr =orowcol[ka,2]  ## up stream stram id 
                        srowcol = np.argwhere(Str==nostr).astype(int)
                        snrow = srowcol.shape[0]
                        iStrchek = np.full((snrow,4),-9999)##### 0 row, 1 col, 2 fac,
                        iStrchek[:,0] = srowcol[:,0]
                        iStrchek[:,1] = srowcol[:,1]
                        iStrchek[:,2] = fac[srowcol[:,0],srowcol[:,1]]
                        iStrchek = iStrchek[iStrchek[:,2].argsort()]
                        noout = Checklake(iStrchek[0,0],iStrchek[0,1],nrows,ncols,lid,lake)
                        Lakinstr = np.full(snrow,-9999)
                        Lakinstr[:] = lake[srowcol[:,0],srowcol[:,1]]
                        d = np.argwhere(Lakinstr==lid).astype(int)  #### the connected stream should not within the lake
                        if len(d) < 1 and noout == 0:
                            GP_cat[orowcol[ka,0],orowcol[ka,1]] = bsid
                            bsid = bsid + 1
################################################################################
################## Part 3Get Lake pourpoint id and remove cat pourpoint that contribute to lake
#    writeraster(outFolder+"Pourpoints_2.asc",GP_cat)
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        lrowcol = np.argwhere(lake==lakeid).astype(int)
        arcatid = np.unique(cat[lrowcol[:,0],lrowcol[:,1]]) ## Get all catchment that intercept with lake
        lakeacc = np.full((len(lrowcol),3),-9999)
        lakeacc[:,0] = lrowcol[:,0]
        lakeacc[:,1] = lrowcol[:,1]
        lakeacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
        lakeacc = lakeacc[lakeacc[:,2].argsort()]
        maxcatid = cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]] ### Get the catchment id that will keeped
        
        ### remove the all pourpoints close to lake pourpoints
        if lakeacc[len(lakeacc)-1,0] != 0 and lakeacc[len(lakeacc)-1,0] != nrows -1 and lakeacc[len(lakeacc)-1,1] != 0 and lakeacc[len(lakeacc)-1,1] != ncols-1:
            GP_cat[lakeacc[len(lakeacc)-1,0]-1,lakeacc[len(lakeacc)-1,1]-1]=-9999 
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]-1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]-1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]+1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]+1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]-1,lakeacc[len(lakeacc)-1,1]]=-9999
            
        ### remove pourpoints that not equal to maxcatid 
        for j in range(0,len(arcatid)):
            if arcatid[j] != maxcatid:
                crowcol = np.argwhere(cat==arcatid[j]).astype(int)
                checkcat = np.full((len(crowcol),4),-9999)
                checkcat[:,0] = crowcol[:,0]
                checkcat[:,1] = crowcol[:,1]
                checkcat[:,2] = GP_cat[crowcol[:,0],crowcol[:,1]]
                checkcat[:,3] = fac[crowcol[:,0],crowcol[:,1]]
                checkcat = checkcat[checkcat[:,3].argsort()]
                dele = checkcat[checkcat[:,2]<blid].astype(int)   #### do not delete the lake and strem ids
                GP_cat[dele[:,0],dele[:,1]]=-9999
                
                ### add keep catchment in case the catchment outlet is at the boundary of the domain
                nrow,ncol = Nextcell(hydir,checkcat[len(checkcat)-1,0],checkcat[len(checkcat)-1,1])
                if nrow > 0 or ncol >0:
                    if nrow >= nrows or ncols >= ncols:
                        continue
                    if cat[nrow,ncol] < 0:
                        GP_cat[checkcat[len(checkcat)-1,0],checkcat[len(checkcat)-1,1]] = bcid
                        bcid = bcid + 1
                        
        #### add lake outlet 
        GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]]= sblid
        sblid = sblid + 1
    return GP_cat
###################################################################3



##################################################################3


###################################################################3

##### out has two column the first column has sub id , the second has down sub id 
def Defcat(out,outletid):
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((len(out),1),-999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((len(out),1),-999)
        poshdid = 0
        for i in range(0,len(otsheds)):
            Shedid[psid] = otsheds[i]
            psid = psid + 1
            irow = np.argwhere(rout[:,1]==otsheds[i]).astype(int)
            for j in range(0,len(irow)):
                noutshd[poshdid] = rout[irow[j],0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd>=0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid>=0]
    return Shedid
###########


##### out has two column the first column has sub id , the second has down sub id 

###########

    
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
               
    cbranch                  = catchmentinfo[catchmentinfo[sub_colnm].isin(Modify_subids)]
    tarinfo                  = catchmentinfo[catchmentinfo[sub_colnm] == subid]   ### define these subs attributes
#    print(subid,Modify_subids) 
    ### average river slope info 
    
    mainriv_merg_info = mainriv.loc[mainriv['SubId'].isin(Modify_subids)]
    idx = tarinfo.index[0]
    if len(mainriv_merg_info) > 0:
        tarinfo.loc[idx,'RivLength'] = np.sum(mainriv_merg_info['RivLength'].values)    
        tarinfo.loc[idx,'RivSlope']  = np.average(mainriv_merg_info['RivSlope'].values ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'FloodP_n']  = np.average(mainriv_merg_info['FloodP_n'].values ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'Q_Mean']    = np.average(mainriv_merg_info['Q_Mean'].values   ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'Ch_n']      = np.average(mainriv_merg_info['Ch_n'].values     ,weights = mainriv_merg_info['RivLength'].values)
        tarinfo.loc[idx,'BkfWidth']  = np.max(mainriv_merg_info['BkfWidth'].values)
        tarinfo.loc[idx,'BkfDepth']  = np.max(mainriv_merg_info['BkfDepth'].values)
        
    tarinfo.loc[idx,'BasArea']       = np.sum(cbranch['BasArea'].values)
    tarinfo.loc[idx,'NonLDArea']     = np.sum(cbranch['NonLDArea'].values)
    tarinfo.loc[idx,'BasSlope']      = np.average(cbranch['BasSlope'].values,  weights = cbranch['BasArea'].values)
    tarinfo.loc[idx,'MeanElev']      = np.average(cbranch['MeanElev'].values,  weights = cbranch['BasArea'].values)
    tarinfo.loc[idx,'BasAspect']     = np.average(cbranch['BasAspect'].values, weights = cbranch['BasArea'].values)
    
    tarinfo.loc[idx,'Max_DEM']       = np.max(cbranch['Max_DEM'].values)
    tarinfo.loc[idx,'Min_DEM']       = np.min(cbranch['Min_DEM'].values)
    
    if Islake == 1:   ## Meger subbasin covered by lakes, Keep lake outlet catchment  DA, stream order info 
        tarinfo.loc[idx,'RivLength'] = 0.0
    elif Islake <0:
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

    mask = mapoldnew_info['SubId'].isin(Modify_subids)
    ### the old downsub id of the dissolved polygon is stored in DowSubId
    for col in tarinfo.columns:
        if col == 'SubId':
            mapoldnew_info.loc[mask,'nsubid']     = tarinfo[col].values[0]
        else:
            mapoldnew_info.loc[mask, col]         = tarinfo[col].values[0]
    return mapoldnew_info

def Modify_Feature_info(Path_feagure,mapoldnew_info):
    sub_colnm = 'SubId'
    layer_cat=QgsVectorLayer(Path_feagure,"")
    Attri_Name = layer_cat.fields().names()     
    features = layer_cat.getFeatures()      
    with edit(layer_cat):
        
        for sf in features:
            Atti_Valu    = sf.attributes()
            sf_subid     = sf[sub_colnm]
            tarinfo      = mapoldnew_info[mapoldnew_info['Old_SubId'] == sf_subid]
            for icolnm in range(0,len(Attri_Name)):     ### copy infomaiton
                if  Attri_Name[icolnm] == 'Obs_NM' or Attri_Name[icolnm] == 'SRC_obs':
                    sf[Attri_Name[icolnm]] = str(tarinfo[Attri_Name[icolnm]].values[0])
                elif Attri_Name[icolnm] == 'cat':
                    continue
                else:    
                    sf[Attri_Name[icolnm]] = float(tarinfo[Attri_Name[icolnm]].values[0])
            layer_cat.updateFeature(sf)
    del layer_cat
    return

####### 
def UpdateNonConnectedLakeCatchmentinfo(Path_Non_ConnL_Cat,mapoldnew_info):
    layer_cat=QgsVectorLayer(Path_Non_ConnL_Cat,"")
    Attri_Name = layer_cat.fields().names()     
    features = layer_cat.getFeatures()      
    with edit(layer_cat):
        for sf in features:
            sf_ocatid_lake  = float(sf['SubId_riv'])
            tarinfo         = mapoldnew_info[mapoldnew_info['Old_SubId'] == sf_ocatid_lake]
            sf['SubId_riv'] = float(tarinfo['SubId'].values[0])
            layer_cat.updateFeature(sf)
    del layer_cat
    return
##########

##########
def UpdateNonConnectedLakeArea_In_Finalcatinfo(Path_Finalcatinfo,Non_ConnL_Cat_info):
    layer_cat=QgsVectorLayer(Path_Finalcatinfo,"")
    Attri_Name = layer_cat.fields().names()     
    features = layer_cat.getFeatures()
    Non_ConnL_Cat_info['SubId_riv'] = Non_ConnL_Cat_info['SubId_riv'].astype(float)
    Non_ConnL_Cat_info['Area_m']    = Non_ConnL_Cat_info['Area_m'].astype(float)      
    with edit(layer_cat):
        for sf in features:
            sf_subid        = float(sf['SubId'])
            tarinfo         = Non_ConnL_Cat_info[Non_ConnL_Cat_info['SubId_riv'] == sf_subid]
            if (len(tarinfo) == 0):
                sf['NonLDArea']      = float(0)
            else:
                total_non_conn_lake_area = 0.0
                for idx in tarinfo.index:
                    total_non_conn_lake_area = total_non_conn_lake_area + float(tarinfo['Area_m'].values[0])
                    
                sf['NonLDArea'] = float(total_non_conn_lake_area)
            layer_cat.updateFeature(sf)
    del layer_cat
    return
#########

def UpdateConnectedLakeArea_In_Finalcatinfo(Path_Finalcatinfo,Conn_Lake_Ids):
    layer_cat=QgsVectorLayer(Path_Finalcatinfo,"")
    Attri_Name = layer_cat.fields().names()     
    features = layer_cat.getFeatures()     
    with edit(layer_cat):
        for sf in features:
            sf_subid        = float(sf['HyLakeId'])
            
            if sf_subid in Conn_Lake_Ids:
                continue 
            sf['HyLakeId']      = float(-1.2345)
            sf['LakeVol']       = float(-1.2345)
            sf['LakeArea']      = float(-1.2345)
            sf['LakeDepth']     = float(-1.2345)
            sf['Laketype']      = float(-1.2345) 
            sf['IsLake']        = float(-1.2345)     
            layer_cat.updateFeature(sf)
    del layer_cat
    return
    
    

#########
def Add_centroid_to_feature(Path_feagure,centroidx_nm = '#',centroidy_nm='#'):
    layer_cat=QgsVectorLayer(Path_feagure,"")
    Attri_Name = layer_cat.fields().names()     
    features = layer_cat.getFeatures()      
    with edit(layer_cat):
        for sf in features:
            centroidxy = sf.geometry().centroid().asPoint()
            sf[centroidx_nm] = centroidxy[0]
            sf[centroidy_nm] = centroidxy[1]
            layer_cat.updateFeature(sf)
    del layer_cat
    return 
    
##########
def Selectfeatureattributes(processing,Input = '#',Output='#',Attri_NM = '#',Values = []):
    exp =Attri_NM + '  IN  (  ' +  str(int(Values[0]))      
    for i in range(1,len(Values)):
        exp = exp + " , "+str(int(Values[i]))        
    exp = exp + ')'
    processing.run("native:extractbyexpression", {'INPUT':Input,'EXPRESSION':exp,'OUTPUT':Output}) 
       
#####
def UpdateTopology(mapoldnew_info,UpdateStreamorder = 1,UpdateSubId = 1):
    idx = mapoldnew_info.index
    
    if UpdateSubId > 0: 
        for i in range(0,len(idx)):
            nsubid     = mapoldnew_info.loc[idx[i],'nsubid']
            subid      = mapoldnew_info.loc[idx[i],'SubId']
            odownsubid = mapoldnew_info.loc[idx[i],'DowSubId']
                
            donsubidinfo = mapoldnew_info.loc[mapoldnew_info['SubId'] == odownsubid]
                
            if (len(donsubidinfo) >0):
                mapoldnew_info.loc[idx[i],'ndownsubid'] = donsubidinfo['nsubid'].values[0]
            else:
                mapoldnew_info.loc[idx[i],'ndownsubid'] = -1
            
        mapoldnew_info['Old_SubId']    = mapoldnew_info['SubId']
        mapoldnew_info['Old_DowSubId'] = mapoldnew_info['DowSubId']
        mapoldnew_info['SubId']        = mapoldnew_info['nsubid']
    
        mapoldnew_info['DowSubId'] = mapoldnew_info['ndownsubid']
    
    if UpdateStreamorder <0:
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
    
#######
def Copyfeature_to_another_shp_by_attribute(Source_shp,Target_shp,Col_NM='SubId',Values=[-1],Attributes = [-1]):
    layer_src=QgsVectorLayer(Source_shp,"")
    layer_trg=QgsVectorLayer(Target_shp,"")
    
    src_features = layer_src.getFeatures()  
    
    Selected_Features = []
    for sf in src_features:
        #centroidxy = sf.geometry().centroid().asPoint()
        Select_value = sf[Col_NM]
        if Select_value in Values:
            src_geometry =  sf.geometry()
            attribute = Attributes.loc[Attributes[Col_NM] == Select_value].values
            temp_feature=QgsFeature()
            temp_feature.setGeometry(src_geometry)
            temp_feature.setAttributes(attribute.tolist()[0])
            Selected_Features.append(temp_feature)
        
    layer_trg.startEditing()
    layer_trg.addFeatures(Selected_Features)
    layer_trg.commitChanges()
    layer_trg.updateExtents()
    del layer_src
    del layer_trg
###########

def ConnectLake_to_NonConnectLake_Updateinfo(NonC_Lakeinfo,finalriv_info,Merged_subids,Connect_Lake_ply_info,ConLakeId):
    
    All_Mergedcatchments = finalriv_info[finalriv_info['SubId'].isin(Merged_subids)]
    All_Mergedcatchments = All_Mergedcatchments.sort_values(["Strahler"], ascending = (True))
    Seg_IDS              = All_Mergedcatchments['Seg_ID'].values
    Seg_IDS              = np.unique(Seg_IDS)   
    routing_info         = All_Mergedcatchments[['SubId','DowSubId']].astype('float').values  

    for iseg in range(0,len(Seg_IDS)):
#            print('#########################################################################################33333')
        i_seg_id        = Seg_IDS[iseg]
        i_seg_info      = All_Mergedcatchments[All_Mergedcatchments['Seg_ID'] == i_seg_id]                      
        i_seg_info      = i_seg_info.sort_values(["Seg_order"], ascending = (True))    
        
        Lakeids_in_seg  = i_seg_info['HyLakeId'].values
        Lakeids_in_seg  = Lakeids_in_seg[Lakeids_in_seg > 0]
        Lakeids_in_seg  = np.unique(Lakeids_in_seg)
        
        
        if len(Lakeids_in_seg) <= 0:  ### No connected lake in this segment
            continue
            
        if ConLakeId > 0:
            Lakeids_in_seg = Lakeids_in_seg[Lakeids_in_seg != ConLakeId]
            
        Con_to_NonCon_Lakeids = Lakeids_in_seg            

        for i in range(0,len(Con_to_NonCon_Lakeids)):
            New_Non_Lakeid       = Con_to_NonCon_Lakeids[i]
            New_Lake_info        = Connect_Lake_ply_info.loc[Connect_Lake_ply_info['Hylak_id'] == New_Non_Lakeid]
            sub_coverd_bylake    = finalriv_info.loc[finalriv_info['HyLakeId'] == New_Non_Lakeid]
            sub_coverd_bylake    = sub_coverd_bylake.sort_values(["DA"], ascending = (True))
            tsubid               = sub_coverd_bylake['SubId'].values[len(sub_coverd_bylake) - 1]
            
            if sub_coverd_bylake['Seg_ID'].values[len(sub_coverd_bylake) - 1] != i_seg_id:
#                print(New_Non_Lakeid,sub_coverd_bylake['Seg_ID'].values[len(sub_coverd_bylake) - 1],i_seg_id)
                continue
                
            New_NonC_Lakeinfo = pd.DataFrame(np.full((1,len(NonC_Lakeinfo.columns)),np.nan), columns = NonC_Lakeinfo.columns)
            
            Upstreamcats         = Defcat(routing_info,tsubid)
            
            New_NonC_Lakeinfo.loc[0,'Gridcode']    = New_Non_Lakeid
            New_NonC_Lakeinfo.loc[0,'Gridcodes']   = New_Non_Lakeid
            New_NonC_Lakeinfo.loc[0,'value']       = New_Non_Lakeid
            New_NonC_Lakeinfo.loc[0,'DA_Area']     = sub_coverd_bylake['DA'].values[len(sub_coverd_bylake) - 1]
        
            New_NonC_Lakeinfo.loc[0,'SubId_riv']   = tsubid
            New_NonC_Lakeinfo.loc[0,'DownLakeID']  = -1.2345
        
            New_NonC_Lakeinfo.loc[0,'HyLakeId']    = New_Lake_info['Hylak_id'].values
            New_NonC_Lakeinfo.loc[0,'LakeVol']     = New_Lake_info['Vol_total'].values
            New_NonC_Lakeinfo.loc[0,'LakeDepth']   = New_Lake_info['Depth_avg'].values
            New_NonC_Lakeinfo.loc[0,'LakeArea']    = New_Lake_info['Lake_area'].values
            New_NonC_Lakeinfo.loc[0,'Laketype']    = New_Lake_info['Lake_type'].values
            
            mask1                                  = NonC_Lakeinfo['SubId_riv'].isin(Upstreamcats)
            mask2                                  = NonC_Lakeinfo['DownLakeID'].values < 0 
            
            Lake_To_Current_Lake                   = np.logical_and(mask1,mask2)
            
            
            NonC_Lakeinfo.loc[Lake_To_Current_Lake,'DownLakeID'] = New_Non_Lakeid
            
            DA_areas                                             = NonC_Lakeinfo.loc[Lake_To_Current_Lake,'DA_Area'].values
            New_NonC_Lakeinfo.loc[0,'Area_m']                    = sub_coverd_bylake['DA'].values[len(sub_coverd_bylake) - 1] - np.sum(DA_areas)
            
            NonC_Lakeinfo = pd.concat([NonC_Lakeinfo, New_NonC_Lakeinfo], ignore_index=True)
    
    
    return NonC_Lakeinfo    

############    
class LRRT:
    def __init__(self, dem_in = '#',dir_in = '#',hyshdply = '#',WidDep = '#',Lakefile = '#'
                                     ,Landuse = '#',Landuseinfo = '#',obspoint = '#',OutHyID = -1 ,OutHyID2 = -1, 
                                     OutputFolder = '#',ProjectNM = '#'):
        self.Path_dem_in = dem_in
        self.Path_dir_in = dir_in
        self.Path_hyshdply_in = hyshdply
        self.Path_WiDep_in = WidDep
        self.Path_Lakefile_in = Lakefile
        self.Path_Landuse_in = Landuse
        self.Path_Landuseinfo_in = Landuseinfo
        self.Path_obspoint_in = obspoint
        self.Path_OutputFolder = OutputFolder
        

        
        self.OutHyID = OutHyID
        self.OutHyID2 = OutHyID2
        
        self.ProjectNM = ProjectNM
        
        self.OutputFolder = os.path.join(self.Path_OutputFolder,self.ProjectNM)
        
        if not os.path.exists(self.OutputFolder):
	           os.makedirs(self.OutputFolder)         
        
        self.Raveinputsfolder = self.OutputFolder + '/'+'RavenInput/'
                
        self.qgisPP = os.environ['QGISPrefixPath']
        self.RoutingToolPath = os.environ['RoutingToolFolder']
                
        self.grassdb =os.path.join(tempfile.gettempdir(), 'grassdata_toolbox',self.ProjectNM)
        if not os.path.exists(self.grassdb):
	           os.makedirs(self.grassdb) 
        os.environ['GISDBASE'] = self.grassdb 
		
        self.grass_location_geo = 'Geographic'
        self.grass_location_pro = 'Projected'
        
        self.tempfolder = os.path.join(tempfile.gettempdir(), 'grassdata_toolbox_temp',self.ProjectNM)

        if not os.path.exists(self.tempfolder):
	           os.makedirs(self.tempfolder)
        
        
        self.sqlpath = os.path.join(self.grassdb,'Geographic\\PERMANENT\\sqlite\\sqlite.db')       
        self.cellSize = -9.9999
        self.SpRef_in = '#'
        self.ncols = -9999
        self.nrows = -9999
        self.Path_Maskply = '#'
        self.Path_dem = os.path.join(self.tempfolder,'dem.tif')
        self.Path_demproj = os.path.join(self.tempfolder,'dem_proj.tif')
        self.Path_allLakeply = os.path.join(self.tempfolder,'Hylake.shp')
        self.Path_WidDepLine = os.path.join(self.tempfolder,'WidDep.shp')
        self.Path_ObsPoint = os.path.join(self.tempfolder,'obspoint.shp')
        self.Path_Landuseinfo = os.path.join(self.tempfolder,'landuseinfo.csv')
        self.Path_allLakeRas = os.path.join(self.tempfolder,'hylakegdal.tif')
        self.Path_finalcatinfo_riv = os.path.join(self.tempfolder,'catinfo_riv.csv')
        self.Path_NonCLakeinfo = os.path.join(self.tempfolder,'NonC_Lakeinfo.csv')
        self.Path_NonCLakeinfo_cat = os.path.join(self.tempfolder,'NonC_Lakeinfo_cat.csv')
        self.Path_finalcatinfo_cat = os.path.join(self.tempfolder,'catinfo_cat.csv')
        self.Path_finalcatinfo = os.path.join(self.tempfolder,'catinfo.csv')
        self.Path_finalcatinfo_riv_type = os.path.join(self.tempfolder,'catinfo_riv.csvt')
        self.Path_finalcatinfo_type = os.path.join(self.tempfolder,'catinfo.csvt')
        self.Path_alllakeinfoinfo = os.path.join(self.tempfolder,'hylake.csv')
        self.Path_Maskply = os.path.join(self.tempfolder, 'HyMask2.shp')
       
########################################################################################
### Remove tempfolders

#########################################################################################
##### function to enerate mask based on the most donwstream polygon id
##### Output: For using hydroshed products the output is a mask polygon: self.Path_Maskply
#####         For using non hydroshed product, the output is a mask based on either dem or 
#####              a rough watershed delineaiton resut. depend on Path_OutletPoint value
#####              at the same time flow directio and accumulation is caculated for this region 
    def Generatmaskregion(self,Path_OutletPoint = '#'):
        #### g
        ### Set up QGIS enviroment 
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

        if self.OutHyID > 0:
            hyinfocsv = self.Path_hyshdply_in[:-3] + "dbf"
            tempinfo = Dbf5(hyinfocsv)
            hyshdinfo = tempinfo.to_dataframe().values 
           
            HydroBasins1 = Defcat(hyshdinfo,self.OutHyID) ### return fid of polygons that needs to be select 
            if self.OutHyID > 0:
                HydroBasins2 = Defcat(hyshdinfo,self.OutHyID2)            
    ###  exculde the Ids in HydroBasins2 from HydroBasins1
                for i in range(len(HydroBasins2)):
                    rows =np.argwhere(HydroBasins1 == HydroBasins2[i])
                    HydroBasins1 = np.delete(HydroBasins1, rows)
                HydroBasins = HydroBasins1            
            else:
                HydroBasins = HydroBasins1
        
    ### Load HydroSHED Layers 
            hyshedl12 = QgsVectorLayer(self.Path_hyshdply_in, "")
        
    ### Build qgis selection expression
            where_clause = '"HYBAS_ID" IN'+ " ("
            for i in range(0,len(HydroBasins)):
                if i == 0:
                    where_clause = where_clause + str(HydroBasins[i])
                else:
                    where_clause = where_clause + "," + str(HydroBasins[i])
            where_clause = where_clause + ")"
    
            req = QgsFeatureRequest().setFlags( QgsFeatureRequest.NoGeometry)
            req.setFilterExpression(where_clause)
            it = hyshedl12.getFeatures( req )
        
        ### obtain all feature id of selected polygons
            selectedFeatureID = []
    
            for feature in it:
    #        print(feature.id())
                selectedFeatureID.append(feature.id())
            
            hyshedl12.select(selectedFeatureID)   ### select with polygon id
        
        # Save selected polygons to output 
            _writer = QgsVectorFileWriter.writeAsVectorFormat(hyshedl12, os.path.join(self.tempfolder, 'HyMask.shp'), "UTF-8", hyshedl12.crs(), "ESRI Shapefile", onlySelected=True)
            processing.run('gdal:dissolve', {'INPUT':os.path.join(self.tempfolder, 'HyMask.shp'),'FIELD':'MAIN_BAS','OUTPUT':self.Path_Maskply})
            
            del hyshedl12
            
        else:
            r_dem_layer = QgsRasterLayer(self.Path_dem_in, "") ### load DEM raster as a  QGIS raster object to obtain attribute        
            self.cellSize = float(r_dem_layer.rasterUnitsPerPixelX())  ### Get Raster cell size
            self.SpRef_in = r_dem_layer.crs().authid()   ### get Raster spatialReference id
            params = {'INPUT': self.Path_dem_in, 'format': 'GTiff', 'OUTPUT': self.Path_dem}
            processing.run('gdal:translate',params)
            import grass.script as grass
            from grass.script import array as garray
            import grass.script.setup as gsetup
            from grass.pygrass.modules.shortcuts import general as g
            from grass.pygrass.modules.shortcuts import raster as r
            from grass.pygrass.modules import Module
            from grass_session import Session
            
            os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='1'))
            PERMANENT = Session()
            PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
            
            grass.run_command("r.import", input = self.Path_dem, output = 'dem', overwrite = True)
            grass.run_command('g.region', raster='dem')
            #### create watershed mask or use dem as mask
            if Path_OutletPoint != '#':
                grass.run_command('r.watershed',elevation = 'dem', drainage = 'dir_grass',accumulation = 'acc_grass2',flags = 's', overwrite = True)
                grass.run_command('r.mapcalc',expression = "acc_grass = abs(acc_grass2@PERMANENT)",overwrite = True)
                grass.run_command('r.water.outlet',input = 'dir_grass', output = 'wat_mask', coordinates  = Path_OutletPoint,overwrite = True)
                grass.run_command('r.mask'  , raster='wat_mask', maskcats = '*',overwrite = True)
                grass.run_command('g.region', raster='wat_mask',overwrite = True)
            else:
                grass.run_command('r.mask'  , raster='dem', maskcats = '*',overwrite = True)
            
            grass.run_command('r.out.gdal', input = 'MASK',output = os.path.join(self.tempfolder, 'Mask1.tif'),format= 'GTiff',overwrite = True)
            processing.run("gdal:polygonize", {'INPUT':os.path.join(self.tempfolder, 'Mask1.tif'),'BAND':1,'FIELD':'DN','EIGHT_CONNECTEDNESS':False,'EXTRA':'','OUTPUT':os.path.join(self.tempfolder, 'HyMask.shp')})
            processing.run('gdal:dissolve', {'INPUT':os.path.join(self.tempfolder, 'HyMask.shp'),'FIELD':'DN','OUTPUT':self.Path_Maskply})
            PERMANENT.close()
            del r_dem_layer
        Qgs.exit()
        
        return 

##################################################################################################  
#### functions to preprocess data, Output:
##  self.Path_dem,self.cellSize,self.SpRef_in
##  
    def Generateinputdata(self):

        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing
        from processing.tools import dataobjects
           
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        context = dataobjects.createContext()
        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)
    
        r_dem_layer = QgsRasterLayer(self.Path_dem_in, "") ### load DEM raster as a  QGIS raster object to obtain attribute        
        self.cellSize = float(r_dem_layer.rasterUnitsPerPixelX())  ### Get Raster cell size
        self.SpRef_in = r_dem_layer.crs().authid()   ### get Raster spatialReference id
        print("Working with a  sptail reference  :   " , r_dem_layer.crs().description(), "      ", self.SpRef_in )
        print("The cell cize is   ",self.cellSize)
    
        if  self.OutHyID > 0:   #### input is using hydroshed  hydroshed polygons 
            params = {'INPUT': self.Path_dem_in,'MASK': self.Path_Maskply,'NODATA': -9999,'ALPHA_BAND': False,'CROP_TO_CUTLINE': True,
                                                                    'KEEP_RESOLUTION': True,
                                                                    'OPTIONS': 'COMPRESS=LZW',
                                                                    'DATA_TYPE': 0,  # Byte
                                                                    'OUTPUT': self.Path_dem}
            dem = processing.run('gdal:cliprasterbymasklayer',params)  #### extract dem
        
        copyfile(os.path.join(self.RoutingToolPath,'catinfo_riv.csvt'),os.path.join(self.tempfolder,'catinfo_riv.csvt')) 
        copyfile(os.path.join(self.RoutingToolPath,'catinfo_riv.csvt'),os.path.join(self.tempfolder,'catinfo_cat.csvt')) 
        copyfile(os.path.join(self.RoutingToolPath,'Nonlakeinfo.csvt'),os.path.join(self.tempfolder,'Nonlakeinfo.csvt'))
        copyfile(os.path.join(self.RoutingToolPath,'Nonlakeinfo.csvt'),os.path.join(self.tempfolder,'NonC_Lakeinfo_cat.csvt'))  
###### set up GRASS environment for translate vector to rasters and clip rasters
        import grass.script as grass
        from grass.script import array as garray
        from grass.script import core as gcore
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session
    
       
        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        ###  hydroshed product use cliped dem to setup g.region and g.mask 
        ###  non hydroshed product the g.region and g.mask is defined in generatemask funciton 
        if self.OutHyID > 0:   
            grass.run_command("r.import", input = self.Path_dem, output = 'dem', overwrite = True)
            grass.run_command('r.mask'  , raster='dem', maskcats = '*',overwrite = True)
            grass.run_command('g.region', raster='dem')

        strtemp_array = garray.array(mapname="dem")
        self.ncols = int(strtemp_array.shape[1])
        self.nrows = int(strtemp_array.shape[0])
        grsregion = gcore.region()
        ### process vector data, clip and import
        processing.run("native:clip", {'INPUT':self.Path_Lakefile_in,'OVERLAY':self.Path_Maskply,'OUTPUT':self.Path_allLakeply},context = context)
        processing.run("native:clip", {'INPUT':self.Path_WiDep_in,'OVERLAY':self.Path_Maskply,'OUTPUT':self.Path_WidDepLine},context = context)
        processing.run("native:clip", {'INPUT':self.Path_obspoint_in,'OVERLAY':self.Path_Maskply,'OUTPUT':self.Path_ObsPoint},context = context)
        grass.run_command("v.import", input = self.Path_WidDepLine, output = 'WidDep', overwrite = True)
        grass.run_command("v.import", input = self.Path_ObsPoint, output = 'obspoint', overwrite = True)
        grass.run_command("v.import", input = self.Path_allLakeply, output = 'Hylake', overwrite = True)


        if self.Path_dir_in != '#':  ### Hydroshed provided with dir 
            grass.run_command("r.external", input = self.Path_dir_in, output = 'dir_in', overwrite = True)
            grass.run_command("r.clip", input = 'dir_in', output = 'dir_Arcgis', overwrite = True, flags = 'r')
            grass.run_command('r.reclass', input='dir_Arcgis',output = 'dir_grass',rules =os.path.join(self.RoutingToolPath,'Arcgis2GrassDIR.txt'), overwrite = True)
        else:  ### non hydroshed if dir has been build 
            dirbuid = grass.find_file('dir_grass', element = 'cell')
            if (len(dirbuid['file']) <=0):
                grass.run_command('r.watershed',elevation = 'dem', drainage = 'dir_grass', accumulation = 'acc_grass2',flags = 's', overwrite = True)
                grass.run_command('r.mapcalc',expression = "acc_grass = abs(acc_grass2@PERMANENT)",overwrite = True)
            grass.run_command('r.reclass', input='dir_grass',output = 'dir_Arcgis',rules = os.path.join(self.RoutingToolPath,'Grass2ArcgisDIR.txt'), overwrite = True)
    
        copyfile(self.Path_Landuseinfo_in, self.Path_Landuseinfo) 
        grass.run_command("r.external", input = self.Path_Landuse_in, output = 'landuse_in', overwrite = True)
        grass.run_command("r.clip", input = 'landuse_in', output = 'landuse', overwrite = True)        
    
        ### change lake vector to lake raster.... with -at option 
        os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(self.cellSize) + "  " +str(self.cellSize)+'  -te   '+ str(grsregion['w'])+"   " +str(grsregion['s'])+"   " +str(grsregion['e'])+"   " +str(grsregion['n'])+"   " + "\"" +  self.Path_allLakeply +"\""+ "    "+ "\""+self.Path_allLakeRas+"\"")
        grass.run_command("r.in.gdal", input = self.Path_allLakeRas, output = 'alllakeraster_in', overwrite = True)
        grass.run_command('r.mapcalc',expression = 'alllake = int(alllakeraster_in)',overwrite = True)
        grass.run_command("r.null", map = 'alllake', setnull = -9999)
        ### rasterize other vectors 
        grass.run_command('v.to.rast',input = 'WidDep',output = 'width',use = 'attr',attribute_column = 'WIDTH',overwrite = True)
        grass.run_command('v.to.rast',input = 'WidDep',output = 'depth',use = 'attr',attribute_column = 'DEPTH',overwrite = True)
        grass.run_command('v.to.rast',input = 'WidDep',output = 'qmean',use = 'attr',attribute_column = 'Q_Mean',overwrite = True)
        grass.run_command('v.to.rast',input = 'obspoint',output = 'obs',use = 'attr',attribute_column = 'Obs_ID',overwrite = True)

        PERMANENT.close()
        del r_dem_layer
        Qgs.exit()
#####################################################################################################

####################################################################################################3
                
    def WatershedDiscretizationToolset(self,accthresold):
        import grass.script as grass
        from grass.script import array as garray
        from grass.script import core as gcore
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)

        grsregion = gcore.region()        
        if self.Path_dir_in == '#':  ### did not provide dir, use dem to generate watershed. recommand !!
            grass.run_command('r.stream.extract',elevation = 'dem',accumulation = 'acc_grass',threshold =accthresold,stream_raster = 'str_grass_r',
                              stream_vector = 'str_grass_v',overwrite = True)
        else:
        ## generate correct stream raster, when the dir is not derived from dem. for Hydroshed Cases 
            grass.run_command('r.accumulate', direction='dir_grass',format = '45degree',accumulation ='acc_grass',
                            stream = 'str_grass_v',threshold = accthresold, overwrite = True)
            grass.run_command('v.to.rast',input = 'str_grass_v',output = 'str_grass_r2',use = 'cat',overwrite = True)
    
            strtemp_array = garray.array(mapname="str_grass_r2")
            acc_array = np.absolute(garray.array(mapname="acc_grass"))
            dirarc_array = garray.array(mapname="dir_Arcgis")
    
    #####Correct stream network
            strids = np.unique(strtemp_array)
            strids = strids[strids >= 0] 

        
            for i in range(0,len(strids)):
                strid = strids[i]

                trow,tcol = Getbasinoutlet(strid,strtemp_array,acc_array,dirarc_array,self.nrows,self.ncols)
                nrow,ncol = Nextcell(dirarc_array,trow,tcol)### get the downstream catchment id
                if nrow <= 0 or ncol <=0 or nrow >=self.nrows or ncol >= self.ncols:
                    continue
                nstrid = strtemp_array[nrow,ncol]
        
                rowcol = np.argwhere(strtemp_array==nstrid).astype(int)
                catacc = np.full((len(rowcol),3),-9999)
                catacc[:,0] = rowcol[:,0]
                catacc[:,1] = rowcol[:,1]
                catacc[:,2] = acc_array[rowcol[:,0],rowcol[:,1]]
                catacc = catacc[catacc[:,2].argsort()]
    
                if nrow != catacc[0,0] or ncol != catacc[0,1]:  ### If the stream 1 did not connect to the begining of the downstream stream (2) 
                    nnrow,nncol = Nextcell(dirarc_array,nrow,ncol) ### The last cell of 2 will be changed. 
                    if nnrow <= 0 or nncol <=0 or nnrow >=self.nrows or nncol >= self.ncols:
                        continue
                    nnstrid = strtemp_array[nnrow,nncol]
                    strtemp_array[nrow,ncol] = nnstrid
        ##### end modify stream grid and store new grids
            temparray = garray.array()
            temparray[:,:] = 0
            temparray[:,:] = strtemp_array[:,:]
            temparray.write(mapname="str_grass_rf", overwrite=True)
            grass.run_command('r.null', map='str_grass_rf',setnull=0)
            grass.run_command('r.mapcalc',expression = 'str_grass_r = int(str_grass_rf)',overwrite = True)
             
#        grass.run_command('r.thin',input = 'str_grass_rfn', output = 'str_grass_r',overwrite = True)
        grass.run_command('r.to.vect',  input = 'str_grass_r',output = 'str', type ='line' ,overwrite = True)

        ##### generate catchment without lakes based on 'str_grass_r'        
        grass.run_command('r.stream.basins',direction = 'dir_grass', stream = 'str_grass_r', basins = 'cat1',overwrite = True)

        
################ check connected lakes  and non connected lakes 
#        grass.run_command('v.out.ogr', input = 'str_grass_v',output = os.path.join(self.tempfolder, "str_grass_v.shp"),format= 'ESRI_Shapefile',overwrite = True)
        grass.run_command('v.select',ainput = 'Hylake',binput = 'str',output = 'lake_str',overwrite = True)
        grass.run_command('v.out.ogr', input = 'lake_str',output = os.path.join(self.tempfolder, "Connect_lake.shp"),format= 'ESRI_Shapefile',overwrite = True)
        os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(self.cellSize) + "  " +str(self.cellSize) +'   -te   '+ str(grsregion['w'])+"   " +str(grsregion['s'])+"   " +str(grsregion['e'])+"   " +str(grsregion['n'])+"   " +"\"" +  os.path.join(self.tempfolder, "Connect_lake.shp") +"\""+ "    "+ "\""+os.path.join(self.tempfolder, "cnhylakegdal.tif")+"\"")
        grass.run_command("r.in.gdal", input = os.path.join(self.tempfolder, "cnhylakegdal.tif"), output = 'cnlakeraster_in', overwrite = True)
        grass.run_command('r.mapcalc',expression = 'Connect_Lake = int(cnlakeraster_in)',overwrite = True)    
        grass.run_command('r.mapcalc',expression = 'Nonconnect_Lake = if(isnull(Connect_Lake),alllake,null())',overwrite = True)
    
#        grass.run_command('v.out.ogr', input = 'finalcat_F',output = os.path.join(self.tempfolder,'finalcat_info1.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
        grass.run_command('r.out.gdal', input = 'cat1',output = os.path.join(self.tempfolder,'cat1.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture')    
        grass.run_command('r.out.gdal', input = 'Connect_Lake',output = os.path.join(self.tempfolder,'Connect_Lake.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
        grass.run_command('r.out.gdal', input = 'Nonconnect_Lake',output = os.path.join(self.tempfolder,'Nonconnect_Lake.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture')    
        grass.run_command('r.out.gdal', input = 'str_grass_r',output = os.path.join(self.tempfolder,'str_grass_r.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture')    

        PERMANENT.close()
###########################################################################################3

############################################################################################
    def AutomatedWatershedsandLakesFilterToolset(self,Thre_Lake_Area_Connect = 0,Thre_Lake_Area_nonConnect = -1,MaximumLakegrids = 1600):

        tempinfo = Dbf5(self.Path_allLakeply[:-3] + "dbf")
        allLakinfo = tempinfo.to_dataframe()
        allLakinfo.to_csv(self.Path_alllakeinfoinfo,index = None, header=True)
        
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        con = sqlite3.connect(self.sqlpath)
        
        VolThreshold =Thre_Lake_Area_Connect ### lake area thresthold for connected lakes 
        NonConLThres = Thre_Lake_Area_nonConnect ### lake area thresthold for non connected lakes 
        nlakegrids = MaximumLakegrids 
    
        temparray = garray.array()
        temparray[:,:] = -9999
        temparray.write(mapname="tempraster", overwrite=True)
        self.ncols = int(temparray.shape[1])
        self.nrows = int(temparray.shape[0])
##### begin processing 

        blid  = 1000000    #### begining of new lake id
        bcid  = 1          ## begining of new cat id of hydrosheds
        bsid  = 2000000    ## begining of new cat id of in inflow of lakes
        blid2 = 3000000
        boid  = 4000000
    
        noncnlake_arr = garray.array(mapname="Nonconnect_Lake")
        conlake_arr   = garray.array(mapname="Connect_Lake")
        cat1_arr      = garray.array(mapname="cat1")
        str_array     = garray.array(mapname="str_grass_r")
        acc_array     = np.absolute(garray.array(mapname="acc_grass"))
        dir_array     = garray.array(mapname="dir_Arcgis")
        obs_array     = garray.array(mapname="obs")

###### generate selected lakes 
        hylake1 = selectlake2(conlake_arr,VolThreshold,allLakinfo) ### remove lakes with lake area smaller than the VolThreshold from connected lake raster 
        if NonConLThres >= 0:
            Lake1 = selectlake(hylake1,noncnlake_arr,NonConLThres,allLakinfo) ### remove lakes with lake area smaller than the NonConLThres from non-connected lake raster 
        else:
            Lake1 = hylake1
        temparray[:,:] = Lake1[:,:]
        temparray.write(mapname="SelectedLakes", overwrite=True)
        grass.run_command('r.null', map='SelectedLakes',setnull=-9999)

####    
        Pourpoints = GenerPourpoint(cat1_arr,Lake1,str_array,self.nrows,self.ncols,blid,bsid,bcid,acc_array,dir_array)
        temparray[:,:] = Pourpoints[:,:]
        temparray.write(mapname="Pourpoints_1", overwrite=True)
        grass.run_command('r.null', map='Pourpoints_1',setnull=-9999)
        grass.run_command('r.to.vect', input='Pourpoints_1',output='Pourpoints_1_F',type='point', overwrite = True)
### cat2        
        grass.run_command('r.stream.basins',direction = 'dir_grass', points = 'Pourpoints_1_F', basins = 'cat2_t',overwrite = True)        
        sqlstat="SELECT cat, value FROM Pourpoints_1_F"
        df_P_1_F = pd.read_sql_query(sqlstat, con)
        df_P_1_F.loc[len(df_P_1_F),'cat'] = '*'
        df_P_1_F.loc[len(df_P_1_F)-1,'value'] = 'NULL'
        df_P_1_F['eq'] = '='
        df_P_1_F.to_csv(os.path.join(self.tempfolder,'rule_cat2.txt'),sep = ' ',columns = ['cat','eq','value'],header = False,index=False)
        grass.run_command('r.reclass', input='cat2_t',output = 'cat2',rules =os.path.join(self.tempfolder,'rule_cat2.txt'), overwrite = True)    
        grass.run_command('v.out.ogr', input = 'Pourpoints_1_F',output = os.path.join(self.tempfolder, "Pourpoints_1_F.shp"),format= 'ESRI_Shapefile',overwrite = True)
#        grass.run_command('r.out.gdal', input = 'cat1',output = os.path.join(self.tempfolder,'cat1.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture')    
# cat3   
        cat2_array =  garray.array(mapname="cat2")
        temcat,outlakeids =CE_mcat4lake(cat2_array,Lake1,acc_array,dir_array,bsid,self.nrows,self.ncols,Pourpoints) 
        temcat2 = CE_Lakeerror(acc_array,dir_array,Lake1,temcat,bsid,blid,boid,self.nrows,self.ncols,cat1_arr)
        temparray[:,:] = temcat2[:,:]
        temparray.write(mapname="cat3", overwrite=True)
        grass.run_command('r.null', map='cat3',setnull=-9999)

        grass.run_command('r.out.gdal', input = 'cat3',output = os.path.join(self.tempfolder,'cat3.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
#        grass.run_command('r.out.gdal', input = 'cat3',output = os.path.join(self.tempfolder,'cat3.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
    
        nPourpoints = GenerateFinalPourpoints(acc_array,dir_array,Lake1,temcat2,bsid,blid,boid,self.nrows,self.ncols,cat1_arr,obs_array)
        temparray[:,:] = nPourpoints[:,:]
        temparray.write(mapname="Pourpoints_2", overwrite=True)
        grass.run_command('r.null', map='Pourpoints_2',setnull=-9999)
        grass.run_command('r.to.vect', input='Pourpoints_2',output='Pourpoints_2_F',type='point', overwrite = True)    
        grass.run_command('v.out.ogr', input = 'Pourpoints_2_F',output = os.path.join(self.tempfolder, "Pourpoints_2_F.shp"),format= 'ESRI_Shapefile',overwrite = True)

        ## with new pour points
# cat 4
        grass.run_command('r.stream.basins',direction = 'dir_grass', points = 'Pourpoints_2_F', basins = 'cat4_t',overwrite = True)        
        sqlstat="SELECT cat, value FROM Pourpoints_2_F"
        df_P_1_F = pd.read_sql_query(sqlstat, con)
        df_P_1_F.loc[len(df_P_1_F),'cat'] = '*'
        df_P_1_F.loc[len(df_P_1_F)-1,'value'] = 'NULL'
        df_P_1_F['eq'] = '='
        df_P_1_F.to_csv(os.path.join(self.tempfolder,'rule_cat4.txt'),sep = ' ',columns = ['cat','eq','value'],header = False,index=False)
        grass.run_command('r.reclass', input='cat4_t',output = 'cat4',rules =os.path.join(self.tempfolder,'rule_cat4.txt'), overwrite = True)
        cat4_array =  garray.array(mapname="cat4")
        grass.run_command('r.out.gdal', input = 'cat4',output = os.path.join(self.tempfolder,'cat4.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
                
        outlakeids= check_lakecatchment(cat4_array,Lake1,acc_array,dir_array,bsid,self.nrows,self.ncols)


    
        ### Modify lake flow directions
        ndir = ChangeDIR(dir_array,Lake1,acc_array,self.ncols,self.nrows,outlakeids,nlakegrids,cat4_array)
        temparray[:,:] = ndir[:,:]
        temparray.write(mapname="ndir_Arcgis", overwrite=True)
        grass.run_command('r.null', map='ndir_Arcgis',setnull=-9999)    
        grass.run_command('r.reclass', input='ndir_Arcgis',output = 'ndir_grass',rules =os.path.join(self.RoutingToolPath,'Arcgis2GrassDIR.txt'),overwrite = True)
# cat5     
        grass.run_command('r.stream.basins',direction = 'ndir_grass', points = 'Pourpoints_2_F', basins = 'cat5_t',overwrite = True)
        sqlstat="SELECT cat, value FROM Pourpoints_2_F"
        df_P_2_F = pd.read_sql_query(sqlstat, con)
        df_P_2_F.loc[len(df_P_2_F),'cat'] = '*'
        df_P_2_F.loc[len(df_P_2_F)-1,'value'] = 'NULL'
        df_P_2_F['eq'] = '='
        df_P_2_F.to_csv(os.path.join(self.tempfolder,'rule_cat5.txt'),sep = ' ',columns = ['cat','eq','value'],header = False,index=False)
        grass.run_command('r.reclass', input='cat5_t',output = 'cat5',rules =os.path.join(self.tempfolder,'rule_cat5.txt'), overwrite = True)
        grass.run_command('r.out.gdal', input = 'cat5',output = os.path.join(self.tempfolder,'cat5.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
        cat5_array =  garray.array(mapname="cat5")
        rowcols = np.argwhere(cat5_array == 0)
        cat5_array[rowcols[:,0],rowcols[:,1]] = -9999
        
        
        finalcat,Non_con_lake_cat = CE_mcat4lake2(cat5_array,Lake1,acc_array,dir_array,bsid,self.nrows,self.ncols,nPourpoints,noncnlake_arr)
        temparray[:,:] = finalcat[:,:]
        temparray.write(mapname="finalcat", overwrite=True)
        grass.run_command('r.null', map='finalcat',setnull=-9999)    
        grass.run_command('r.to.vect', input='finalcat',output='finalcat_F1',type='area', overwrite = True)    
        temparray[:,:] = Non_con_lake_cat[:,:]
        temparray.write(mapname="Non_con_lake_cat", overwrite=True)
        grass.run_command('r.null', map='Non_con_lake_cat',setnull=-9999)   
         
####   dissolve final catchment polygons    
        grass.run_command('v.db.addcolumn', map= 'finalcat_F1', columns = "Gridcode VARCHAR(40)")
        grass.run_command('v.db.update', map= 'finalcat_F1', column = "Gridcode",qcol = 'value')
#    grass.run_command('v.reclass', input= 'finalcat_F1', column = "Gridcode",output = 'finalcat_F',overwrite = True)
        grass.run_command('v.dissolve', input= 'finalcat_F1', column = "Gridcode",output = 'finalcat_F',overwrite = True)   
#    grass.run_command('v.select',ainput = 'str',binput = 'finalcat_F',output = 'str_finalcat',overwrite = True)
        grass.run_command('v.overlay',ainput = 'str',binput = 'finalcat_F1',operator = 'and',output = 'str_finalcat',overwrite = True)  
        grass.run_command('v.out.ogr', input = 'str_finalcat',output = os.path.join(self.tempfolder,'str_finalcat.shp'),format= 'ESRI_Shapefile',overwrite = True)
        grass.run_command('v.out.ogr', input = 'finalcat_F1',output = os.path.join(self.tempfolder,'finalcat_F1.shp'),format= 'ESRI_Shapefile',overwrite = True)                
        grass.run_command('v.out.ogr', input = 'str',output = os.path.join(self.tempfolder,'str.shp'),format= 'ESRI_Shapefile',overwrite = True)    
        grass.run_command('r.to.vect', input='SelectedLakes',output='SelectedLakes_F',type='area', overwrite = True) 
        grass.run_command('v.out.ogr', input = 'str_finalcat',output = os.path.join(self.tempfolder,'str_finalcat.shp'),format= 'ESRI_Shapefile',overwrite = True)
#    grass.run_command('v.db.join', map= 'SelectedLakes_F',column = 'value', other_table = 'result',other_column ='SubId', overwrite = True)
        con.close()
        PERMANENT.close()
                
############################################################################3
    def RoutingNetworkTopologyUpdateToolset_riv(self,projection = 'default'):
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing
        from processing.tools import dataobjects
           
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        context = dataobjects.createContext()
        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)


        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        
        con = sqlite3.connect(self.sqlpath)
      
        #### generate new stream raster based on str_grass_r from acc thresthold and finalcat.
         
        finalcat_arr      = garray.array(mapname="finalcat")          # final catchment raster generated by lake filter tool set
        str_grass_r_array = garray.array(mapname="str_grass_r")  # river networkt raster, each river has unique id
        acc_array         = np.absolute(garray.array(mapname="acc_grass"))
        NonCL_array       = garray.array(mapname="Non_con_lake_cat")
        conlake_arr       = garray.array(mapname="Connect_Lake")
        self.ncols        = int(finalcat_arr.shape[1])                  # obtain rows and cols
        self.nrows        = int(finalcat_arr.shape[0])  
        
        nstrarray = garray.array()
        nstrarray[:,:] = -9999        
        ## loop for each str_grass_r
        
        ## obtain all stream id 
        strids = np.unique(str_grass_r_array)#
        strids = strids[strids>0]
        #obtain all cat id 
        fcatids = np.unique(finalcat_arr)#### cat all catchment idd
        fcatids = fcatids[fcatids>0]
        ###A array store new stream id , old stream id and catchemtn id 
        nstrinfo = np.full((len(strids)*len(fcatids),5),-9999)   ### maximum possible number of new stream segments 
        nstrinfodf = pd.DataFrame(nstrinfo, columns = ['No_Lake_SubId', "Finalcat_SubId","With_Lake_SubId","MinAcc","MaxAcc"])
        
        ### assign a unique stream id based on each combination of str_grass_r and final cat
        nstrid = 1
        idx    = 0
        for i in range(0,len(strids)):
            strid = strids[i]
            strmask = str_grass_r_array == strid  ### mask array, true at location of str == strid 
            catsinstr = finalcat_arr[strmask]  ### get all cat id overlay with this stream
            catsinstrids = np.unique(catsinstr)  ### get unique cat id overlay with this stream
            catsinstrids = catsinstrids[catsinstrids>0]
            for j in range(0,len(catsinstrids)):   ### loop for each cat id, and assgin for a new stream id 
                finalcatmask = finalcat_arr == catsinstrids[j]
                mask = np.logical_and(finalcatmask, strmask)   ### build mask for grids belong to current stream and current cat
                nstrarray[mask] = nstrid   ## assgin new stream id to mask region
                stracc          = acc_array[mask]
                nstrinfodf.loc[idx,'No_Lake_SubId']   = strid
                nstrinfodf.loc[idx,'Finalcat_SubId']  = catsinstrids[j]
                nstrinfodf.loc[idx,'With_Lake_SubId'] = nstrid
                nstrinfodf.loc[idx,'MinAcc']          = np.amin(stracc)
                nstrinfodf.loc[idx,'MaxAcc']          = np.amax(stracc)
                nstrid = nstrid + 1
                idx    = idx    + 1
            
            i_nstrinfodf = nstrinfodf.loc[nstrinfodf['No_Lake_SubId'] == strid]
            i_nstrinfodf = i_nstrinfodf.sort_values(["MinAcc"], ascending = (True))
#            print("#####################################################################3")
#            print(len(catsinstrids))
#            print(i_nstrinfodf)
            for j in range(0,len(i_nstrinfodf) - 1):
                if i_nstrinfodf['MaxAcc'].values[j] > i_nstrinfodf['MaxAcc'].values[j +1]:
                    mask = nstrarray == i_nstrinfodf['With_Lake_SubId'].values[j+1]   ### build mask for grids belong to current stream and current cat
                    nstrarray[mask]   = i_nstrinfodf['With_Lake_SubId'].values[j]   ## assgin new stream id to mask region                    
                    connected_Lakeids = conlake_arr[mask]
                    connected_Lakeids = connected_Lakeids[connected_Lakeids > 0]
                    connected_Lakeids = np.unique(connected_Lakeids)
                    if(len(connected_Lakeids) == 1):
                        catid                       = i_nstrinfodf['Finalcat_SubId'].values[j+1]
                        finalcatmask2               = finalcat_arr == catid
                        NonCL_array[finalcatmask2]  = connected_Lakeids
        
            
        temparray = garray.array()
        temparray[:,:] = -9999
        temparray[:,:] = NonCL_array[:,:]
        temparray.write(mapname="Non_con_lake_cat_t", overwrite=True)
        grass.run_command('r.null', map='Non_con_lake_cat_t',setnull=-9999)
        grass.run_command('r.mapcalc',expression = 'Non_con_lake_cat = int(Non_con_lake_cat_t)',overwrite = True)
        grass.run_command('r.to.vect', input='Non_con_lake_cat',output='Non_con_lake_cat_t',type='area', overwrite = True)
        grass.run_command('v.db.addcolumn', map= 'Non_con_lake_cat_t', columns = "Gridcodes VARCHAR(40)")
        grass.run_command('v.db.update', map= 'Non_con_lake_cat_t', column = "Gridcodes",qcol = 'value')
        grass.run_command('v.dissolve', input= 'Non_con_lake_cat_t', column = "Gridcodes",output = 'Non_con_lake_cat_1',overwrite = True)    
        grass.run_command('v.out.ogr', input = 'Non_con_lake_cat_1',output = os.path.join(self.tempfolder,'Non_con_lake_cat_1.shp'),format= 'ESRI_Shapefile',overwrite = True)

                
        
        temparray = garray.array()
        temparray[:,:] = -9999
        temparray[:,:] = nstrarray[:,:]
        temparray.write(mapname="nstr_seg_t", overwrite=True)  #### write new stream id to a grass raster 
        grass.run_command('r.null', map='nstr_seg_t',setnull=-9999)
        grass.run_command('r.mapcalc',expression = 'nstr_seg = int(nstr_seg_t)',overwrite = True)
        grass.run_command('r.out.gdal', input = 'nstr_seg',output =os.path.join(self.tempfolder,'nstr_seg.tif'),format= 'GTiff',overwrite = True)  
        
        ### Generate new catchment based on new stream   
        grass.run_command('r.stream.basins',direction = 'ndir_grass', stream = 'nstr_seg', basins = 'Net_cat',overwrite = True)        
        grass.run_command('r.out.gdal', input = 'Net_cat',output =os.path.join(self.tempfolder,'Net_cat.tif'),format= 'GTiff',overwrite = True)
        grass.run_command('r.to.vect', input='Net_cat',output='Net_cat_F1',type='area', overwrite = True)    ## save to vector 
        grass.run_command('v.db.addcolumn', map= 'Net_cat_F1', columns = "GC_str VARCHAR(40)")
        grass.run_command('v.db.addcolumn', map= 'Net_cat_F1', columns = "Area_m double")        
        grass.run_command('v.db.update', map= 'Net_cat_F1', column = "GC_str",qcol = 'value')
        grass.run_command('v.dissolve', input= 'Net_cat_F1', column = "GC_str",output = 'Net_cat_F',overwrite = True) ### dissolve based on gridcode 
        grass.run_command('v.out.ogr', input = 'Net_cat_F',output = os.path.join(self.tempfolder,'Net_cat_F.shp'),format= 'ESRI_Shapefile',overwrite = True)
        grass.run_command('v.db.addcolumn', map= 'Net_cat_F', columns = "Gridcode INT")        
        grass.run_command('v.db.update', map= 'Net_cat_F', column = "Gridcode",qcol = 'GC_str')
        

        # create a river network shpfile that have river segment based on   'Net_cat_F.shp'
#        grass.run_command('r.mapcalc',expression = 'str1 = if(isnull(str_grass_r),null(),1)',overwrite = True) ## a rive 
#        grass.run_command('r.to.vect',  input = 'str1',output = 'str1', type ='line' ,overwrite = True)

        ## obtain a stream vector, segmentation based on new catchment polygon 
        if self.Path_dir_in == '#':
            grass.run_command('v.overlay',ainput = 'str_grass_v',alayer = 2, atype = 'line', binput = 'Net_cat_F',operator = 'and',output = 'nstr_nfinalcat',overwrite = True) 
        else:
            grass.run_command('v.overlay',ainput = 'str_grass_v', binput = 'Net_cat_F',operator = 'and',output = 'nstr_nfinalcat',overwrite = True)  
        grass.run_command('v.out.ogr', input = 'nstr_nfinalcat',output = os.path.join(self.tempfolder,'nstr_nfinalcat.shp'),format= 'ESRI_Shapefile',overwrite = True)
        processing.run('gdal:dissolve', {'INPUT':os.path.join(self.tempfolder, 'nstr_nfinalcat.shp'),'FIELD':'b_GC_str','OUTPUT':os.path.join(self.tempfolder, "nstr_nfinalcat_F.shp")}) 
        grass.run_command("v.import", input =os.path.join(self.tempfolder, "nstr_nfinalcat_F.shp"), output = 'nstr_nfinalcat_F', overwrite = True)               
  

        
        grass.run_command('v.db.addcolumn', map= 'nstr_nfinalcat_F', columns = "Gridcode INT")
        grass.run_command('v.db.addcolumn', map= 'nstr_nfinalcat_F', columns = "Length_m double")
        grass.run_command('v.db.addcolumn', map = 'Non_con_lake_cat_1',columns = "Area_m double")
        grass.run_command('v.db.addcolumn', map= 'Non_con_lake_cat_1', columns = "Gridcode INT")
        grass.run_command('v.db.update', map= 'Non_con_lake_cat_1', column = "Gridcode",qcol = 'Gridcodes')
        grass.run_command('v.db.addcolumn', map= 'Non_con_lake_cat_1', columns = "value INT")
        grass.run_command('v.db.update', map= 'Non_con_lake_cat_1', column = "value",qcol = 'Gridcodes')        
        grass.run_command('v.db.update', map= 'nstr_nfinalcat_F', column = "Gridcode",qcol = 'b_GC_str')
        grass.run_command('v.db.dropcolumn', map= 'nstr_nfinalcat_F', columns = ['b_GC_str'])       
#        grass.run_command('v.out.ogr', input = 'nstr_nfinalcat_F',output = os.path.join(self.tempfolder,'nstr_nfinalcat_F_Final.shp'),format= 'ESRI_Shapefile',overwrite = True)
#        Qgs.exit()       
        PERMANENT.close()
        
        
        #### 
        
        ### Calulate catchment area slope river length under projected coordinates system. 
        if projection != 'default':

            os.system('gdalwarp ' + "\"" + self.Path_dem +"\""+ "    "+ "\""+self.Path_demproj+"\""+ ' -t_srs  ' + "\""+projection+"\"")

            project = Session()
            project.open(gisdb=self.grassdb, location=self.grass_location_pro,create_opts=projection)
            grass.run_command("r.import", input = self.Path_demproj, output = 'dem_proj', overwrite = True)
            grass.run_command('g.region', raster='dem_proj')  
        
            grass.run_command('v.proj', location=self.grass_location_geo,mapset = 'PERMANENT', input = 'nstr_nfinalcat_F',overwrite = True)
            grass.run_command('v.proj', location=self.grass_location_geo,mapset = 'PERMANENT', input = 'Net_cat_F',overwrite = True)
            grass.run_command('v.proj', location=self.grass_location_geo,mapset = 'PERMANENT', input = 'Non_con_lake_cat_1',overwrite = True)
        
            grass.run_command('v.to.db', map= 'Net_cat_F',option = 'area',columns = "Area_m", units = 'meters') 
            grass.run_command('v.to.db', map= 'nstr_nfinalcat_F',option = 'length', columns = "Length_m",units = 'meters')
            grass.run_command('v.to.db', map = 'Non_con_lake_cat_1',option = 'area',columns = "Area_m", units = 'meters')
            
            grass.run_command('r.slope.aspect', elevation= 'dem_proj',slope = 'slope',aspect = 'aspect',precision = 'DCELL',overwrite = True)
            grass.run_command('r.mapcalc',expression = 'tanslopedegree = tan(slope) ',overwrite = True) 
            project.close 
        
            PERMANENT = Session()
            PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
            grass.run_command('g.region', raster='dem')
            grass.run_command('v.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'nstr_nfinalcat_F',overwrite = True)
            grass.run_command('v.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'Net_cat_F',overwrite = True) 
            grass.run_command('v.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'Non_con_lake_cat_1',overwrite = True) 
            grass.run_command('r.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'tanslopedegree',overwrite = True) 
            grass.run_command('r.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'slope',overwrite = True) 
            grass.run_command('r.proj', location=self.grass_location_pro,mapset = 'PERMANENT', input = 'aspect',overwrite = True)

        ### Add attributes to each catchments 
        # read raster arrays 
        Lake1_arr = garray.array(mapname="SelectedLakes")  
        noncnlake_arr = garray.array(mapname="Nonconnect_Lake")
        NonCL_array = garray.array(mapname="Non_con_lake_cat")
        conlake_arr = garray.array(mapname="Connect_Lake")  
        dem_array = garray.array(mapname="dem")
        width_array = garray.array(mapname="width")
        depth_array = garray.array(mapname="depth")
        obs_array = garray.array(mapname="obs")
        slope_array = garray.array(mapname="tanslopedegree")
        aspect_array = garray.array(mapname="aspect")
        slope_deg_array = garray.array(mapname="slope")
        nstr_seg_array = garray.array(mapname="nstr_seg")
        acc_array = np.absolute(garray.array(mapname="acc_grass"))
        dir_array = garray.array(mapname="ndir_Arcgis")#ndir_Arcgis
        Q_mean_array = garray.array(mapname="qmean")
        landuse_array = garray.array(mapname="landuse")
        Netcat_array = garray.array(mapname="Net_cat")
        grass.run_command('r.out.gdal', input = 'ndir_Arcgis',output = os.path.join(self.tempfolder,'ndir.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 

        grass.run_command('r.out.gdal', input = 'dir_Arcgis',output = os.path.join(self.tempfolder,'dir.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
        grass.run_command('r.out.gdal', input = 'acc_grass',output = os.path.join(self.tempfolder,'acc.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
        ## read landuse and lake infomation data 
        tempinfo = Dbf5(self.Path_allLakeply[:-3] + "dbf")
        allLakinfo = tempinfo.to_dataframe()
        landuseinfo = pd.read_csv(self.Path_Landuseinfo_in,sep=",",low_memory=False)
    ### read length and area
        sqlstat="SELECT Gridcode, Length_m FROM nstr_nfinalcat_F"
        rivleninfo = pd.read_sql_query(sqlstat, con)
        sqlstat="SELECT Gridcode, Area_m FROM Net_cat_F"
        catareainfo = pd.read_sql_query(sqlstat, con)
        sqlstat="SELECT Gridcode, Area_m FROM Non_con_lake_cat_1"
        NonConcLakeInfo = pd.read_sql_query(sqlstat, con)
        NonConcLakeInfo['SubId_riv'] = -9999
        NonConcLakeInfo['HyLakeId']  = -9999
        NonConcLakeInfo['LakeVol']   = -9999
        NonConcLakeInfo['LakeDepth'] = -9999
        NonConcLakeInfo['LakeArea']  = -9999
        NonConcLakeInfo['Laketype']  = -9999
        NonConcLakeInfo['DownLakeID'] = -9999
        NonConcLakeInfo['DA_Area'] = NonConcLakeInfo['Area_m'].values
        
        sqlstat="SELECT Obs_ID, DA_obs, STATION_NU, SRC_obs FROM obspoint"
        obsinfo = pd.read_sql_query(sqlstat, con)
        obsinfo['Obs_ID'] = obsinfo['Obs_ID'].astype(float) 
        
                  
        allcatid = np.unique(nstr_seg_array)
        allcatid = allcatid[allcatid > 0]
        catinfo2 = np.full((len(allcatid),31),-9999.00000)    
        catinfodf = pd.DataFrame(catinfo2, columns = ['SubId', "DowSubId",'RivSlope','RivLength','BasSlope','BasAspect','BasArea',
                            'BkfWidth','BkfDepth','IsLake','HyLakeId','LakeVol','LakeDepth',
                             'LakeArea','Laketype','IsObs','MeanElev','FloodP_n','Q_Mean','Ch_n','DA','Strahler','Seg_ID','Seg_order'
                             ,'Max_DEM','Min_DEM','DA_Obs','DA_error','Obs_NM','SRC_obs','NonLDArea'])
        catinfodf['Obs_NM']   =catinfodf['Obs_NM'].astype(str)         
        catinfodf['SRC_obs']  =catinfodf['SRC_obs'].astype(str)   
                   
        catinfo,NonConcLakeInfo= Generatecatinfo_riv(nstr_seg_array,acc_array,dir_array,Lake1_arr,dem_array,
             catinfodf,allcatid,width_array,depth_array,obs_array,slope_array,aspect_array,landuse_array,
             slope_deg_array,Q_mean_array,Netcat_array,landuseinfo,allLakinfo,self.nrows,self.ncols,
             rivleninfo.astype(float),catareainfo.astype(float),obsinfo,NonConcLakeInfo,NonCL_array) 
        catinfo = Streamorderanddrainagearea(catinfo)         
        catinfo.to_csv(self.Path_finalcatinfo_riv, index = None, header=True)
        NonConcLakeInfo.to_csv(self.Path_NonCLakeinfo, index = None, header=True)
        
        
        ### add lake info to selected laeks 
        grass.run_command('db.in.ogr', input=self.Path_alllakeinfoinfo,output = 'alllakeinfo',overwrite = True)
        grass.run_command('v.db.join', map= 'SelectedLakes_F',column = 'value', other_table = 'alllakeinfo',other_column ='Hylak_id', overwrite = True)
        
        grass.run_command('db.in.ogr', input=self.Path_NonCLakeinfo,output = 'result_nonlake',overwrite = True)
        grass.run_command('v.db.join', map= 'Non_con_lake_cat_1',column = 'Gridcode', other_table = 'result_nonlake',other_column ='Gridcode', overwrite = True)
        grass.run_command('v.out.ogr', input = 'Non_con_lake_cat_1',output = os.path.join(self.tempfolder,'Non_con_lake_cat_info.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')

        ### add catchment info to all river segment 
        grass.run_command('db.in.ogr', input=self.Path_finalcatinfo_riv,output = 'result_riv',overwrite = True)
        grass.run_command('v.db.join', map= 'nstr_nfinalcat_F',column = 'Gridcode', other_table = 'result_riv',other_column ='SubId', overwrite = True)
        grass.run_command('v.db.dropcolumn', map= 'nstr_nfinalcat_F', columns = ['Length_m','Gridcode'])
        grass.run_command('v.out.ogr', input = 'nstr_nfinalcat_F',output = os.path.join(self.tempfolder,'finalriv_info1.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
        
        grass.run_command('v.db.join', map= 'Net_cat_F',column = 'Gridcode', other_table = 'result_riv',other_column ='SubId', overwrite = True)
        grass.run_command('v.db.dropcolumn', map= 'Net_cat_F', columns = ['Area_m','Gridcode','GC_str'])
        grass.run_command('v.out.ogr', input = 'Net_cat_F',output = os.path.join(self.tempfolder,'finalriv_catinfo1.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')

        PERMANENT.close()
                
        processing.run("native:dissolve", {'INPUT':os.path.join(self.tempfolder,'finalriv_catinfo1.shp'),'FIELD':['SubId'],'OUTPUT':os.path.join(self.tempfolder,'finalriv_catinfo_dis.shp')},context = context)
        processing.run("native:dissolve", {'INPUT':os.path.join(self.tempfolder,'finalriv_info1.shp'),'FIELD':['SubId'],'OUTPUT':os.path.join(self.tempfolder,'finalriv_info_dis.shp')},context = context)
        
        processing.run("native:centroids", {'INPUT':os.path.join(self.tempfolder,'finalriv_catinfo_dis.shp'),'ALL_PARTS':False,'OUTPUT':os.path.join(self.tempfolder,'Centerpoints.shp')},context = context)
        processing.run("native:addxyfields", {'INPUT':os.path.join(self.tempfolder,'Centerpoints.shp'),'CRS':QgsCoordinateReferenceSystem(self.SpRef_in),'OUTPUT':os.path.join(self.tempfolder,'ctwithxy.shp')},context = context)
        processing.run("native:joinattributestable",{'INPUT':os.path.join(self.tempfolder,'finalriv_catinfo_dis.shp'),'FIELD':'SubId','INPUT_2':os.path.join(self.tempfolder,'ctwithxy.shp'),'FIELD_2':'SubId',
                          'FIELDS_TO_COPY':['x','y'],'METHOD':0,'DISCARD_NONMATCHING':False,'PREFIX':'centroid_','OUTPUT':os.path.join(self.OutputFolder,'finalriv_info_ply.shp')},context = context)
        processing.run("native:joinattributestable",{'INPUT':os.path.join(self.tempfolder,'finalriv_info_dis.shp'),'FIELD':'SubId','INPUT_2':os.path.join(self.tempfolder,'ctwithxy.shp'),'FIELD_2':'SubId',
                          'FIELDS_TO_COPY':['x','y'],'METHOD':0,'DISCARD_NONMATCHING':False,'PREFIX':'centroid_','OUTPUT':os.path.join(self.OutputFolder,'finalriv_info.shp')},context = context)
        
        
        ##### export  lakes
        ## Non_Connected Lakes
        NonConcLakeInfo = NonConcLakeInfo[NonConcLakeInfo['SubId_riv'] > 0]
        NonCL_Lakeids = NonConcLakeInfo['Gridcode'].values
        NonCL_Lakeids = np.unique(NonCL_Lakeids)
        NonCL_Lakeids = NonCL_Lakeids[NonCL_Lakeids > 0]
        
        exp ='Hylak_id' + '  IN  (  ' +  str(int(NonCL_Lakeids[0]))      
        for i in range(1,len(NonCL_Lakeids)):
            exp = exp + " , "+str(int(NonCL_Lakeids[i]))        
        exp = exp + ')'
        processing.run("native:extractbyexpression", {'INPUT':self.Path_allLakeply,'EXPRESSION':exp,'OUTPUT':os.path.join(self.OutputFolder, 'Non_Con_Lake_Ply.shp')})


        exp ='Gridcode' + '  IN  (  ' +  str(int(NonCL_Lakeids[0]))      
        for i in range(1,len(NonCL_Lakeids)):
            exp = exp + " , "+str(int(NonCL_Lakeids[i]))        
        exp = exp + ')'
        processing.run("native:extractbyexpression", {'INPUT':os.path.join(self.tempfolder,'Non_con_lake_cat_info.shp'),'EXPRESSION':exp,'OUTPUT':os.path.join(self.OutputFolder,'Non_con_lake_cat_info.shp')})
        
        
        ### Non_Connected Lakes
        CL_Lakeids = catinfo['HyLakeId'].values
        CL_Lakeids = np.unique(CL_Lakeids)
        CL_Lakeids = CL_Lakeids[CL_Lakeids > 0]
        exp ='Hylak_id' + '  IN  (  ' +  str(int(CL_Lakeids[0]))      
        for i in range(1,len(CL_Lakeids)):
            exp = exp + " , "+str(int(CL_Lakeids[i]))        
        exp = exp + ')'
        processing.run("native:extractbyexpression", {'INPUT':self.Path_allLakeply,'EXPRESSION':exp,'OUTPUT':os.path.join(self.OutputFolder, 'Con_Lake_Ply.shp')})
        
        Qgs.exit()  

        
    def RoutingNetworkTopologyUpdateToolset_cat_default(self,projection = 'default'):
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session
        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
        
        Netcat_array = garray.array(mapname="Net_cat")
        catinfo = pd.read_csv(self.Path_finalcatinfo_riv)
        NonConcLakeInfo= pd.read_csv(self.Path_NonCLakeinfo)
        NonConcLakeInfo['SubId_cat'] = -9999
        
        lakeids = catinfo['HyLakeId'].values
        lakeids = np.unique(lakeids)
        lakeids = lakeids[lakeids > 0]
        maxsubid = np.max(catinfo['SubId'].values)
        for i in range(0,len(lakeids)):
            lakeid = lakeids[i]
            nlakecatid = maxsubid + 10 + i
            lakeindex = catinfo['HyLakeId'] == lakeid   #### rows in catinfo that has this lake id
            row_lakeindex = np.argwhere(catinfo['HyLakeId'].values == lakeid)   #### rows in catinfo that has this lake id 
            lakecat_info = catinfo[lakeindex]
            
            catscoverlays = catinfo[lakeindex]['SubId'].values  #### catid that covered by this lake
            catscoverlays = np.unique(catscoverlays) 
            #### find lake catchment downstream id 
            downcatid = -1 
            for icat in range(0,len(lakecat_info)):
                if len(np.argwhere(catscoverlays == lakecat_info['DowSubId'].values[icat])) <= 0:
                    downcatid = lakecat_info['DowSubId'].values[icat]
    
            ### change upstream catchment's downstream catchment id  to newlakecatid 
            upcats_to_lakecats = catinfo['DowSubId'].isin(catscoverlays)#### rows that downstream id is the lake catchment
            upcatrows  = np.argwhere(upcats_to_lakecats == True)
            for iupcatrowidx in range(0,len(upcatrows)):
                iupcatrow =upcatrows[iupcatrowidx] 
                catinfo.loc[iupcatrow,'DowSubId'] = nlakecatid
                
            ####change non connected lake's subbasin id 
            NoncLakes_inthesecats = NonConcLakeInfo['SubId_riv'].isin(catscoverlays)
            NoncLakes_rows  = np.argwhere(NoncLakes_inthesecats == True)
            for nonclrowidx in range(0,len(NoncLakes_rows)):
                nonclrow =NoncLakes_rows[nonclrowidx] 
                NonConcLakeInfo.loc[nonclrow,'SubId_cat'] = nlakecatid
                            

            ## change catchment raster merge lake catchments 
            for j in range(0,len(lakecat_info)):
                jcat = lakecat_info['SubId'].values[j]
                jcatrowcol = np.argwhere(Netcat_array==jcat).astype(int)
                Netcat_array[jcatrowcol[:,0],jcatrowcol[:,1]] = nlakecatid

            #### change lake catchemt infomation 
            for irow in range(0,len(row_lakeindex)):
                irow_idx = row_lakeindex[irow]
                catinfo.loc[irow_idx,'SubId'] = nlakecatid
                catinfo.loc[lakeindex,'DowSubId'] = downcatid
                catinfo.loc[lakeindex,'RivSlope'] = np.average(lakecat_info['RivSlope'].values,weights = lakecat_info['RivLength'].values)
                catinfo.loc[lakeindex,'RivLength'] =0.0
                catinfo.loc[lakeindex,'BasAspect'] = np.average(lakecat_info['BasAspect'].values,weights = lakecat_info['BasArea'].values)
                catinfo.loc[lakeindex,'BasArea'] = np.sum(lakecat_info['BasArea'].values)
                catinfo.loc[lakeindex,'BasSlope'] = np.average(lakecat_info['BasSlope'].values,weights = lakecat_info['BasArea'].values)
                catinfo.loc[lakeindex,'BkfWidth'] = np.max(lakecat_info['BkfWidth'].values)
                catinfo.loc[lakeindex,'BkfDepth'] = np.max(lakecat_info['BkfDepth'].values)
                catinfo.loc[lakeindex,'MeanElev'] = np.average(lakecat_info['MeanElev'].values,weights = lakecat_info['BasArea'].values)
                if len(lakecat_info.loc[lakecat_info['IsObs'] > 0]) > 0:
                    catinfo.loc[lakeindex,'IsObs'] = lakecat_info.loc[lakecat_info['IsObs'] > 0]['IsObs'].values[0]
                    catinfo.loc[lakeindex,'Obs_NM'] = lakecat_info.loc[lakecat_info['IsObs'] > 0]['Obs_NM'].values[0]
                    catinfo.loc[lakeindex,'DA_Obs'] = lakecat_info.loc[lakecat_info['IsObs'] > 0]['DA_Obs'].values[0]
                    catinfo.loc[lakeindex,'DA_error'] = lakecat_info.loc[lakecat_info['IsObs'] > 0]['DA_error'].values[0]
                    catinfo.loc[lakeindex,'SRC_obs'] = lakecat_info.loc[lakecat_info['IsObs'] > 0]['SRC_obs'].values[0]
                else:
                    catinfo.loc[lakeindex,'IsObs'] = -1.2345
                    catinfo.loc[lakeindex,'Obs_NM'] = -1.2345
                    catinfo.loc[lakeindex,'DA_Obs'] = -1.2345
                    catinfo.loc[lakeindex,'DA_error'] = -1.2345  
                    catinfo.loc[lakeindex,'SRC_obs'] = -1.2345                    
                catinfo.loc[lakeindex,'FloodP_n'] = np.average(lakecat_info['FloodP_n'].values,weights = lakecat_info['RivLength'].values)
                catinfo.loc[lakeindex,'Q_Mean'] = np.average(lakecat_info['Q_Mean'].values,weights = lakecat_info['RivLength'].values)
                catinfo.loc[lakeindex,'Ch_n'] = np.average(lakecat_info['Ch_n'].values,weights = lakecat_info['RivLength'].values)
                catinfo.loc[lakeindex,'DA'] = np.max(lakecat_info['DA'].values)
                catinfo.loc[lakeindex,'Strahler'] =np.max(lakecat_info['Strahler'].values)
                catinfo.loc[lakeindex,'Seg_ID'] =  -1.2345
                catinfo.loc[lakeindex,'Seg_order'] =  np.max(lakecat_info['Seg_order'].values)
                catinfo.loc[lakeindex,'Max_DEM'] = np.max(lakecat_info['Max_DEM'].values)
                catinfo.loc[lakeindex,'Min_DEM'] = np.min(lakecat_info['Min_DEM'].values)
        
        for inoncl in range(0,len(NonConcLakeInfo)):
            if NonConcLakeInfo['SubId_cat'].values[inoncl] < 0:
                NonConcLakeInfo.loc[inoncl,'SubId_cat'] = NonConcLakeInfo['SubId_riv'].values[inoncl]
    
        catinfo.to_csv(self.Path_finalcatinfo_cat, index = None, header=True)
        NonConcLakeInfo.to_csv(self.Path_NonCLakeinfo_cat, index = None, header=True)
        temparray = garray.array()
        temparray[:,:] = -9999
        temparray[:,:] = Netcat_array[:,:]
        temparray.write(mapname="Cat_Lake_combined", overwrite=True)  #### write new stream id to a grass raster 
        grass.run_command('r.null', map='Cat_Lake_combined',setnull=-9999)
        grass.run_command('r.mapcalc',expression = 'Cat_Lake_combined = int(Cat_Lake_combined)',overwrite = True)
        grass.run_command('r.out.gdal', input = 'Cat_Lake_combined',output = os.path.join(self.tempfolder,'Cat_Lake_combined.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture') 
        grass.run_command('r.to.vect', input='Cat_Lake_combined',output='Cat_Lake_combined_F1',type='area', overwrite = True)    ## save to vector 
        grass.run_command('v.db.addcolumn', map= 'Cat_Lake_combined_F1', columns = "GC VARCHAR(40)")     
        grass.run_command('v.db.update', map= 'Cat_Lake_combined_F1', column = "GC",qcol = 'value')
        grass.run_command('v.dissolve', input= 'Cat_Lake_combined_F1', column = "GC",output = 'Cat_Lake_combined_F',overwrite = True) ### dissolve based on gridcode
        grass.run_command('v.db.addcolumn', map= 'Cat_Lake_combined_F', columns = "Gridcode INT")        
        grass.run_command('v.db.update', map= 'Cat_Lake_combined_F', column = "Gridcode",qcol = 'GC')
        
        
        grass.run_command('db.in.ogr', input=self.Path_NonCLakeinfo_cat,output = 'result_nonlake2',overwrite = True)
        grass.run_command('v.db.join', map= 'Non_con_lake_cat_1',column = 'Gridcode', other_table = 'result_nonlake2',other_column ='Gridcode', subset_columns = 'SubId_cat', overwrite = True)
        grass.run_command('v.out.ogr', input = 'Non_con_lake_cat_1',output = os.path.join(self.OutputFolder,'Non_con_lake.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')

        
        
        grass.run_command('db.in.ogr', input=self.Path_finalcatinfo_cat,output = 'result_cat',overwrite = True)
        grass.run_command('v.db.join', map= 'Cat_Lake_combined_F',column = 'Gridcode', other_table = 'result_cat',other_column ='SubId', overwrite = True)
        grass.run_command('v.db.dropcolumn', map= 'Net_cat_F', columns = ['Area_m','Gridcode','GC'])
        grass.run_command('v.out.ogr', input = 'Cat_Lake_combined_F',output = os.path.join(self.tempfolder,'finalcatinfo_cat1.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
        PERMANENT.close()
        
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing
        from processing.tools import dataobjects
           
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        context = dataobjects.createContext()
        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)
        processing.run("native:centroids", {'INPUT':os.path.join(self.tempfolder,'finalcatinfo_cat1.shp'),'ALL_PARTS':False,'OUTPUT':os.path.join(self.tempfolder,'Centerpoints2.shp')},context = context)
        processing.run("native:addxyfields", {'INPUT':os.path.join(self.tempfolder,'Centerpoints2.shp'),'CRS':QgsCoordinateReferenceSystem(self.SpRef_in),'OUTPUT':os.path.join(self.tempfolder,'ctwithxy2.shp')},context = context)
        processing.run("native:joinattributestable",{'INPUT':os.path.join(self.tempfolder,'finalcatinfo_cat1.shp'),'FIELD':'SubId','INPUT_2':os.path.join(self.tempfolder,'ctwithxy2.shp'),'FIELD_2':'SubId',
                          'FIELDS_TO_COPY':['x','y'],'METHOD':0,'DISCARD_NONMATCHING':False,'PREFIX':'centroid_','OUTPUT':os.path.join(self.OutputFolder,'finalcat_info.shp')},context = context)
        Qgs.exit()
             
###########################################################################3
    def Output_Clean(self,Out = 'Simple',clean = 'True'):
        
        import grass.script as grass
        from grass.script import array as garray
        import grass.script.setup as gsetup
        from grass.pygrass.modules.shortcuts import general as g
        from grass.pygrass.modules.shortcuts import raster as r
        from grass.pygrass.modules import Module
        from grass_session import Session

        os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD',GRASS_VERBOSE='-1'))
        PERMANENT = Session()
        PERMANENT.open(gisdb=self.grassdb, location=self.grass_location_geo,create_opts=self.SpRef_in)
                
        if Out == 'Simple':
#            grass.run_command('v.out.ogr', input = 'finalcat_F',output = os.path.join(self.OutputFolder,'finalcat_info.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
#            grass.run_command('v.out.ogr', input = 'SelectedLakes_F',output = os.path.join(self.OutputFolder,'SelectedLakes.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
#            grass.run_command('v.out.ogr', input = 'nstr_nfinalcat_F',output = os.path.join(self.OutputFolder,'finalriv_info.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
#            grass.run_command('v.out.ogr', input = 'Cat_Lake_combined_F',output = os.path.join(self.OutputFolder,'finalcat_info.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
#            grass.run_command('v.out.ogr', input = 'Hylake',output = os.path.join(self.OutputFolder,'AllLakes.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
            grass.run_command('v.out.ogr', input = 'obspoint',output = os.path.join(self.OutputFolder,'obspoint.shp'),format= 'ESRI_Shapefile',overwrite = True,quiet = 'Ture')
        if Out == 'All':
            grass.run_command('r.out.gdal', input = 'SelectedLakes',output = os.path.join(self.OutputFolder,'SelectedLakes.tif'),format= 'GTiff',overwrite = True,quiet = 'Ture')
        if clean == 'True':
            shutil.rmtree(self.grassdb,ignore_errors=True)
            shutil.rmtree(self.tempfolder,ignore_errors=True)


    def GenerateRavenInput(self,Finalcat_NM = 'finalcat_info',lenThres = 1,iscalmanningn = -1,Nonconnectlake = -1,NonconLakeinfo = 'Non_con_lake',Startyear = 1980,EndYear = 2010
                          ,CA_HYDAT = '#',WarmUp = 0,Template_Folder = '#'):
        
        
        Model_Folder     = os.path.join(self.OutputFolder,'Model')
        Raveinputsfolder = os.path.join(Model_Folder,'RavenInput')
        Obs_Folder       = os.path.join(Raveinputsfolder,'obs')
        
        if not os.path.exists(Model_Folder):
            os.makedirs(Model_Folder)
        
        ### check if there is a model input template provided    
        if Template_Folder != '#':
            fromDirectory = Template_Folder
            toDirectory   = Model_Folder
            copy_tree(fromDirectory, toDirectory)
            
        if not os.path.exists(Raveinputsfolder):
            os.makedirs(Raveinputsfolder)
        if not os.path.exists(Obs_Folder):
            os.makedirs(Obs_Folder)
            


        finalcatchpath = os.path.join(self.OutputFolder,Finalcat_NM)
        
        tempinfo = Dbf5(finalcatchpath + ".dbf")#np.genfromtxt(hyinfocsv,delimiter=',')
        ncatinfo = tempinfo.to_dataframe()
        ncatinfo2 = ncatinfo.drop_duplicates('SubId', keep='first')
        ncatinfo2 = ncatinfo2[ncatinfo2['SubId'] > 0]
        
        if Nonconnectlake > 0:
            NonconLakeinfopath = os.path.join(self.OutputFolder,NonconLakeinfo)
            tempinfo = Dbf5( NonconLakeinfopath + ".dbf")#np.genfromtxt(hyinfocsv,delimiter=',')
            nclakeinfo = tempinfo.to_dataframe()
            nclakeinfo['Gridcode'] = nclakeinfo['Gridcode'].astype(float)
            nclakeinfo['SubId_cat'] = nclakeinfo['SubId_cat'].astype(float)
            nclakeinfo['Area_m'] = nclakeinfo['Area_m'].astype(float)
            nclakeinfo = nclakeinfo.drop_duplicates('Gridcode', keep='first')
        else:
            nclakeinfo = pd.DataFrame(np.full((1,4),-9999), columns = ['Gridcode', "SubId_riv","SubId_cat","Area_m"])
#            print(nclakeinfo)
            
        nclakeinfo = Writervhchanl(ncatinfo2,Raveinputsfolder,lenThres,iscalmanningn,nclakeinfo)
        writelake(ncatinfo2,Raveinputsfolder,nclakeinfo)
        nclakeinfo.to_csv(os.path.join(self.OutputFolder,'Non_connect_Lake_routing_info.csv'),index = None, header=True)
        WriteObsfiles(ncatinfo2,Raveinputsfolder,Obs_Folder,Startyear + WarmUp,EndYear,CA_HYDAT,Template_Folder)
        
        

    def Locate_subid_needsbyuser(self,Path_Points = '#',Guage_Col_Name = 'Obs_NM',Guage_NMS = '#',subid_col_Name='SubId',Path_products='#'):
        # obtain subbasin ID based on either points or guage names
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        
        SubId_Selected = -1 
        if Guage_NMS[0] != '#':
            hyinfocsv  = Path_products[:-3] + "dbf"
            tempinfo   = Dbf5(hyinfocsv)
            hyshdinfo2 = tempinfo.to_dataframe().drop_duplicates(subid_col_Name, keep='first')
            hyshdinfo2 = hyshdinfo2.loc[hyshdinfo2[Guage_Col_Name] != '-9999.0']
            hyshdinfo2 = hyshdinfo2.loc[hyshdinfo2[Guage_Col_Name].isin(Guage_NMS)]
            hyshdinfo2 = hyshdinfo2[[Guage_Col_Name,subid_col_Name]]
            hyshdinfo2.to_csv(os.path.join(self.OutputFolder,'SubIds_Selected.csv'),sep=',', index = None)
            SubId_Selected = hyshdinfo2[subid_col_Name].values
        
        if Path_Points != '#':
            processing.run("saga:addpolygonattributestopoints", {'INPUT':Path_Points,'POLYGONS':Path_products,'FIELDS':subid_col_Name,'OUTPUT':os.path.join(self.tempfolder,'Sub_Selected_by_Points.shp')})
            hyinfocsv  = os.path.join(self.tempfolder,'Sub_Selected_by_Points.shp')[:-3] + "dbf"
            tempinfo   = Dbf5(hyinfocsv)
            hyshdinfo2 = tempinfo.to_dataframe()
            hyshdinfo2.to_csv(os.path.join(self.OutputFolder,'SubIds_Selected.csv'),sep=',', index = None)
            SubId_Selected = hyshdinfo2[subid_col_Name].values
        
        Qgs.exit()
        
#        self.selectfeaturebasedonID(Path_shpfile = Path_products,sub_colnm = 'SubId',down_colnm = 'DowSubId',mostdownid = SubId_Selected,mostupstreamid = np.full(len(SubId_Selected),-1),OutBaseName='finalcat_info')
        return SubId_Selected
            
    def Select_Routing_product_based_SubId(self,Path_final_riv = '#',Path_final_riv_ply = '#',Path_Con_Lake_ply = '#',Path_NonCon_Lake_ply='#',Path_NonClakeinfo = '#',
                                           sub_colnm = 'SubId',down_colnm = 'DowSubId',mostdownid = [-1],mostupstreamid = [-1]):

        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

        ##3
        hyinfocsv = Path_final_riv[:-3] + "dbf"
        tempinfo = Dbf5(hyinfocsv)
        hyshdinfo2 = tempinfo.to_dataframe().drop_duplicates(sub_colnm, keep='first')
        hyshdinfo =  hyshdinfo2[[sub_colnm,down_colnm]].astype('float').values
                        
        #### incase  mostupstreamid did not procvided, automatically assgin -1 
        if mostupstreamid[0] == -1:
            mostupstreamid = np.full(len(mostdownid),-1) 
            
        ### Loop for each downstream id     
        for isub in range(0,len(mostdownid)):
            OutHyID  = mostdownid[isub]
            OutHyID2 = mostupstreamid[isub]
                                    
            OutputFolder_isub = os.path.join(self.OutputFolder,'SubId_'+str(OutHyID))
            if not os.path.exists(OutputFolder_isub):
                os.makedirs(OutputFolder_isub)
                
            ## find all subid control by this subid 
            HydroBasins1 = Defcat(hyshdinfo,OutHyID) ### return fid of polygons that needs to be select 
            if OutHyID2 > 0:
                HydroBasins2 = Defcat(hyshdinfo,OutHyID2)            
                ###  exculde the Ids in HydroBasins2 from HydroBasins1
                for i in range(len(HydroBasins2)):
                    rows =np.argwhere(HydroBasins1 == HydroBasins2[i])
                    HydroBasins1 = np.delete(HydroBasins1, rows)
                HydroBasins = HydroBasins1            
            else:
                HydroBasins = HydroBasins1
            
            #### extract final_riv and final_riv_ply
            Outputfilename_riv = os.path.join(OutputFolder_isub,os.path.basename(Path_final_riv))
            Selectfeatureattributes(processing,Input = Path_final_riv,Output=Outputfilename_riv,Attri_NM = 'SubId',Values = HydroBasins)
            Outputfilename_riv_ply = os.path.join(OutputFolder_isub,os.path.basename(Path_final_riv_ply))
            Selectfeatureattributes(processing,Input = Path_final_riv_ply,Output=Outputfilename_riv_ply,Attri_NM = 'SubId',Values = HydroBasins)
            
            #### extract lakes 
            finalcat_csv     = Outputfilename_riv_ply[:-3] + "dbf"
            finalcat_info    = Dbf5(finalcat_csv)
            finalcat_info    = finalcat_info.to_dataframe().drop_duplicates('SubId', keep='first')
        
            Connect_Lakeids  = np.unique(finalcat_info['HyLakeId'].values)
            Connect_Lakeids  = Connect_Lakeids[Connect_Lakeids > 0]
        
            SubIds           = np.unique(finalcat_info['SubId'].values)
            SubIds           = SubIds[SubIds > 0]        
        
            NonConn_Lakes    = Path_NonClakeinfo[:-3] + "dbf"
            NonConn_Lakes    = Dbf5(NonConn_Lakes)
            NonConn_Lakes    = NonConn_Lakes.to_dataframe()
            NonConn_Lakes['SubId_riv']    = pd.to_numeric(NonConn_Lakes['SubId_riv'], downcast='float')
        
            NonConn_Lakes_p  = NonConn_Lakes.loc[NonConn_Lakes['SubId_riv'].isin(SubIds)]
            NonCL_Lakeids    = NonConn_Lakes_p['value'].values
        
        
            Selectfeatureattributes(processing,Input = Path_Con_Lake_ply,Output=os.path.join(OutputFolder_isub,os.path.basename(Path_Con_Lake_ply)),Attri_NM = 'Hylak_id',Values = Connect_Lakeids)
        
            Selectfeatureattributes(processing,Input = Path_NonCon_Lake_ply,Output=os.path.join(OutputFolder_isub,os.path.basename(Path_NonCon_Lake_ply)),Attri_NM = 'Hylak_id',Values = NonCL_Lakeids)
        
            Selectfeatureattributes(processing,Input = Path_NonClakeinfo,Output=os.path.join(OutputFolder_isub,os.path.basename(Path_NonClakeinfo)),Attri_NM = 'value',Values = NonCL_Lakeids)
        
        Qgs.exit()
    
    def Customize_Routing_Topology(self,DataFolder = '#',finalrvi_ply_NM = 'finalriv_info_ply.shp',Non_ConnL_Cat_NM = 'Non_con_lake_cat_info.shp',Non_ConnL_ply_NM='Non_Con_Lake_Ply.shp',
                                   ConnL_ply_NM='Con_Lake_Ply.shp',finalriv_NM = 'finalriv_info.shp', Area_Min = -1,
                                   sub_colnm='SubId',down_colnm='DowSubId',DA_colnm = 'DA',SegID_colnm = 'Seg_ID'):
        #### generate river catchments based on minmum area.
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing
        from processing.tools import dataobjects
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        context = dataobjects.createContext()
        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)
                
        
        # obtain river segment ids based on area thresthold 
        
        Path_final_rviply = os.path.join(DataFolder,finalrvi_ply_NM)
        Path_final_riv    = os.path.join(DataFolder,finalriv_NM)
        Path_Non_ConL_cat = os.path.join(DataFolder,Non_ConnL_Cat_NM)
        Path_Conl_ply     = os.path.join(DataFolder,ConnL_ply_NM)
        Path_Non_ConL_ply = os.path.join(DataFolder,Non_ConnL_ply_NM)            
        
        
        finalriv_csv     = Path_final_rviply[:-3] + "dbf"
        finalriv_info    = Dbf5(finalriv_csv)
        finalriv_info    = finalriv_info.to_dataframe().drop_duplicates(sub_colnm, keep='first')
              
        Selected_riv     = finalriv_info.loc[finalriv_info[DA_colnm] >= Area_Min*1000*1000] # find river with drainage area larger than area thresthold 
        Selected_riv     = UpdateTopology(Selected_riv,UpdateSubId = -1)
        Selected_riv     = Selected_riv.sort_values(["Strahler"], ascending = (True))   ###sort selected river by Strahler stream order 
        Subid_main       = Selected_riv[sub_colnm].values
                
        mapoldnew_info   = finalriv_info.copy(deep = True)
        mapoldnew_info['nsubid']     = -1
        mapoldnew_info['ndownsubid'] = -1   
        mapoldnew_info.reset_index(drop=True, inplace=True)  
        
        ### Obtain connected lakes based on current river segment
        Connected_Lake_Mainriv = Selected_riv['HyLakeId'].values
        Connected_Lake_Mainriv = np.unique(Connected_Lake_Mainriv[Connected_Lake_Mainriv>0])
        Lakecover_riv          = finalriv_info.loc[finalriv_info['HyLakeId'].isin(Connected_Lake_Mainriv)]
        Subid_lakes            = Lakecover_riv[sub_colnm].values
        
        
        #####
        NonConn_Lakes               = Path_Non_ConL_cat[:-3] + "dbf"
        NonConn_Lakes               = Dbf5(NonConn_Lakes)
        NonConn_Lakes               = NonConn_Lakes.to_dataframe()
        NonConn_Lakes['SubId_riv']  = pd.to_numeric(NonConn_Lakes['SubId_riv'], downcast='float')
        NonConn_Lakes['DownLakeID'] = pd.to_numeric(NonConn_Lakes['DownLakeID'], downcast='float')
        NonConn_Lakes['Area_m']     = pd.to_numeric(NonConn_Lakes['Area_m'], downcast='float')
        NonConn_Lakes['DA_Area']    = pd.to_numeric(NonConn_Lakes['DA_Area'], downcast='float')
        
        New_NonConn_Lakes           = NonConn_Lakes.copy(deep = True)
        
        Conn_Lakes_ply              = Path_Conl_ply[:-3] + "dbf"
        Conn_Lakes_ply              = Dbf5(Conn_Lakes_ply)
        Conn_Lakes_ply              = Conn_Lakes_ply.to_dataframe()
        Conn_Lakes_ply              = Conn_Lakes_ply.drop_duplicates('Hylak_id', keep='first')
        All_Conn_Lakeids            = Conn_Lakes_ply['Hylak_id'].values        
        mask                        = np.in1d(All_Conn_Lakeids, Connected_Lake_Mainriv)
        Conn_To_NonConlakeids       = All_Conn_Lakeids[np.logical_not(mask)]
                
        
#        Select_SubId_Link= Select_SubId_Lakes(ConLakeId,finalriv_info,Subid_main)
        ### obtain rivsegments that covered by remaining lakes  
        Selected_riv_ids  = np.unique(Subid_main) #np.unique(np.concatenate([Subid_main,Subid_lakes]))
        routing_info      = finalriv_info[['SubId','DowSubId']].astype('float').values
        Seg_IDS           = Selected_riv['Seg_ID'].values
        Seg_IDS           = np.unique(Seg_IDS)
        for iseg in range(0,len(Seg_IDS)):
#            print('#########################################################################################33333')
            i_seg_id        = Seg_IDS[iseg]
            i_seg_info      = Selected_riv[Selected_riv['Seg_ID'] == i_seg_id]                      
            i_seg_info      = i_seg_info.sort_values(["Seg_order"], ascending = (True))

            modifysubids = []
            seg_order = 1
            for iorder in range(0,len(i_seg_info)):
                tsubid              = i_seg_info['SubId'].values[iorder]
                iorder_Lakeid       = i_seg_info['HyLakeId'].values[iorder]
                modifysubids.append(tsubid)
                processed_subid = np.unique(mapoldnew_info.loc[mapoldnew_info['nsubid'] > 0][sub_colnm].values)
                
                ### two seg has the same HyLakeId id, can be merged 
                if iorder == len(i_seg_info) - 1:
                    seg_sub_ids   = np.asarray(modifysubids)
                    ## if needs to add lake sub around the main stream 
                    if iorder_Lakeid > 0:
                        subid_cur_lake_info        = finalriv_info.loc[finalriv_info['HyLakeId'] ==iorder_Lakeid]  
                        routing_info_lake          = subid_cur_lake_info[['SubId','DowSubId']].astype('float').values
                        UpstreamLakeids            = Defcat(routing_info_lake,tsubid)     
                        seg_sub_ids                = np.unique(np.concatenate([seg_sub_ids,UpstreamLakeids]))
                        seg_sub_ids                = seg_sub_ids[seg_sub_ids>0]
                        
                    ### merge all subbasin not connected to the main river but drainarge to this tsubid     
                    All_up_subids         = Defcat(routing_info,tsubid)
                    All_up_subids         = All_up_subids[All_up_subids > 0]
                    mask1                 = np.in1d(All_up_subids, Subid_main)  ### exluced ids that belongs to main river stream 
                    All_up_subids_no_main = All_up_subids[np.logical_not(mask1)]
                    
                    seg_sub_ids   = np.unique(np.concatenate([seg_sub_ids,All_up_subids_no_main]))
                    seg_sub_ids   = seg_sub_ids[seg_sub_ids>0]            
                    mask          = np.in1d(seg_sub_ids, processed_subid)
                    seg_sub_ids   = seg_sub_ids[np.logical_not(mask)]
                    mapoldnew_info = New_SubId_To_Dissolve(subid = tsubid,catchmentinfo = finalriv_info,mapoldnew_info = mapoldnew_info,ismodifids = 1,modifiidin = seg_sub_ids,mainriv = Selected_riv,Islake = 2,seg_order = seg_order) 
                    New_NonConn_Lakes = ConnectLake_to_NonConnectLake_Updateinfo(NonC_Lakeinfo = New_NonConn_Lakes,finalriv_info = finalriv_info ,Merged_subids = seg_sub_ids,Connect_Lake_ply_info = Conn_Lakes_ply,ConLakeId = iorder_Lakeid)
                    modifysubids   = []
                    seg_order      = seg_order + 1
                        
                elif i_seg_info['HyLakeId'].values[iorder] == i_seg_info['HyLakeId'].values[iorder + 1] and i_seg_info['IsObs'].values[iorder] < 0:
                    continue                    
                else: 
                    seg_sub_ids    = np.asarray(modifysubids)
                    ## if needs to add lake sub around the main stream                     
                    if iorder_Lakeid > 0:
                        subid_cur_lake_info        = finalriv_info.loc[finalriv_info['HyLakeId'] ==iorder_Lakeid]  
                        routing_info_lake          = subid_cur_lake_info[['SubId','DowSubId']].astype('float').values
                        UpstreamLakeids            = Defcat(routing_info_lake,tsubid)     
                        seg_sub_ids                = np.unique(np.concatenate([seg_sub_ids,UpstreamLakeids]))
                        seg_sub_ids                = seg_sub_ids[seg_sub_ids>0]
                        
                    All_up_subids         = Defcat(routing_info,tsubid)
                    All_up_subids         = All_up_subids[All_up_subids > 0]
                    mask1                 = np.in1d(All_up_subids, Subid_main)  ### exluced ids that belongs to main river stream 
                    All_up_subids_no_main = All_up_subids[np.logical_not(mask1)]
                    
                    seg_sub_ids   = np.unique(np.concatenate([seg_sub_ids,All_up_subids_no_main]))
                    seg_sub_ids   = seg_sub_ids[seg_sub_ids>0]            
                    mask          = np.in1d(seg_sub_ids, processed_subid)
                    seg_sub_ids   = seg_sub_ids[np.logical_not(mask)]
                    mapoldnew_info = New_SubId_To_Dissolve(subid = tsubid,catchmentinfo = finalriv_info,mapoldnew_info = mapoldnew_info,ismodifids = 1,modifiidin = seg_sub_ids,mainriv = Selected_riv,Islake = 2,seg_order = seg_order) 
                    New_NonConn_Lakes = ConnectLake_to_NonConnectLake_Updateinfo(NonC_Lakeinfo = New_NonConn_Lakes,finalriv_info = finalriv_info ,Merged_subids = seg_sub_ids,Connect_Lake_ply_info = Conn_Lakes_ply,ConLakeId = iorder_Lakeid)
                    modifysubids   = []
                    seg_order      = seg_order + 1

        ######################################################################
        ## select river based on river ids 
        Path_Temp_final_rviply = os.path.join(self.tempfolder,'temp1_finalriv_ply.shp')
        Path_Temp_final_rvi    = os.path.join(self.tempfolder,'temp1_finalriv.shp')
        
            
        Selectfeatureattributes(processing,Input = Path_final_riv    ,Output=Path_Temp_final_rvi    ,Attri_NM = 'SubId',Values = Selected_riv_ids)
        processing.run("native:dissolve", {'INPUT':Path_final_rviply,'FIELD':['SubId'],'OUTPUT':Path_Temp_final_rviply},context = context)   
         
        UpdateTopology(mapoldnew_info)        
        Modify_Feature_info(Path_Temp_final_rviply,mapoldnew_info)
        Modify_Feature_info(Path_Temp_final_rvi,mapoldnew_info)

        ######################################################################################################3
        ## create output folder
        outputfolder_subid = os.path.join(self.OutputFolder,'SubArea' +'_'+str(Area_Min))
        if not os.path.exists(outputfolder_subid):
            os.makedirs(outputfolder_subid)
        
        #### export lake polygons 
        
        New_NonConn_Lakes.to_csv(os.path.join(outputfolder_subid,'Non_connect_lake_info.csv'),index='False')            
        Selectfeatureattributes(processing,Input =Path_Conl_ply ,Output=os.path.join(outputfolder_subid,ConnL_ply_NM),Attri_NM = 'Hylak_id',Values = Connected_Lake_Mainriv)        
        Selectfeatureattributes(processing,Input =Path_Non_ConL_ply ,Output=os.path.join(outputfolder_subid,Non_ConnL_ply_NM),Attri_NM = 'Hylak_id',Values = New_NonConn_Lakes['value'].values)
        
        ###
        
        Copyfeature_to_another_shp_by_attribute(Source_shp = Path_Conl_ply,Target_shp =os.path.join(outputfolder_subid,Non_ConnL_ply_NM),Col_NM='Hylak_id',Values=Conn_To_NonConlakeids,Attributes = Conn_Lakes_ply)
        
        
        Path_out_final_rviply = os.path.join(outputfolder_subid,'finalriv_ply.shp')
        Path_out_final_rvi    = os.path.join(outputfolder_subid,'finalriv.shp')            
        processing.run("native:dissolve", {'INPUT':Path_Temp_final_rviply,'FIELD':['SubId'],'OUTPUT':Path_out_final_rviply},context = context)
        processing.run("native:dissolve", {'INPUT':Path_Temp_final_rvi,'FIELD':['SubId'],'OUTPUT':Path_out_final_rvi},context = context)
        
#        Path_final_rviply = os.path.join(DataFolder,finalrvi_ply_NM)
#        Path_final_riv    = os.path.join(DataFolder,finalriv_NM)
#        Path_Non_ConL_cat = os.path.join(DataFolder,Non_ConnL_Cat_NM)
#        Path_Conl_ply     = os.path.join(DataFolder,ConnL_ply_NM)
#        Path_Non_ConL_ply = os.path.join(DataFolder,Non_ConnL_ply_NM)          
        
                    
        Qgs.exit() 
        
        ######################################################################
        
        



        
     
#        Path_final_rviply = os.path.join(Datafolder,finalrvi_ply_NM)
#        Path_final_riv    = os.path.join(Datafolder,finalriv_NM)
#        Path_Non_ConL_cat = os.path.join(Datafolder,Non_ConnL_Cat_NM)
#        Path_Conl_ply     = os.path.join(Datafolder,ConnL_ply_NM)
#        Path_Non_ConL_ply = os.path.join(Datafolder,Non_ConnL_ply_NM)   
             
###########################################################################3
    def SelectLakes(self,Datafolder,finalrvi_ply_NM = 'finalriv_info_ply.shp',Non_ConnL_Cat_NM = 'Non_con_lake_cat_info.shp',Non_ConnL_ply_NM='Non_Con_Lake_Ply.shp',
                    ConnL_ply_NM='Con_Lake_Ply.shp',finalriv_NM = 'finalriv_info.shp',Thres_Area_Conn_Lakes = -1,Thres_Area_Non_Conn_Lakes = -1,
                    Selection_Method = 'ByArea',sub_colnm = 'SubId',SelectionName = 'All'):

        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing
        from processing.tools import dataobjects
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        context = dataobjects.createContext()
        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)
        
                
        Path_finalcat    = os.path.join(Datafolder,finalrvi_ply_NM)
        finalcat_info    = Path_finalcat[:-3] + "dbf"
        finalcat_info    = Dbf5(finalcat_info)
        finalcat_info    = finalcat_info.to_dataframe()
        
        Path_Non_ConnL_info  = os.path.join(Datafolder,Non_ConnL_Cat_NM)
        Non_ConnL_info       = Path_Non_ConnL_info[:-3] + "dbf"
        Non_ConnL_info       = Dbf5(Non_ConnL_info)
        Non_ConnL_info       = Non_ConnL_info.to_dataframe()
        
        Non_ConnL_info['LakeArea'] = Non_ConnL_info['LakeArea'].astype(float)
        Non_ConnL_info['HyLakeId'] = Non_ConnL_info['HyLakeId'].astype(int)
        finalcat_info['LakeArea']  = finalcat_info['LakeArea'].astype(float)
        finalcat_info['HyLakeId']  = finalcat_info['HyLakeId'].astype(int)  
        finalcat_info['Seg_ID']  = finalcat_info['Seg_ID'].astype(int)        
        
        if Selection_Method == 'ByArea':
            ### process connected lakes first 
            Selected_ConnLakes = finalcat_info[finalcat_info['LakeArea'] >= Thres_Area_Conn_Lakes]['HyLakeId'].values
            Selected_ConnLakes = np.unique(Selected_ConnLakes) 
            Un_Selected_ConnLakes_info  =finalcat_info.loc[(finalcat_info['LakeArea'] < Thres_Area_Conn_Lakes) & (finalcat_info['LakeArea'] > 0)]
            ### process non connected selected lakes 
            if Thres_Area_Non_Conn_Lakes >= 0:
                Selected_Non_ConnLakes = Non_ConnL_info[Non_ConnL_info['LakeArea'] >= Thres_Area_Non_Conn_Lakes]['HyLakeId'].values
                Selected_Non_ConnLakes = np.unique(Selected_Non_ConnLakes)
                Selected_Non_ConnL_info = Non_ConnL_info[Non_ConnL_info['LakeArea'] >= Thres_Area_Non_Conn_Lakes]
            else:
                Selected_Non_ConnLakes = np.full(1,-1) 
                Un_Selected_Non_ConnLakes =  Non_ConnL_info['HyLakeId'].values 
                Selected_Non_ConnL_info = Non_ConnL_info[Non_ConnL_info['LakeArea'] <-1000000]         
        else:
            print(todo)
            
        OutFolderSelectedLakes = os.path.join(Datafolder,SelectionName)    
        if not os.path.exists(OutFolderSelectedLakes):
            os.makedirs(OutFolderSelectedLakes)     
        
        
               
        Selectfeatureattributes(processing,Input = os.path.join(Datafolder,Non_ConnL_Cat_NM) ,Output=os.path.join(OutFolderSelectedLakes,Non_ConnL_Cat_NM),Attri_NM = 'value',Values = Selected_Non_ConnLakes)

        Selectfeatureattributes(processing,Input = os.path.join(Datafolder,Non_ConnL_ply_NM),Output=os.path.join(OutFolderSelectedLakes,Non_ConnL_ply_NM),Attri_NM = 'Hylak_id',Values = Selected_Non_ConnLakes)
        
        Selectfeatureattributes(processing,Input = os.path.join(Datafolder,ConnL_ply_NM),Output=os.path.join(OutFolderSelectedLakes,ConnL_ply_NM),Attri_NM = 'Hylak_id',Values = Selected_ConnLakes)
        
        
        
        #### upate NoncalnUpdateNonConnectedLakeArea_In_Finalcatinfo
        
        Path_Temp_final_rviply = os.path.join(self.tempfolder,'temp_finalriv_ply_selectlake.shp')
        Path_Temp_final_rvi    = os.path.join(self.tempfolder,'temp_finalriv_selectlake.shp')
        
        processing.run("native:dissolve", {'INPUT':os.path.join(Datafolder,finalriv_NM),'FIELD':['SubId'],'OUTPUT':Path_Temp_final_rvi},context = context)
        processing.run("native:dissolve", {'INPUT':os.path.join(Datafolder,finalrvi_ply_NM),'FIELD':['SubId'],'OUTPUT':Path_Temp_final_rviply},context = context)
        
        ####disolve catchment that are covered by non selected connected lakes 
                
        UpdateConnectedLakeArea_In_Finalcatinfo(Path_Temp_final_rviply,Selected_ConnLakes)
        UpdateConnectedLakeArea_In_Finalcatinfo(Path_Temp_final_rvi,Selected_ConnLakes)
        
        
        finalcat_info_temp    = Path_Temp_final_rviply[:-3] + "dbf"
        finalcat_info_temp    = Dbf5(finalcat_info_temp)
        finalcat_info_temp    = finalcat_info_temp.to_dataframe()
        
        mapoldnew_info      = finalcat_info_temp.copy(deep = True)
        mapoldnew_info['nsubid'] = mapoldnew_info['SubId']        
        ### Loop each unselected lake stream seg 
        
        Seg_IDS = Un_Selected_ConnLakes_info['Seg_ID'].values
        Seg_IDS = np.unique(Seg_IDS)
        for iseg in range(0,len(Seg_IDS)):
#            print('#########################################################################################33333')
            i_seg_id   = Seg_IDS[iseg]
            i_seg_info = finalcat_info_temp[finalcat_info_temp['Seg_ID'] == i_seg_id]                      
            i_seg_info = i_seg_info.sort_values(["Seg_order"], ascending = (True))
#            i_seg_info2 = finalcat_info_temp[finalcat_info_temp['Seg_ID'] == i_seg_id] 
#            i_seg_info2 = i_seg_info2.sort_values(["Seg_order"], ascending = (True))
#            print(i_seg_info[['SubId' ,'DowSubId','HyLakeId','Seg_ID','Seg_order']])
#            print(i_seg_info2[['SubId','DowSubId','HyLakeId','Seg_ID','Seg_order']])
            
            ###each part of the segment are not avaiable to be merged 
            N_Hylakeid = np.unique(i_seg_info['HyLakeId'].values)
            if len(i_seg_info) == len(N_Hylakeid):
                continue 
            
            ### All lakes in this segment are removed 
            if np.max(N_Hylakeid)  < 0:   ##
                tsubid        = i_seg_info['SubId'].values[len(i_seg_info) - 1]
                seg_sub_ids   = i_seg_info['SubId'].values
                mapoldnew_info = New_SubId_To_Dissolve(subid = tsubid,catchmentinfo = finalcat_info_temp,mapoldnew_info = mapoldnew_info,ismodifids = 1,modifiidin = seg_sub_ids,mainriv = finalcat_info_temp,Islake = -1,seg_order = 1) 
            
            ### loop from the first order of the current segment
            modifysubids = []
            seg_order = 1
            for iorder in range(0,len(i_seg_info)):
                tsubid = i_seg_info['SubId'].values[iorder]
                modifysubids.append(tsubid)
                
                ### two seg has the same HyLakeId id, can be merged 
                if iorder == len(i_seg_info) - 1:
                    seg_sub_ids   = np.asarray(modifysubids)
                    mapoldnew_info = New_SubId_To_Dissolve(subid = tsubid,catchmentinfo = finalcat_info_temp,mapoldnew_info = mapoldnew_info,ismodifids = 1,modifiidin = seg_sub_ids,mainriv = finalcat_info_temp,Islake = i_seg_info['HyLakeId'].values[iorder],seg_order = seg_order) 
                    modifysubids = []
                    seg_order    = seg_order + 1
                        
                elif i_seg_info['HyLakeId'].values[iorder] == i_seg_info['HyLakeId'].values[iorder + 1]:
                    continue                    
                else: 
                    seg_sub_ids   = np.asarray(modifysubids)
                    mapoldnew_info = New_SubId_To_Dissolve(subid = tsubid,catchmentinfo = finalcat_info_temp,mapoldnew_info = mapoldnew_info,ismodifids = 1,modifiidin = seg_sub_ids,mainriv = finalcat_info_temp,Islake = i_seg_info['HyLakeId'].values[iorder],seg_order = seg_order) 
                    modifysubids = []
                    seg_order    = seg_order + 1
                    
#            i_seg_info2 = mapoldnew_info[mapoldnew_info['Seg_ID'] == i_seg_id] 
#            i_seg_info2 = i_seg_info2.sort_values(["Seg_order"], ascending = (True))
#            print(i_seg_info2[['SubId','DowSubId','HyLakeId','Seg_ID','Seg_order','nsubid']])            

        UpdateTopology(mapoldnew_info,UpdateStreamorder = -1)          
        mapoldnew_info.to_csv( os.path.join(Datafolder,'mapoldnew.csv'),sep=',',index=None)    

        Modify_Feature_info(Path_Temp_final_rviply,mapoldnew_info)
        Modify_Feature_info(Path_Temp_final_rvi,mapoldnew_info)        
        

        UpdateNonConnectedLakeCatchmentinfo(os.path.join(OutFolderSelectedLakes,Non_ConnL_Cat_NM),mapoldnew_info)

        Non_ConL_cat_csv         = os.path.join(OutFolderSelectedLakes,Non_ConnL_Cat_NM)[:-3] + "dbf"
        Non_ConL_cat_info        = Dbf5(Non_ConL_cat_csv)
        Non_ConL_cat_info        = Non_ConL_cat_info.to_dataframe()
        
        UpdateNonConnectedLakeArea_In_Finalcatinfo(Path_Temp_final_rvi,Non_ConL_cat_info)
        UpdateNonConnectedLakeArea_In_Finalcatinfo(Path_Temp_final_rviply,Non_ConL_cat_info)
        
        processing.run("native:dissolve", {'INPUT':Path_Temp_final_rvi,'FIELD':['SubId'],'OUTPUT':os.path.join(OutFolderSelectedLakes,finalriv_NM)},context = context)
        processing.run("native:dissolve", {'INPUT':Path_Temp_final_rviply,'FIELD':['SubId'],'OUTPUT':os.path.join(OutFolderSelectedLakes,finalrvi_ply_NM)},context = context)
                
        return 
        

    def Define_Final_Catchment(self,Datafolder,finalrvi_ply_NM = 'finalriv_info_ply.shp',finalriv_NM = 'finalriv_info.shp',Non_ConnL_Cat_NM = 'Non_con_lake_cat_info.shp',sub_colnm = 'SubId'):
        
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing
        from processing.tools import dataobjects
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        context = dataobjects.createContext()
        context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)
        
        Path_final_rviply = os.path.join(Datafolder,finalrvi_ply_NM)
        Path_final_riv    = os.path.join(Datafolder,finalriv_NM)
        Path_Non_ConL_cat = os.path.join(Datafolder,Non_ConnL_Cat_NM)
         
         
        Path_Temp_final_rviply = os.path.join(self.tempfolder,'temp_finalriv_ply.shp')
        Path_Temp_final_rvi    = os.path.join(self.tempfolder,'temp_finalriv.shp')
        
        processing.run("native:dissolve", {'INPUT':Path_final_rviply,'FIELD':['SubId'],'OUTPUT':Path_Temp_final_rviply},context = context)
        processing.run("native:dissolve", {'INPUT':Path_final_riv,'FIELD':['SubId'],'OUTPUT':Path_Temp_final_rvi},context = context)
        
        ### read riv ply info
        finalrivply_csv     = Path_Temp_final_rviply[:-3] + "dbf"
        finalrivply_info    = Dbf5(finalrivply_csv)
        finalrivply_info    = finalrivply_info.to_dataframe().drop_duplicates(sub_colnm, keep='first')
        mapoldnew_info      = finalrivply_info.copy(deep = True)
        mapoldnew_info['nsubid'] = mapoldnew_info['SubId']
        AllConnectLakeIDS   = finalrivply_info['HyLakeId'].values
        AllConnectLakeIDS   = AllConnectLakeIDS[AllConnectLakeIDS > 0]
        AllConnectLakeIDS   = np.unique(AllConnectLakeIDS)
            
        ### process connected lakes  merge polygons 
        for i in range(0,len(AllConnectLakeIDS)):
            lakeid       = AllConnectLakeIDS[i]
            Lakesub_info = finalrivply_info.loc[finalrivply_info['HyLakeId'] == lakeid]
            Lakesub_info = Lakesub_info.sort_values(["DA"], ascending = (False))
            tsubid       = Lakesub_info[sub_colnm].values[0]
            lakesubids   = Lakesub_info[sub_colnm].values
            mapoldnew_info = New_SubId_To_Dissolve(subid = tsubid,catchmentinfo = finalrivply_info,mapoldnew_info = mapoldnew_info,ismodifids = 1,modifiidin = lakesubids,mainriv = finalrivply_info,Islake = 1) 
        
        UpdateTopology(mapoldnew_info,UpdateStreamorder = -1)          
        mapoldnew_info.to_csv( os.path.join(Datafolder,'mapoldnew.csv'),sep=',',index=None)    

        Modify_Feature_info(Path_Temp_final_rviply,mapoldnew_info)
        Modify_Feature_info(Path_Temp_final_rvi,mapoldnew_info)        

        ## process Non connected lakes
        
        UpdateNonConnectedLakeCatchmentinfo(Path_Non_ConL_cat,mapoldnew_info)

        Non_ConL_cat_csv         = Path_Non_ConL_cat[:-3] + "dbf"
        Non_ConL_cat_info        = Dbf5(Non_ConL_cat_csv)
        Non_ConL_cat_info        = Non_ConL_cat_info.to_dataframe()
        
        UpdateNonConnectedLakeArea_In_Finalcatinfo(Path_final_rviply,Non_ConL_cat_info)
#        Copyfeature_to_another_shp_by_attribute(Path_final_rviply,Path_Non_ConL_cat)

        Path_final_rviply = os.path.join(Datafolder,'finalcat_info.shp')
        Path_final_rvi    = os.path.join(Datafolder,'finalcat_info_riv.shp')
        processing.run("native:dissolve", {'INPUT':Path_Temp_final_rvi,'FIELD':['SubId'],'OUTPUT':Path_final_rvi},context = context)
        processing.run("native:dissolve", {'INPUT':Path_Temp_final_rviply,'FIELD':['SubId'],'OUTPUT':Path_final_rviply},context = context)
        Add_centroid_to_feature(Path_final_rviply,'centroid_x','centroid_y')
        

    def PlotHydrography_Raven(self,Path_rvt_Folder = '#',Path_Hydrographs_output_file=['#'],Scenario_NM = ['#','#']):
        
        
        Obs_rvt_NMS = []
        ###obtain obs rvt file name 
        for file in os.listdir(Path_rvt_Folder):
            if file.endswith(".rvt"):
                Obs_rvt_NMS.append(file)
                
        
        Obs_subids = []                 
        for i in range(0,len(Obs_rvt_NMS)):
            ###find subID
            obs_nm =  Obs_rvt_NMS[i]
            ifilepath = os.path.join(Path_rvt_Folder,obs_nm)
            f = open(ifilepath, "r")
            
            for line in f:
                firstline_info = line.split()
                if firstline_info[0] == ':ObservationData':
                    obssubid  = firstline_info[2]
                    break  ### only read first line
                else:
                    obssubid  = '#' 
                    break     ### only read first line
                    
            Obs_subids.append(obssubid)
            ## this is not a observation rvt file 
            if obssubid =='#':
                continue 
            ####assign column name in the hydrography.csv
            
            colnm_obs = 'sub'+obssubid+' (observed) [m3/s]'
            colnm_sim = 'sub'+obssubid+' [m3/s]'
            colnm_Date = 'date'
            colnm_hr   = 'hour'
            
            ##obtain data from all provided hydrograpy csv output files each hydrograpy csv need has a coorespond scenario name
            Initial_data_frame = 1
            readed_data_correc = 1
            for j in range(0,len(Path_Hydrographs_output_file)):
                
                Path_Hydrographs_output_file_j = Path_Hydrographs_output_file[j]
            
                i_simresult = pd.read_csv(Path_Hydrographs_output_file[0],sep=',')
                colnames = i_simresult.columns
            
                ## check if obs name exist in the hydrograpy csv output files
                if colnm_obs in colnames:
                    
                    ## Initial lize the reaed in data frame 
                    if Initial_data_frame == 1:
                        Readed_Data = i_simresult[[colnm_Date,colnm_hr]]
                        Readed_Data['Obs']  = i_simresult[colnm_obs] 
                        Readed_Data[Scenario_NM[j]] = i_simresult[colnm_sim] 
                        Readed_Data['Date'] = pd.to_datetime(i_simresult[colnm_Date] + ' ' + i_simresult[colnm_hr])
                        Initial_data_frame = -1 
                    else:
                        Readed_Data[Scenario_NM[j]] = i_simresult[colnm_sim].values
                else:
                    readed_data_correc = -1
                    continue 
                    
            if readed_data_correc == -1:
                continue 
                
   
            Readed_Data = Readed_Data.drop(columns=[colnm_Date,colnm_hr])     
            Readed_Data = Readed_Data.set_index('Date') 
            
            Readed_Data = Readed_Data.resample('D').sum()
            Readed_Data['ModelTime'] = Readed_Data.index.strftime('%m-%d-%Y')
            plotGuagelineobs(Scenario_NM,Readed_Data,os.path.join(self.OutputFolder,obs_nm + '.pdf'))
            
            
###########################################################################################33
# Individule functions  not used in 
###################################################################################
    def selectfeaturebasedonID(self,Path_shpfile = '#',sub_colnm = 'SubId',down_colnm = 'DowSubId',mostdownid = [-1],mostupstreamid = [-1],OutputFileName = '#'):
        #selection feature based on subbasin ID
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())

        
        
        hyinfocsv = Path_shpfile[:-3] + "dbf"
        tempinfo = Dbf5(hyinfocsv)
        hyshdinfo2 = tempinfo.to_dataframe().drop_duplicates(sub_colnm, keep='first')
        hyshdinfo =  hyshdinfo2[[sub_colnm,down_colnm]].astype('float').values
        
        
        if mostupstreamid[0] == -1:
            mostupstreamid = np.full(len(mostdownid),-1)
            
        for isub in range(0,len(mostdownid)):
            OutHyID  = mostdownid[isub]
            OutHyID2 = mostupstreamid[isub]
            HydroBasins1 = Defcat(hyshdinfo,OutHyID) ### return fid of polygons that needs to be select 
            if OutHyID2 > 0:
                HydroBasins2 = Defcat(hyshdinfo,OutHyID2)            
    ###  exculde the Ids in HydroBasins2 from HydroBasins1
                for i in range(len(HydroBasins2)):
                    rows =np.argwhere(HydroBasins1 == HydroBasins2[i])
                    HydroBasins1 = np.delete(HydroBasins1, rows)
                HydroBasins = HydroBasins1            
            else:
                HydroBasins = HydroBasins1
        
            exp =sub_colnm + '  IN  (  ' +  str(HydroBasins[0]) #'SubId  IN  ( \'1\',\'1404\',\'1851\') '   
            for i in range(1,len(HydroBasins)):
                exp = exp + " , "+str(HydroBasins[i])  
            exp = exp + ')'
    ### Load feature layers
        
            outfilename = os.path.join(self.OutputFolder,OutBaseName+'_'+str(OutHyID)+'.shp')
            processing.run("native:extractbyexpression", {'INPUT':Path_shpfile,'EXPRESSION':exp,'OUTPUT':outfilename})
        Qgs.exit()
        return 

    def ExtractLakesForGivenWatersheds(self,Path_plyfile = '#',Path_Con_Lake_ply = '#',Path_NonCon_Lake_ply='#',Path_NonClakeinfo = '#'):
        #selection feature based on subbasin ID
        QgsApplication.setPrefixPath(self.qgisPP, True)
        Qgs = QgsApplication([],False)
        Qgs.initQgis()
        from qgis import processing
        from processing.core.Processing import Processing   
        feedback = QgsProcessingFeedback()
        Processing.initialize()
        QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
        
        finalcat_csv     = Path_plyfile[:-3] + "dbf"
        finalcat_info    = Dbf5(finalcat_csv)
        finalcat_info    = finalcat_info.to_dataframe().drop_duplicates('SubId', keep='first')
        
        Connect_Lakeids  = np.unique(finalcat_info['HyLakeId'].values)
        Connect_Lakeids  = Connect_Lakeids[Connect_Lakeids > 0]
        
        SubIds           = np.unique(finalcat_info['SubId'].values)
        SubIds           = SubIds[SubIds > 0]        
        
        NonConn_Lakes    = Path_NonClakeinfo[:-3] + "dbf"
        NonConn_Lakes    = Dbf5(NonConn_Lakes)
        NonConn_Lakes    = NonConn_Lakes.to_dataframe()
        NonConn_Lakes['SubId_riv']    = pd.to_numeric(NonConn_Lakes['SubId_riv'], downcast='float')
        
        NonConn_Lakes_p  = NonConn_Lakes.loc[NonConn_Lakes['SubId_riv'].isin(SubIds)]
        NonCL_Lakeids    = NonConn_Lakes_p['value'].values
        
        
        Selectfeatureattributes(processing,Input = Path_Con_Lake_ply,Output=os.path.join(self.OutputFolder,'Con_Lake_Ply.shp'),Attri_NM = 'Hylak_id',Values = Connect_Lakeids)
        
        Selectfeatureattributes(processing,Input = Path_NonCon_Lake_ply,Output=os.path.join(self.OutputFolder,'Non_Con_Lake_Ply.shp'),Attri_NM = 'Hylak_id',Values = NonCL_Lakeids)
        
        Selectfeatureattributes(processing,Input = Path_NonClakeinfo,Output=os.path.join(self.OutputFolder,'Non_con_lake_cat_info.shp'),Attri_NM = 'value',Values = NonCL_Lakeids)
    
        Qgs.exit()  
        return 



        