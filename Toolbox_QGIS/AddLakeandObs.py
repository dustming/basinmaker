
# coding: utf-8

# In[1]:

########## Write two dimension array into ASCII format
# w_filname  : output file name 
# nraster    : the two dimension array that will be saved in ASCII format
# dataset    : a sample raster dataset. the artribute from this sample dataset 
#              will be used to transfer two dimension array into ASCII 
#return      : None 

def writeraster(w_filname,nraster,dataset):
    orvh = open(w_filname,"w")
    ncols = arcpy.GetRasterProperties_management(dataset, "COLUMNCOUNT")
    nrows = arcpy.GetRasterProperties_management(dataset, "ROWCOUNT")
    xllcorner = arcpy.GetRasterProperties_management(dataset, "LEFT")
    yllcorner = arcpy.GetRasterProperties_management(dataset, "BOTTOM")
    orvh.write("ncols      "+str(ncols) + "\n")
    orvh.write("nrows      "+ str(nrows) + "\n")
    orvh.write("xllcorner    "+str(xllcorner) + "\n")
    orvh.write("yllcorner    "+str(yllcorner) + "\n")
    orvh.write("cellsize     "+str(cellSize) + "\n")
    orvh.write("NODATA_value  -9999" + "\n")
    orvh.close()
    f_handle = open(w_filname, 'a')
    np.savetxt(f_handle,nraster,fmt='%i')
    f_handle.close()


########## Find the downstream grid row and columon number based on flow direction 
# N_dir  : flow direction two dimension array 
# N_row  : The current Row id 
# N_col  : The current column id 
# return N_nrow : row id of the next downstream cell 
# return N_ncol : column id of the next downstream cell  

def Nextcell(N_dir,N_row,N_col):
    
    if N_dir[N_row,N_col] == 1:
        N_nrow = N_row + 0
        N_ncol = N_col + 1
    elif N_dir[N_row,N_col] == 2:
        N_nrow = N_row + 1
        N_ncol = N_col + 1
    elif N_dir[N_row,N_col] == 4:
        N_nrow = N_row + 1
        N_ncol = N_col + 0
    elif N_dir[N_row,N_col] == 8:
        N_nrow = N_row + 1
        N_ncol = N_col - 1
    elif N_dir[N_row,N_col] == 16:
        N_nrow = N_row + 0
        N_ncol = N_col - 1
    elif N_dir[N_row,N_col] == 32:
        N_nrow = N_row - 1
        N_ncol = N_col - 1
    elif N_dir[N_row,N_col] == 64:
        N_nrow = N_row - 1
        N_ncol = N_col + 0
    elif N_dir[N_row,N_col] == 128:
        N_nrow = N_row - 1
        N_ncol = N_col + 1
    else:
        N_nrow = -9999
        N_ncol = -9999
    return N_nrow,N_ncol


########## Return row and column number of a specific catchment's outlet  
# ID     : The specific catchment id 
# basin  : The catchment 2D array  
# fac    : the flow accumulation 2D array  
# dir    : the flow direction 2D array
# nrows, ncols  : maximum number of rows and columns in the basin 2D array   
# return crow,ccol  : row and column id of the catchment outlet.  

def Getbasinoutlet(ID,basin,fac,dir,nrows,ncols):
    catrowcol = np.argwhere(basin==ID).astype(int)
    catacc = np.full((len(catrowcol),3),-9999)
    catacc[:,0] = catrowcol[:,0]
    catacc[:,1] = catrowcol[:,1]
    catacc[:,2] = fac[catrowcol[:,0],catrowcol[:,1]]
    catacc = catacc[catacc[:,2].argsort()]
    ### check if it is a real basin outlet 
    crow = catacc[len(catrowcol)-1,0]
    ccol = catacc[len(catrowcol)-1,1]
          
    nrow,ncol =  Nextcell(dir,crow,ccol)
#    arcpy.AddMessage(str(nrow) + '   '+str(ncol) + '   '+str(crow) + '   '+str(ccol) + '   ')
    if nrow < 0 or ncol < 0:
        return crow, ccol
    elif nrow >= nrows or ncol >= ncols:
        return crow, ccol
    elif basin[nrow,ncol] < 0:
        return crow, ccol
    elif basin[nrow,ncol] != ID:   #  all above means the outlet is the real loutlet 
        return crow, ccol
    else:
        crow = nrow 
        ccol = ncol
        ifound = 0
        for i in range(0,1000): #### find next 1000 grids, to find the basin outlet 
            nrow,ncol =  Nextcell(dir,crow,ccol)
#            arcpy.AddMessage("processing    " + str(nrow) + '   '+str(ncol) + '   '+str(crow) + '   '+str(ccol) + '   ') 
            if nrow < 0 or ncol < 0:
                ifound = 1
                break
            elif nrow >= nrows or ncol >= ncols:
                ifound = 1
                break
            elif basin[nrow,ncol] < 0:
                ifound = 1
                break
            elif basin[nrow,ncol] != ID:
                ifound =  1 #     all above means the outlet is the real loutlet 
                break
            else:
                crow = nrow
                ccol = ncol
                continue
        if ifound == 0: 
            arcpy.AddMessage(" true basin outlet not found for ID...."+ str(ID))
        return crow,ccol        


########## Return row and column number of a specific catchment's outlet  
# hylake       : a dataframe from HydroLAKE database containing all infromation about lakes 
# noncnlake    : a two dimension array of lakes that are not connected by river network 
# NonConLThres : the flow accumulation 2D array  
# dir    : the flow direction 2D array
# nrows, ncols  : maximum number of rows and columns in the basin 2D array   
# return crow,ccol  : row and column id of the catchment outlet.  

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

def Addobspoints(obs,pourpoints,boid,cat):
    obsids = np.unique(obs)
    obsids = obsids[obsids>=0]
    for i in range(0,len(obsids)):
        rowcol = np.argwhere(obs==obsids[i]).astype(int)
        if cat[rowcol[0,0],rowcol[0,1]] > 0:
            if pourpoints[rowcol[0,0],rowcol[0,1]] < 0:
                pourpoints[rowcol[0,0],rowcol[0,1]] = boid + obsids[i]
    return pourpoints

#############################################33

##############################################3
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
        lid = arlakeid[i]
        rowcol = np.argwhere(lake==lid.astype(int))
        nrow = rowcol.shape[0]
        Stridinlake = np.full(nrow,-9999)
        Stridinlake[:] = Str[rowcol[:,0],rowcol[:,1]]
        Strid_L  = np.unique(Stridinlake[np.argwhere(Stridinlake > 0).astype(int)]) ### Get all stream that intercept with lake
        ##### find the intercept point of stream and lake
        for j in range(0,len(Strid_L)):  #### loop for each stream intercept with lake
            strid = Strid_L[j]
            strrowcol = np.argwhere(Str == strid).astype(int)
            nstrrow = strrowcol.shape[0]
            Strchek = np.full((nstrrow,4),-9999)##### 0 row, 1 col, 2 fac,
            Strchek[:,0] = strrowcol[:,0]
            Strchek[:,1] = strrowcol[:,1]
            Strchek[:,2] = fac[strrowcol[:,0],strrowcol[:,1]]
            Strchek = Strchek[Strchek[:,2].argsort()].astype(int)
            ibg=-99
            for irowst in range(nstrrow):#### search from smallest acc stream cell
                if lake[Strchek[irowst,0],Strchek[irowst,1]] == lid and ibg==-99: ### if the begining of stream in lake the stream is ingnored
                    ibg = 1
                    if irowst != 0:
                        if lake[Strchek[irowst-1,0],Strchek[irowst-1,1]] == -9999:
                        ##### this means the stream connect two lakes, so must assign an pourpoints
                            if len(np.unique(lake[Strchek[:,0],Strchek[:,1]])) >= 3:
                                GP_cat[Strchek[irowst-1,0],Strchek[irowst-1,1]] = bsid
                                bsid = bsid + 1
                        ### double check if the head stream cell is nearby the lake, if near by the lake the stream was ignored
                        if Strchek[0,0] != 0 and Strchek[0,0] != nrows -1 and Strchek[0,1] != 0 and Strchek[0,1] != ncols-1:
                            noout = Checklake(Strchek[0,0],Strchek[0,1],nrows,ncols,lid,lake)
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
                                    nostr =orowcol[ka,2]
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
        if lakeacc[len(lakeacc)-1,0] != 0 and lakeacc[len(lakeacc)-1,0] != nrows -1 and lakeacc[len(lakeacc)-1,1] != 0 and lakeacc[len(lakeacc)-1,1] != ncols-1:
            GP_cat[lakeacc[len(lakeacc)-1,0]-1,lakeacc[len(lakeacc)-1,1]-1]=-9999 ### remove the all pourpoints close to lake pourpoints
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]-1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]-1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]+1,lakeacc[len(lakeacc)-1,1]+1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]+1]=-9999
            GP_cat[lakeacc[len(lakeacc)-1,0]-1,lakeacc[len(lakeacc)-1,1]]=-9999
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
                nrow,ncol = Nextcell(hydir,checkcat[len(checkcat)-1,0],checkcat[len(checkcat)-1,1])
                if nrow > 0 or ncol >0:
                    if nrow >= nrows or ncols >= ncols:
                        continue
                    if cat[nrow,ncol] < 0:
                        GP_cat[checkcat[len(checkcat)-1,0],checkcat[len(checkcat)-1,1]] = bcid
                        bcid = bcid + 1
        GP_cat[lakeacc[len(lakeacc)-1,0],lakeacc[len(lakeacc)-1,1]]= sblid
        sblid = sblid + 1
    return GP_cat
###################################################################3

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
######################################################

def CE_mcat4lake(cat1,lake,fac,fdir,bsid,nrows,ncols,Pourpoints):
    #####adjust for some lakes was divided into two catchment beacuse of flow direction and in stream. double check each lake
    ##### and decide if need to merge these catchment into one lake catchment.
    cat = copy.copy(cat1)
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
    outlakeids = np.full(1000000,-99999)
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
            outlakeids[outi] = lakeid
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
    outlakeids= outlakeids[outlakeids > 0]
    return cat,outlakeids
###################################################33
def CE_mcat4lake2(cat1,lake,fac,fdir,bsid,nrows,ncols,Pourpoints):
    cat = copy.copy(cat1)
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>=0]
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
        if len(pp) == 1:
            if arclakeid < 0:
                cat[lrowcol[:,0],lrowcol[:,1]] = pp
            else:
                cat[lrowcol[:,0],lrowcol[:,1]] = arclakeid
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
    return cat
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
    obsids = obsids[obsids>=0]
    for i in range(0,len(obsids)):
        rowcol = np.argwhere(obs==obsids[i]).astype(int)
        if Poups[rowcol[0,0],rowcol[0,1]] < 0 and lake[rowcol[0,0],rowcol[0,1]] < 0:
            Poups[rowcol[0,0],rowcol[0,1]] = ncatid
            ncatid = ncatid + 1
    return Poups
#######
####
# ####################################################33




###################################33




########################################################3
##################################333
#################################################
######################################################
################################################################################33
def Dirpoints(N_dir,p_row,p_col,lake1,lid,goodpoint,k):
    ndir = copy.copy(N_dir)
    ip = copy.copy(k) + 1
#    arcpy.AddMessage(str(p_row) + "    " + str(p_col) + "     " + str(ip))
    if lake1[p_row + 0,p_col + 1] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row + 0,]
        if len(tt[tt[:,1] == p_col + 1,]) < 1:
            ndir[p_row + 0,p_col + 1] = 16
            goodpoint[ip,0] = p_row + 0
            goodpoint[ip,1] = p_col + 1
            ip = ip + 1
    if lake1[p_row + 1,p_col + 1] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row + 1,]
        if len(tt[tt[:,1] == p_col + 1,]) < 1:
            ndir[p_row + 1,p_col + 1] = 32
            goodpoint[ip,0] = p_row + 1
            goodpoint[ip,1] = p_col + 1
            ip = ip + 1
    if lake1[p_row + 1,p_col + 0] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row + 1,]
        if len(tt[tt[:,1] == p_col + 0,]) < 1:
            ndir[p_row + 1,p_col + 0] = 64
            goodpoint[ip,0] = p_row + 1
            goodpoint[ip,1] = p_col + 0
            ip = ip + 1
    if lake1[p_row + 1,p_col - 1] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row + 1,]
        if len(tt[tt[:,1] == p_col - 1,]) < 1:
            ndir[p_row + 1,p_col - 1] = 128
            goodpoint[ip,0] = p_row + 1
            goodpoint[ip,1] = p_col - 1
            ip = ip + 1
    if lake1[p_row + 0,p_col - 1] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row + 0,]
        if len(tt[tt[:,1] == p_col - 1,]) < 1:
            ndir[p_row + 0,p_col - 1] = 1
            goodpoint[ip,0] = p_row + 0
            goodpoint[ip,1] = p_col - 1
            ip = ip + 1
    if lake1[p_row - 1,p_col - 1] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row - 1,]
        if len(tt[tt[:,1] == p_col - 1,]) < 1:
            ndir[p_row - 1,p_col - 1] = 2
            goodpoint[ip,0] = p_row - 1
            goodpoint[ip,1] = p_col - 1
            ip = ip + 1
    if lake1[p_row - 1,p_col + 0] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row - 1,]
        if len(tt[tt[:,1] == p_col + 0,]) < 1:
            ndir[p_row - 1,p_col + 0] = 4
            goodpoint[ip,0] = p_row - 1
            goodpoint[ip,1] = p_col + 0
            ip = ip + 1
    if lake1[p_row - 1,p_col + 1] == lid:
        tt = goodpoint[goodpoint[:,0] == p_row - 1,]
        if len(tt[tt[:,1] == p_col + 1,]) < 1:
            ndir[p_row - 1,p_col + 1] = 8
            goodpoint[ip,0] = p_row - 1
            goodpoint[ip,1] = p_col + 1
            ip = ip + 1
#    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir,goodpoint,ip

def checklakeboundary(lake1,lid,p_row,p_col):
    numnonlake = 0.0
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
    

def Dirpoints2(N_dir,p_row,p_col,lake1,lid,goodpoint,k,ncols,nrows):
    ndir = copy.copy(N_dir)
    ip = copy.copy(k) + 1
#dir 1
    if lake1[p_row + 0,p_col + 1] == lid:  #### it is a lake cell
        if p_row + 0 == nrows-1 or p_col + 1 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row + 1 == nrows-1 or p_col + 1 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row + 1 == nrows-1 or p_col + 0 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row + 1 == nrows-1 or p_col - 1 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row + 0 == nrows-1 or p_col - 1 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row - 1 == nrows-1 or p_col - 1 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row - 1 == nrows-1 or p_col + 0 == ncols - 1:   ### lake not at the boundary of domain
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
        if p_row - 1 == nrows-1 or p_col + 1 == ncols - 1:   ### lake not at the boundary of domain
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

def ChangeDIR(dir,lake1,acc,ncols,nrows,outlakeids,nlakegrids):
    ndir = copy.copy(dir)
    for i in range(0,len(outlakeids)):
        lid = outlakeids[i]
        goodpoint = np.full((20000,2),-99999)
        lrowcol = np.argwhere(lake1==lid).astype(int)
        if len(lrowcol) > nlakegrids:
            continue
#        arcpy.AddMessage(str(lid) + "    " +  str(len(lrowcol)) + "     " + str(i))
        arcpy.AddMessage("start modify lake flow direction, the lake id is    " + str(int(lid)))
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
#                    arcpy.AddMessage("start of checking:    " + str(ip) + "     "+ str(ipo)+ "   ")
#                    arcpy.AddMessage(goodpoint[0:ip,0])
                    trow = goodpoint[i,0]
                    tcol = goodpoint[i,1]
                    if trow >= nrows - 1 or tcol == ncols  - 1:
                        continue
                    ndir,goodpoint,k1= Dirpoints2(ndir,trow,tcol,lake1,lid,goodpoint,k,ncols,nrows)
#                    arcpy.AddMessage("start of checking:    " + str(k1) + "     "+ str(len(goodpoint[goodpoint[:,0]>0,]))+ "   ")
                    k = k1 - 1
            ipo = ip
            ip = len(goodpoint[goodpoint[:,0]>0,]) - 1
            k = ip
#            arcpy.AddMessage("start of checking:    " + str(ip) + "     "+ str(ipo)+ "   " + str(len(lrowcol)))
#        goodpoint[goodpoint[:,0]>0,]
#        ndir[goodpoint[goodpoint[:,0]>0,0].astype(int),goodpoint[goodpoint[:,0]>0,1].astype(int)] = 9
    return ndir
import numpy as np
from scipy.optimize import curve_fit
import copy
import sys
import shutil
import os
import csv
from simpledbf import Dbf5
import pandas as pd
from shutil import copyfile

def AutomatedWatershedsandLakesFilterToolset(OutputFolder = '#',Thre_Lake_Area_Connect = 0,Thre_Lake_Area_nonConnect = -1,
MaximumLakegrids = 10000):
    import os
    import sys
    import tempfile
    import shutil
    import os
    import sys
    import tempfile
    import shutil

    tempinfo = Dbf5(OutputFolder + "/"+'Hylake.dbf')#np.genfromtxt(hyinfocsv,delimiter=',')
    allLakinfo = tempinfo.to_dataframe()
    RoutingToolPath = os.environ['RoutingToolFolder']
    
    
    gisdb =os.path.join(tempfile.gettempdir(), 'grassdata_toolbox')# "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/"
    os.environ['GISDBASE'] = gisdb    
        
    import grass.script as grass
    from grass.script import array as garray
    import grass.script.setup as gsetup
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules import Module
    from grass_session import Session

    os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD'))
    PERMANENT = Session()
    PERMANENT.open(gisdb=gisdb, location='mytest')
    grass.run_command('g.region', raster='dem')
#    grass.run_command('r.mask'  , raster='dem', maskcats = '*',overwrite = True)
    

    OutputFolder = OutputFolder ### the project folder
    VolThreshold =Thre_Lake_Area_Connect ### lake area thresthold for connected lakes 
    NonConLThres = Thre_Lake_Area_nonConnect ### lake area thresthold for non connected lakes 
    nlakegrids = MaximumLakegrids 
    



    temparray = garray.array()
    temparray[:,:] = -9999
    temparray.write(mapname="tempraster", overwrite=True)
##### begin processing 

    blid = 1000000    #### begining of new lake id
    bcid = 1          ## begining of new cat id of hydrosheds
    bsid = 2000000    ## begining of new cat id of in inflow of lakes
    blid2 = 3000000
    boid = 4000000
    
    noncnlake_arr = garray.array(mapname="Nonconnect_Lake")
    conlake_arr = garray.array(mapname="Connect_Lake")
    cat1_arr = garray.array(mapname="cat1")
    str_array = garray.array(mapname="str_grass_r")
    acc_array = garray.array(mapname="acc_grass")
    dir_array = garray.array(mapname="dir_Arcgis")
    obs_array = garray.array(mapname="obs")

    hylake1 = selectlake2(conlake_arr,VolThreshold,allLakinfo) ### remove lakes with lake area smaller than the VolThreshold from connected lake raster 
    if NonConLThres >= 0:
        Lake1 = selectlake(hylake1,noncnlake_arr,NonConLThres,allLakinfo) ### remove lakes with lake area smaller than the NonConLThres from non-connected lake raster 
    else:
        Lake1 = hylake1
    temparray[:,:] = Lake1[:,:]
    temparray.write(mapname="SelectedLakes", overwrite=True)
    grass.run_command('r.null', map='SelectedLakes',setnull=-9999)
    
    ncols = int(temparray.shape[0])
    nrows = int(temparray.shape[1])
    print(ncols,nrows)

    Pourpoints = GenerPourpoint(cat1_arr,Lake1,str_array,nrows,ncols,blid,bsid,bcid,acc_array,dir_array)
#    (cat,Lake1,Str100,nrows,ncols,blid,bsid,bcid,fac,OutputFolder + "/",hydir)
    temparray[:,:] = Pourpoints[:,:]
    temparray.write(mapname="Pourpoints_1", overwrite=True)
    grass.run_command('r.null', map='Pourpoints_1',setnull=-9999)
    grass.run_command('r.to.vect', input='Pourpoints_1',output='Pourpoints_1_F',type='point', overwrite = True)
    grass.run_command('r.stream.basins',direction = 'dir_Grass', points = 'Pourpoints_1_F', basins = 'cat2',overwrite = True)
    
    cat2_array =  garray.array(mapname="cat2")
    temcat,outlakeids =CE_mcat4lake(cat2_array,Lake1,acc_array,dir_array,bsid,nrows,ncols,Pourpoints) 
    temcat2 = CE_Lakeerror(acc_array,dir_array,Lake1,temcat,bsid,blid,boid,nrows,ncols,cat1_arr)
    temparray[:,:] = temcat2[:,:]
    temparray.write(mapname="cat3", overwrite=True)
    grass.run_command('r.null', map='cat3',setnull=-9999)
    
    nPourpoints = GenerateFinalPourpoints(acc_array,dir_array,Lake1,temcat2,bsid,blid,boid,nrows,ncols,cat1_arr,obs_array)
    temparray[:,:] = Pourpoints[:,:]
    temparray.write(mapname="Pourpoints_2", overwrite=True)
    grass.run_command('r.null', map='Pourpoints_2',setnull=-9999)
    grass.run_command('r.to.vect', input='Pourpoints_2',output='Pourpoints_2_F',type='point', overwrite = True)    
    
    ndir = ChangeDIR(dir_array,Lake1,acc_array,ncols,nrows,outlakeids,nlakegrids)
    temparray[:,:] = ndir[:,:]
    temparray.write(mapname="ndir_Arcgis", overwrite=True)
    grass.run_command('r.null', map='ndir_Arcgis',setnull=-9999)    
    grass.run_command('r.reclass', input='ndir_Arcgis',output = 'ndir_Grass',rules = RoutingToolPath + '/'+ 'Arcgis2GrassDIR.txt',overwrite = True)
    
    
    
    grass.run_command('r.stream.basins',direction = 'ndir_Grass', points = 'Pourpoints_2_F', basins = 'cat4',overwrite = True)
    cat4_array =  garray.array(mapname="cat4")
    rowcols = np.argwhere(cat4_array == 0)
    cat4_array[rowcols[:,0],rowcols[:,1]] = -9999
    finalcat = CE_mcat4lake2(cat4_array,Lake1,acc_array,dir_array,bsid,nrows,ncols,nPourpoints)
    temparray[:,:] = finalcat[:,:]
    temparray.write(mapname="finalcat", overwrite=True)
    grass.run_command('r.null', map='finalcat',setnull=-9999)    
    grass.run_command('r.to.vect', input='finalcat',output='finalcat_F',type='area', overwrite = True)    



    PERMANENT.close()

 