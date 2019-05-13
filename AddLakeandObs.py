
# coding: utf-8

# In[1]:


def Defcat(out,outletid):
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((10000000,1),-99999999999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((10000000,1),-999999999999999)
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

##################################################################3
def dbftocsv(filename,outname):
    if filename.endswith('.dbf'):
        print "Converting %s to csv" % filename
        csv_fn = outname
        with open(csv_fn,'wb') as csvfile:
            in_db = dbf.Dbf(filename)
            out_csv = csv.writer(csvfile)
            names = []
            for field in in_db.header.fields:
                names.append(field.name)
            out_csv.writerow(names)
            for rec in in_db:
                out_csv.writerow(rec.fieldData)
            in_db.close()
            print "Done..."
    else:
        print "Filename does not end with .dbf"

##################################################################3
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

##################################################################3
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

##################################################################3
def Getbasinoutlet(ID,basin,fac):
    catrowcol = np.argwhere(basin==ID).astype(int)
    catacc = np.full((len(catrowcol),3),-9999)
    catacc[:,0] = catrowcol[:,0]
    catacc[:,1] = catrowcol[:,1]
    catacc[:,2] = fac[catrowcol[:,0],catrowcol[:,1]]
    catacc = catacc[catacc[:,2].argsort()]
    return catacc[len(catrowcol)-1,0],catacc[len(catrowcol)-1,1]

##################################################################3

##################################################################3
def selectlake(hylake,noncnlake,NonConLThres,hylakeinfo):
#    arcpy.AddMessage("123asdfasdfasdfasd")
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(noncnlake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(noncnlake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['HYLAK_ID'] == sl_lid]
        if len(slakeinfo) <=0:
            continue
        if slakeinfo.iloc[0]['LAKE_AREA'] > NonConLThres:
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
        slakeinfo = hylakeinfo.loc[hylakeinfo['HYLAK_ID'] == sl_lid]
        if len(slakeinfo)<=0:
            arcpy.AddMessage("Lake excluded     asdfasd " + str(sl_lid))
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999
            continue
        if slakeinfo.iloc[0]['LAKE_AREA'] < Lakehres:
            arcpy.AddMessage("Lake excluded     due to area " + str(sl_lid))
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
def GenerPourpoint(cat,lake,Str,nrows,ncols,blid,bsid,bcid,fac,outFolder,hydir):
    GP_cat = copy.copy(cat)
    sblid = copy.copy(blid)
    ############### Part 1 Get all pourpoints of hydroshed catchment
    arcatid = np.unique(cat)#### cat all catchment idd
    arcatid = arcatid[arcatid>=0]
    catoutloc = np.full((len(arcatid),3),-9999)
    for i in range(0,len(arcatid)):
        catid = arcatid[i]
        catrowcol = np.argwhere(cat==catid).astype(int)
        catacc = np.full((len(catrowcol),3),-9999)
        catacc[:,0] = catrowcol[:,0]
        catacc[:,1] = catrowcol[:,1]
        catacc[:,2] = fac[catrowcol[:,0],catrowcol[:,1]]
        catacc = catacc[catacc[:,2].argsort()].astype(int)
        GP_cat[catacc[:,0],catacc[:,1]]=-9999   ### set the catment cells into null
        GP_cat[catacc[len(catrowcol)-1,0],catacc[len(catrowcol)-1,1]]=bcid #### change the outlet of catchment into wid
        bcid = bcid + 1
        catoutloc[i,0] = catid  ## catchment id
        catoutloc[i,1] = catacc[len(catrowcol)-1,0]  #### catchment pourpont row
        catoutloc[i,2] = catacc[len(catrowcol)-1,1]  #### catchment pourpont col
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
def Getcatrivlenslope(catrow,catcol,rivlen,dem,fac,hydir,finalcat,trow,tcol,nrows,ncols,slope):
    rivs = rivlen[catrow,catcol]
    rivs = np.unique(rivs)
    rivs = rivs[rivs > 0]
    rivcinfo = np.full((len(catrow),4),-999999999999)
    rivcinfo[:,0] = rivlen[catrow,catcol]
    rivcinfo[:,1] = fac[catrow,catcol]
    rivcinfo[:,2] = catrow
    rivcinfo[:,3] = catcol
    rivout =  np.full((len(rivs),4),-999999999999.00)
    for i in range(0,len(rivs)):
        rivsid = rivs[i]
        rivcinfo2 = rivcinfo[rivcinfo[:,0]==rivsid,]
        rivcinfo2 = rivcinfo2[rivcinfo2[:,1].argsort()]
        prow = rivcinfo2[0,2].astype(int)
        pcol = rivcinfo2[0,3].astype(int)
        lid = finalcat[prow,pcol]
        nout = Checkcat(prow,pcol,nrows,ncols,lid,finalcat)
        rivtemp = np.full((len(catrow),4),-9999999999.99)
        icell = 0
        if nout > 0:
            nrow = prow
            ncol = pcol
            while finalcat[nrow,ncol] == finalcat[trow,tcol]:
                flen_orow,flen_ocol = nrow,ncol
                if flen_orow < 0 or flen_ocol<0:
                    break
                rivtemp[icell,0] = rivlen[nrow,ncol]
                rivtemp[icell,1] = dem[nrow,ncol]
                rivtemp[icell,3] = slope[nrow,ncol]
                if icell > 0:
                    if rivtemp[icell,0] != rivtemp[icell - 1,0]:
                        rivtemp[icell,2] = rivtemp[icell,0] + rivtemp[icell - 1,2]
                    else:
                        rivtemp[icell,2] = rivtemp[icell-1,2]
                else:
                    rivtemp[icell,2] = rivtemp[icell,0]
                icell = icell + 1
                nrow,ncol = Nextcell(hydir,int(flen_orow),int(flen_ocol))
                if nrow < 0 or ncol < 0:
                    nrow,ncol = Nextcell(hydir,int(trow),int(tcol))
                    print "warning : check river system for catchment: ",finalcat[nrow,ncol],finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol
                    nrow,ncol = 0,0
            rivtemp = rivtemp[rivtemp[:,0]>0,]
            if icell > 0:
                rivout[i,0] = rivtemp[icell-1,2]
                rivout[i,1] = np.mean(rivtemp[:,3])
    rivout = rivout[rivout[:,0]>0,]
    if len(rivout) > 0:
        rivout = rivout[rivout[:,0].argsort()]
        outrivlen = rivout[len(rivout)-1,0]
        outrivslp = rivout[len(rivout)-1,1]
    else:
        outrivlen = -9999.00
        outrivslp = -9999.00
    return outrivlen, outrivslp
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
        lakenrow,lakencol = Nextcell(fdir,lorow,locol)
        if lakenrow < 0 or lakencol < 0:
            lakedowcatid = -999
        else:
            if lakenrow >= nrows or lakencol >= ncols:
                continue
            lakedowcatid = cat[lakenrow,lakencol]
        if not arclakeid < bsid and arclakeid > blid:
            continue
        arcatid,catcounts = np.unique(cat[lrowcol[:,0],lrowcol[:,1]],return_counts=True) ###### all catchment id containing this lake
        tarid = 0
        ### if there are more than 1 catchment in cat1, determine if they need to be combined
        ### check if these catchment flow into the lake if it is true, change catchment id into lake catchment id
        if len(arcatid)>1:  #
            if float(max(catcounts))/float(len(lrowcol)) < 0.9: #and len(lrowcol) < 10000: #: float(max(catcounts))/float(len(lrowcol)) < 0.8 and 
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
        trow,tcol = Getbasinoutlet(GWatids[i],GWat,fac)
        Poups[trow,tcol] = ncatid
        ncatid = ncatid + 1
    OWat = copy.copy(cat)
    OWatids = np.unique(cat)
    OWatids = OWatids[OWatids>=0]
    for i in range(0,len(OWatids)):
        trow,tcol = Getbasinoutlet(OWatids[i],OWat,fac)
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


def Addnlinklakes(fcat,alllake,lake1,fac,sbslid):
    alllakeid = np.unique(alllake)
    sllid = copy.copy(sbslid)
    alllakeid = alllakeid[alllakeid>=0]
    for i in range(0,len(alllakeid)):
        lid = alllakeid[i]
        ibglake = np.argwhere(lake1==lid).astype(int)
        if len(ibglake) == 0: ## this lake is not big lakes
            lrowcol = np.argwhere(alllake==lid).astype(int)
            lcatacc = np.full((len(lrowcol),3),-9999)
            lcatacc[:,0] = lrowcol[:,0]
            lcatacc[:,1] = lrowcol[:,1]
            lcatacc[:,2] = fac[lrowcol[:,0],lrowcol[:,1]]
            lcatacc = lcatacc[lcatacc[:,2].argsort()]
            loutrow = lcatacc[len(lcatacc)-1,0] ### get lake outlet row and col
            loutcol = lcatacc[len(lcatacc)-1,1]
            loutcatids = fcat[lcatacc[:,0],lcatacc[:,1]]
            loutcatids = np.unique(loutcatids)
            if len(loutcatids) == 1:
                for j in range(0,len(lcatacc)):
                    fcat[lcatacc[j,0],lcatacc[j,1]] = sllid + 1
                sllid = sllid + 1
    return fcat


###################################33


def Generatecatinfo(Watseds,fac,fdir,lake,dem,area,hycat,hycatinfo,catinfo,allcatid,lakeinfo,width,depth,
                    rivlen,obs,nrows,ncols,slope):
    finalcat = copy.copy(Watseds)
    for i in range(0,len(allcatid)):
        catid = allcatid[i].astype(int)
        catinfo[i,0] = catid
        rowcol = np.argwhere(finalcat==catid).astype(int)
        trow,tcol = Getbasinoutlet(catid,finalcat,fac)
        nrow,ncol = Nextcell(fdir,trow,tcol)### get the downstream catchment id
        if nrow < 0 or ncol < 0:
            catinfo[i,1] = -1
        elif nrow >= nrows or ncol >= ncols:
            catinfo[i,1] = -1
        elif finalcat[nrow,ncol] < 0:
            catinfo[i,1] = -1
        else:
            catinfo[i,1] = finalcat[nrow,ncol]
        catinfo[i,2] = trow
        catinfo[i,3] = tcol
################################## Get lake information
        lakeid = lake[trow,tcol]
        if lakeid > 0:
            slakeinfo = lakeinfo.loc[lakeinfo['HYLAK_ID'] == lakeid]
            catinfo[i,4] = lakeid
            catinfo[i,5] = slakeinfo.iloc[0]['VOL_TOTAL']
            catinfo[i,6] = slakeinfo.iloc[0]['LAKE_AREA']
            catinfo[i,7] = slakeinfo.iloc[0]['DEPTH_AVG']
            catinfo[i,8] = slakeinfo.iloc[0]['SLOPE_100']
            catinfo[i,9] = slakeinfo.iloc[0]['WSHD_AREA']
            catinfo[i,10] = slakeinfo.iloc[0]['LAKE_TYPE']
########Check if it is observation points
        if obs[trow,tcol]  > 0:
            catinfo[i,23] =  obs[trow,tcol]
########Got basin width and depth
        catwidth,catdepth = Getcatwd(rowcol[:,0],rowcol[:,1],width,depth,-1) ### width depth in m
        catinfo[i,12] = float(sum(dem[rowcol[:,0],rowcol[:,1]])/float(len(rowcol))) ### average elevation
#        catinfo[i,13] = float(sum(area[rowcol[:,0],rowcol[:,1]]))/1000/1000  #### maximum area in km^2
        catinfo[i,14] = max(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
        catinfo[i,15] = min(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
        catinfo[i,16] = dem[trow,tcol] #### outlet elevation
        catinfo[i,17] = max(catwidth,1)
        catinfo[i,18] = max(catdepth,1)
        catinfo[i,19] = 0.030
#######Got basin area and rivlen
        catinfo[i,11] = np.mean(area[rowcol[:,0],rowcol[:,1]])
        catrivlen,catrivslp = Getcatrivlenslope(rowcol[:,0],rowcol[:,1],rivlen,dem,fac,fdir,finalcat,
                                                trow,tcol,nrows,ncols,slope)
        catinfo[i,20] = catrivlen
        catinfo[i,21] = catrivslp
        slopet = slope[rowcol[:,0],rowcol[:,1]]
        slopet = slopet[slopet>0,]
        catinfo[i,22] = np.mean(slopet)
    return catinfo


########################################################3
def Getcatwd(catrow,catcol,width,depth,DA):
    wds = width[catrow,catcol]
    dps = depth[catrow,catcol]
    if max(wds) > 0:
        catwd = max(wds)
        catdps = max(dps)
    else:
        if DA > 0:
            Q = 0.025*DA**0.9302
            catwd = 7.2 *Q **(0.5)
            catdps = 0.27*Q**(0.30)
        else:
            catwd = 15
            catdps = 7.5
    return catwd,catdps
############################################################

def Writervhchanl(ocatinfo,outFolder,nrows,ncols):
    catinfo = copy.copy(ocatinfo)
#    print int(catinfo.iloc[0]['SUBID']),len(catinfo.index)
    ochn = open(outFolder+"modelchannel.rvp","w")
##################3
    orvh = open(outFolder+"test.rvh","w")
    orvh.write("# --------------------------------------------"+"\n")
    orvh.write("# Raven HRU Input file"+"\n")
    orvh.write("#  lake catchment emulation"+"\n")
    orvh.write("# --------------------------------------------"+"\n")
    orvh.write(":SubBasins"+"\n")
    orvh.write("  :Attributes   NAME  DOWNSTREAM_ID       PROFILE REACH_LENGTH  GAUGED"+"\n")
    orvh.write("  :Units        none           none          none           km    none"+"\n")
    tab = "     "
    for i in range(0,len(catinfo.index)):
        ### Get catchment width and dpeth
        catid = int(catinfo.iloc[i]['SUBID'])
        temp = catinfo.iloc[i]['RIVLEN']
        if (temp >= 3000):
            catlen = float(temp)/1000 #### in km
            strRlen = str(catlen)
        else:
            catlen = -9999
            strRlen = 'ZERO-'
        #####################################################3
        Strcat = str(catid)
        StrDid = str(int(catinfo.iloc[i]['DOWSUBID']))
        pronam = 'Chn_'+ Strcat
        chslope = catinfo.iloc[i]['RIVSLOPE']
        if chslope < 0:
            chslope = catinfo.iloc[i]['BASINSLOPE']
        writechanel(pronam,catinfo.iloc[i]['BKFWIDTH'],catinfo.iloc[i]['BKFDEPTH'],
                    chslope,ochn,catinfo.iloc[i]['MEANELEV'])
        if catinfo.iloc[i]['ISOBS'] >= 0 :
            Guage = '1'
        else:
            Guage = '0'
        orvh.write("  "+Strcat+tab+'sub'+Strcat+tab+StrDid+tab+pronam+tab+strRlen+tab+Guage+"\n")
    orvh.write(":EndSubBasins"+"\n")
    orvh.write("\n")
##########################################
    orvh.write(":HRUs"+"\n")
    orvh.write("  :Attributes AREA ELEVATION  LATITUDE  LONGITUDE   BASIN_ID  LAND_USE_CLASS  VEG_CLASS   SOIL_PROFILE  AQUIFER_PROFILE   TERRAIN_CLASS   SLOPE   ASPECT"+"\n")
    orvh.write("  :Units       km2         m       deg        deg       none            none       none           none             none            none     deg      deg"+"\n")
    for i in range(0,len(catinfo.index)):
        catid = int(catinfo.iloc[i]['SUBID'])
        catslope = catinfo.iloc[i]['BASINSLOPE']
        catarea2 = float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00
        StrGid =  str(catid)+tab
        StrGidarea = str(catarea2)+tab
        StrGidelev = str(catinfo.iloc[i]['MEANELEV'])+tab
        lat = str(catinfo.iloc[i]['INSIDE_Y'])+tab
        lon = str(catinfo.iloc[i]['INSIDE_X'])+tab
        LAND_USE_CLASS = 'FOREST'+tab
        VEG_CLASS = 'FOREST'+tab
        SOIL_PROFILE ='SOILPROF'+tab
        AQUIFER_PROFILE ='[NONE]'+tab
        TERRAIN_CLASS ='[NONE]'+tab
        SLOPE = str(catslope)+tab
        ASPECT = '200'+tab
        orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+StrGid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
    orvh.write(":EndHRUs"+"\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return catinfo
##############################

#########################################################
def writechanel(chname,chwd,chdep,chslope,orchnl,elev):
    ### Following SWAT instructions, assume a trapezoidal shape channel, with channel sides has depth and width ratio of 2. zch = 2
    zch = 2
    sidwd = zch * chdep ###river side width
    tab = "          "
    botwd = chwd - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*chwd
        sidwd = 0.5*0.5*chwd
    mann = "0.035"
    zfld = 4 + elev
    zbot = elev - chdep
    sidwdfp = 4/0.25
    Channame = ":ChannelProfile"+tab+chname+tab
    orchnl.write(Channame+"\n")
    Chanslop = "  :Bedslope"+tab+str(chslope)
    orchnl.write(Chanslop+"\n")
    orchnl.write("  :SurveyPoints"+"\n")
    orchnl.write("    0"+tab+str(zfld)+"\n")
    orchnl.write("    "+str(sidwdfp)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd)+tab+str(zbot)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + sidwd + botwd)+tab+str(zbot)+"\n")
    orchnl.write("    "+str(sidwdfp + 2*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(elev)+"\n")
    orchnl.write("    "+str(2*sidwdfp + 4*chwd + 2*sidwd + botwd)+tab+str(zfld)+"\n")
    orchnl.write("  :EndSurveyPoints"+"\n")
    orchnl.write("  :RoughnessZones"+"\n")
    orchnl.write("    0" + tab + mann +"\n")
    orchnl.write("  :EndRoughnessZones"+"\n")
    orchnl.write(":EndChannelProfile"+"\n")
    orchnl.write("\n")
    orchnl.write("##############new channel ##############################\n")
#########################################################################################################33

def writelake(catinfo,outFolderraven):
    f2 = open(outFolderraven+"TestLake.rvh","w")
    tab = '       '
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HYLAKEID'] > 0:
            lakeid = int(catinfo.iloc[i]['HYLAKEID'])
            catid = catinfo.iloc[i]['SUBID']
            A = catinfo.iloc[i]['LAKEAREA']*1000*1000
            h0 = catinfo.iloc[i]['LAKEDEPTH']
            WeirCoe = 0.6
            Crewd = 10
#            if slakeinfo.iloc[0]['Wshd_area'] < 6000 and slakeinfo.iloc[0]['Wshd_area'] > 0:
        ######write lake information to file
            f2.write(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## " +"\n")
            f2.write("  :SubBasinID  "+str(int(catid))+ "\n")
            f2.write("  :HRUID   "+str(int(catid))+ "\n")
            f2.write("  :Type RESROUTE_STANDARD   "+"\n")
            f2.write("  :WeirCoefficient  "+str(WeirCoe)+ "\n")
            f2.write("  :CrestWidth "+str(Crewd)+ "\n")
            f2.write("  :MaxDepth "+str(h0)+ "\n")
            f2.write("  :LakeArea    "+str(A)+ "\n")
            f2.write(":EndReservoir   "+"\n")
            f2.write("#############################################"+"\n")
            f2.write("###New Lake starts"+"\n")
    f2.close()

#################################################################################################################3
def Writecatinfotodbf(OutputFoldersub,catinfo):
    dbfile = OutputFoldersub+ 'finalcat.shp'
    inFeatures = dbfile
    fieldPrecision = 10
    field_scale = 3
    # Execute AddField twice for two new fields
    arcpy.AddField_management(dbfile, "SubId", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "DowSubId", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Area2", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Rivlen", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "RivSlope", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "BasinSlope", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "BkfWidth", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "BkfDepth", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "IsLake", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "HyLakeId", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "LakeVol", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "LakeDepth", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "LakeArea", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Laketype", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "IsObs", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "MeanElev", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    rows = arcpy.UpdateCursor(dbfile)
    for row in rows:
        gridcode = row.gridcode
        sinfo = catinfo[catinfo[:,0] == gridcode,]
        if len(sinfo) > 0:
            row.SubId = sinfo[0,0]
            row.DowSubId = sinfo[0,1]
            row.Area2 = sinfo[0,11]
            row.Rivlen = sinfo[0,20]
            row.RivSlope = sinfo[0,21]
            row.BasinSlope = sinfo[0,22]
            row.BkfWidth = sinfo[0,17]
            row.BkfDepth = sinfo[0,18]
            if(sinfo[0,6] > 0):
                row.IsLake = 1
                row.HyLakeId = sinfo[0,4]
                row.LakeVol = sinfo[0,5]
                row.LakeDepth = sinfo[0,7]
                row.LakeArea = sinfo[0,6]
                row.Laketype = sinfo[0,10]
            else:
                row.IsLake = -9999.99
                row.HyLakeId =  -9999.99
                row.LakeVol =  -9999.99
                row.LakeDepth =  -9999.99
                row.LakeArea =  -9999.99
                row.Laketype = -9999.99
            if sinfo[0,23] > 0:
                row.IsObs = sinfo[0,23]
            else:
                row.IsObs = -9999.99
            row.MeanElev = sinfo[0,12]
        rows.updateRow(row)
    del row
    del rows
##################################333
#################################################
def Maphru2force(orank,cat,catinfo,fnrows,fncols,outfolder,InputsFolder,outFolderraven):
    arcpy.RasterToPoint_conversion(outfolder + "finalcat.asc", outfolder + "Finalcat_Point.shp", "VALUE")
    ExtractValuesToPoints(outfolder + "Finalcat_Point.shp", InputsFolder + "forcingncgrid", outfolder + "MapForcing.shp",
                      "NONE", "VALUE_ONLY")
    dbftocsv(outfolder + "MapForcing.dbf",outfolder + "MapForcing.csv")
    Mapforcing = pd.read_csv(outfolder + "MapForcing.csv",sep=",",low_memory=False)
    ogridforc = open(outFolderraven+"GriddedForcings2.txt","w")
    ogridforc.write(":GridWeights" +"\n")
    ogridforc.write("   #      " +"\n")
    ogridforc.write("   # [# HRUs]"+"\n")
    sNhru = len(catinfo)
    ogridforc.write("   :NumberHRUs       "+ str(sNhru) + "\n")
    sNcell = fnrows*fncols
    ogridforc.write("   :NumberGridCells  "+str(sNcell)+"\n")
    ogridforc.write("   #            "+"\n")
    ogridforc.write("   # [HRU ID] [Cell #] [w_kl]"+"\n")
    ncncols = fncols
    ncnrows = fnrows
    tab = '       '
    for i in range(0,len(catinfo.index)):
        catid = int(catinfo.iloc[i]['SUBID'])
        catmapf = Mapforcing.loc[Mapforcing['GRID_CODE'] == catid]
        rankids = catmapf.values[:,2]
        ids = np.unique(rankids)
        sumwt = 0.0
        for j in range(0,len(ids)):
            StrGid = str(int(catid))+tab
            ncrowcol = np.argwhere(orank==ids[j])
            Strcellid = str((ncrowcol[0,0] * ncncols + ncrowcol[0,1]))+tab
            if len(ids) == 1:
                pesr = 1
            else:
                if j < len(ids) - 1:
                    pesr = float(len(rankids[np.argwhere(rankids == ids[j])]))/float(len(rankids))
                    sumwt = sumwt + pesr
#                print j,pesr,sumwt,float(len(rankids[np.argwhere(rankids == ids[j])])),float(len(rankids))
                else:
                    pesr = 1 - sumwt
#                print j,pesr,sumwt
            ogridforc.write("    "+StrGid+Strcellid+str(pesr) + "\n")
    ogridforc.write(":EndGridWeights")
    ogridforc.close()
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
    

def Dirpoints2(N_dir,p_row,p_col,lake1,lid,goodpoint,k):
    ndir = copy.copy(N_dir)
    ip = copy.copy(k) + 1
#dir 1
    if lake1[p_row + 0,p_col + 1] == lid:  #### it is a lake cell 
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

def ChangeDIR(dir,lake1,acc,ncols,nrows,outlakeids):
    ndir = copy.copy(dir)
    for i in range(0,len(outlakeids)):
        lid = outlakeids[i]
        arcpy.AddMessage("start modify lake flow direction, the lake id is    " + str(int(lid)))
        goodpoint = np.full((20000,2),-99999)
        lrowcol = np.argwhere(lake1==lid).astype(int)
#        if len(lrowcol) > 100:
#            continue
#        arcpy.AddMessage(str(lid) + "    " +  str(len(lrowcol)) + "     " + str(i))
        prow,pcol = Getbasinoutlet(lid,lake1,acc)
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
                    ndir,goodpoint,k1= Dirpoints2(ndir,trow,tcol,lake1,lid,goodpoint,k)
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
import arcpy
from arcpy import env
from arcpy.sa import *
import copy
import sys
import shutil
import os
import csv
from simpledbf import Dbf5
from dbfpy import dbf
import pandas as pd
from shutil import copyfile
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

##### Readed inputs
OutputFolder = sys.argv[1]
VolThreshold = float(sys.argv[2])
NonConLThres = float(sys.argv[3])

arcpy.env.workspace =OutputFolder
os.chdir(OutputFolder)

cellSize = float(arcpy.GetRasterProperties_management(OutputFolder + "/" + "dir", "CELLSIZEX").getOutput(0))
SptailRef = arcpy.Describe(OutputFolder + "/" + "dir").spatialReference
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84

arcpy.env.XYTolerance = cellSize
arcpy.arcpy.env.cellSize = cellSize
arcpy.env.extent = arcpy.Describe( "dir").extent
arcpy.env.snapRaster =  "dir"
###### Read inputs
Null = -9999
blid = 1000000    #### begining of new lake id
bcid = 1       ## begining of new cat id of hydrosheds
bsid = 2000000   ## begining of new cat id of in inflow of lakes
blid2 = 3000000
boid = 4000000
hylake =  arcpy.RasterToNumPyArray(OutputFolder + "/"+'cnlake.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+'cnlake.asc',dtype = 'i4',skiprows = 6) # raster of hydro lake
nchylake =arcpy.RasterToNumPyArray(OutputFolder + "/"+'noncnlake.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+'noncnlake.asc',dtype = 'i4',skiprows = 6) # raster of hydro lake
cat = nchylake =arcpy.RasterToNumPyArray(OutputFolder + "/"+'hybasinfid.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+'hybasinfid.asc',dtype = 'i4',skiprows = 6)   #### raster of hydroshed basin fid
hylakeinfo = pd.read_csv(OutputFolder + "/"+"lakeinfo.csv",sep=",",low_memory=False)       # dataframe of hydrolake database
fac = arcpy.RasterToNumPyArray(OutputFolder + "/"+'acc.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+'acc.asc',dtype = 'i4',skiprows = 6)   # raster of hydrolakes
hydir = arcpy.RasterToNumPyArray(OutputFolder + "/"+'dir.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+'dir.asc',dtype = 'i4',skiprows = 6)   #### raster of hydroshed basin fid
hydem = arcpy.RasterToNumPyArray(OutputFolder + "/"+'dem.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+"dem.asc",dtype = 'i4',skiprows = 6) #### raster of hydroshed dem
obs = arcpy.RasterToNumPyArray(OutputFolder + "/"+'obs.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+"obs.asc",dtype = 'i4',skiprows = 6)
width = arcpy.RasterToNumPyArray(OutputFolder + "/"+'width.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+"width.asc",dtype = 'i4',skiprows = 6)
depth = arcpy.RasterToNumPyArray(OutputFolder + "/"+'depth.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+"depth.asc",dtype = 'i4',skiprows = 6)
allsubinfo = pd.read_csv(OutputFolder + "/"+'hybinfo.csv',sep=",",low_memory=False)
allsubinfo['FID'] = pd.Series(allsubinfo['HYBAS_ID'], index=allsubinfo.index)
allLakinfo = pd.read_csv(OutputFolder + "/"+'lakeinfo.csv',sep=",",low_memory=False)
dataset = "dir"
Str100 = arcpy.RasterToNumPyArray(OutputFolder + "/"+'strlink.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+'strlink.asc',dtype = 'i4',skiprows = 6)
hylake1 = selectlake2(hylake,VolThreshold,allLakinfo)
if NonConLThres > 0:
    Lake1 = selectlake(hylake1,nchylake,NonConLThres,allLakinfo)
else:
    Lake1 = hylake1
writeraster(OutputFolder + "/"+"Lake1.asc",Lake1,dataset)
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"Lake1.prj")
ncols = int(arcpy.GetRasterProperties_management(dataset, "COLUMNCOUNT").getOutput(0))
nrows = int(arcpy.GetRasterProperties_management(dataset, "ROWCOUNT").getOutput(0))
Pourpoints = GenerPourpoint(cat,Lake1,Str100,nrows,ncols,blid,bsid,bcid,fac,OutputFolder + "/",hydir)
writeraster(OutputFolder + "/"+"Pourpoints.asc",Pourpoints,dataset)
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"Pourpoints.prj")
arcpy.ASCIIToRaster_conversion(OutputFolder + "/"+"Pourpoints.asc", OutputFolder + "/"+"Pourpoints1","INTEGER")
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"Pourpoints1.prj")
outWatershed = Watershed("dir", OutputFolder + "/"+"Pourpoints1", "VALUE")
outSetNull = SetNull(outWatershed, outWatershed, "VALUE < 1")
outSetNull.save(OutputFolder + "/"+"temcat1")
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"temcat1.prj")
arcpy.RasterToASCII_conversion("temcat1",  "temcat1.asc")
cat1 =  arcpy.RasterToNumPyArray(OutputFolder + "/"+'temcat1.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+"temcat1.asc",dtype = 'i4',skiprows = 6)
temcat,outlakeids =CE_mcat4lake(cat1,Lake1,fac,hydir,bsid,nrows,ncols,Pourpoints)
temcat2 = CE_Lakeerror(fac,hydir,Lake1,temcat,bsid,blid,boid,nrows,ncols,cat)
writeraster(OutputFolder + "/"+"temcat2.asc",temcat2,dataset)
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"temcat2.prj")
nPourpoints = GenerateFinalPourpoints(fac,hydir,Lake1,temcat2,bsid,blid,boid,nrows,ncols,cat,obs)
writeraster(OutputFolder + "/"+"fPourpoints.asc",nPourpoints,dataset)
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"fPourpoints.prj")
arcpy.ASCIIToRaster_conversion(OutputFolder + "/"+"fPourpoints.asc", OutputFolder + "/"+"fPourpoints1","INTEGER")
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"fPourpoints1.prj")
################ modifiy lake flow directions
#arcpy.AddMessage(len(outlakeids))
ndir = ChangeDIR(hydir,Lake1,fac,ncols,nrows,outlakeids)
writeraster(OutputFolder + "/"+"ndir.asc",ndir,dataset)
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"ndir.prj")
arcpy.ASCIIToRaster_conversion(OutputFolder + "/"+"ndir.asc", OutputFolder + "/"+"ndir","INTEGER")
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"ndir.prj")
##################################3
outWatershed = Watershed(OutputFolder + "/"+"ndir", OutputFolder + "/"+"fPourpoints1", "VALUE")
outSetNull = SetNull(outWatershed, outWatershed, "VALUE < 1")
outExtractByMask = ExtractByMask(outSetNull, OutputFolder+ "/" +"dir")
outExtractByMask.save(OutputFolder + "/"+"temcat3")
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"temcat3.prj")
arcpy.RasterToASCII_conversion(OutputFolder + "/"+'temcat3', OutputFolder + "/"+'temcat3.asc')
temcat3 = arcpy.RasterToNumPyArray(OutputFolder + "/"+'temcat3.asc',nodata_to_value=-9999)#np.loadtxt(OutputFolder + "/"+"temcat3.asc",skiprows = 6)
rowcols = np.argwhere(temcat3 == 0)
temcat3[rowcols[:,0],rowcols[:,1]] = -9999
finalcat = CE_mcat4lake2(temcat3,Lake1,fac,hydir,bsid,nrows,ncols,nPourpoints)
writeraster(OutputFolder + "/"+"finalcat.asc",finalcat,dataset)
copyfile( OutputFolder + "/"+"dir.prj" ,  OutputFolder + "/"+"finalcat.prj")
arcpy.RasterToPolygon_conversion(OutputFolder + "/"+"finalcat.asc",OutputFolder + "/"+"cattemp3.shp", "NO_SIMPLIFY", "VALUE")
copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"cattemp3.prj")
arcpy.Dissolve_management(OutputFolder + "/"+"cattemp3.shp", OutputFolder + "/"+"finalcat.shp", ["gridcode"])
copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"finalcat.prj")
arcpy.AddMessage("The generated catchment with lake is located at OutputFolder with name finalcat.shp")
