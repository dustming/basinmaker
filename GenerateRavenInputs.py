
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
def Generaterivnetwork(hydir,cat,allsubinfo,fac,OutputFoldersub):
    flenriv = copy.copy(hydir)
    flenriv[:,:] = -9999   ##### generate empty river raster
    arcatid = np.unique(cat) #### cat all cat id in target small basin
    arcatid = arcatid[arcatid>=0]
    for i in range(0,len(arcatid)):  #### loop for each catchmant in small basin
        lfid = arcatid[i] ### get the fid in large cat file
        lcatinfo = allsubinfo.loc[allsubinfo['FID'] == lfid] ### Get the cat info in large basin info file
        hyid = lcatinfo['HYBAS_ID'].iloc[0]
        Inhyid = allsubinfo.loc[allsubinfo['NEXT_DOWN'] == hyid]
        if len(Inhyid) > 0:
            for in_i in range(0,len(Inhyid)):
                in_FID = Inhyid['FID'].iloc[in_i]
                pp = np.argwhere(cat == in_FID)
                if len(pp) <= 0:
                    continue
                orow,ocol = Getbasinoutlet(in_FID,cat,fac)
                nrow,ncol = Nextcell(hydir,orow,ocol)
                rowcol = np.full((10000,2),-9999) ### creat two dimension array to store route form beginning to outlet of target catchment
                rowcol [0,0] = nrow
                rowcol [0,1] = ncol
                flen_k = 0
                trow,tcol = Getbasinoutlet(lfid,cat,fac)
                while nrow != trow or ncol != tcol:
                    flen_orow,flen_ocol = nrow,ncol
                    if flen_orow < 0 or flen_ocol<0:
                        break
                    nrow,ncol = Nextcell(hydir,int(flen_orow),int(flen_ocol))
                    flen_k = flen_k + 1
                    rowcol [flen_k,0] = nrow
                    rowcol [flen_k,1] = ncol
                rowcol [flen_k+1,0] = trow
                rowcol [flen_k+1,1] = tcol
                rowcol = rowcol[rowcol[:,0]>=0].astype(int)
                flenriv[rowcol[:,0],rowcol[:,1]] = 1
        else: ### for head watersheds
            if lcatinfo['COAST'].iloc[0] == 1:
                continue
            in_FID = lfid
            trow,tcol = Getbasinoutlet(lfid,cat,fac)
            catrowcol = np.argwhere(cat==in_FID).astype(int)
            catacc = np.full((len(catrowcol),6),-9999)
            catacc[:,0] = catrowcol[:,0]
            catacc[:,1] = catrowcol[:,1]
            catacc[:,2] = fac[catrowcol[:,0],catrowcol[:,1]]
            catacc[:,3] = trow
            catacc[:,4] = tcol
            catacc = catacc[catacc[:,2] > 100]
            if len(catacc) > 0:
                catacc[:,5] = (catacc[:,0] - catacc[:,3])*(catacc[:,0] - catacc[:,3]) + (catacc[:,1] - catacc[:,4])*(catacc[:,1] - catacc[:,4])
                catacc = catacc[catacc[:,5].argsort()]
                nrow,ncol = catacc[len(catacc) - 1,0],catacc[len(catacc) - 1,1]
                rowcol = np.full((10000,2),-9999) ### creat two dimension array to store route form beginning to outlet of target catchment
                rowcol [0,0] = nrow
                rowcol [0,1] = ncol
                flen_k = 0
                while nrow != trow or ncol != tcol:
                    orow,ocol = nrow,ncol
                    if orow < 0 or ocol<0:
                        break
                    nrow,ncol = Nextcell(hydir,orow,ocol)
                    flen_k = flen_k + 1
                    rowcol [flen_k,0] = nrow
                    rowcol [flen_k,1] = ncol
                rowcol [flen_k+1,0] = trow
                rowcol [flen_k+1,1] = tcol
                rowcol = rowcol[rowcol[:,0]>=0].astype(int)
                flenriv[rowcol[:,0],rowcol[:,1]] = 1
    return flenriv

##################################



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
            noout=1
        if not lake[prow-1,pcol-1] == lid:
            noout=1
        if not lake[prow-1,pcol] == lid:
            noout=1
        if not lake[prow,pcol+1] == lid:
            noout=1
        if not lake[prow,pcol-1] == lid:
            noout=1
        if not lake[prow+1,pcol-1] == lid:
            noout=1
        if not lake[prow+1,pcol+1] == lid:
            noout=1
        if not lake[prow+1,pcol] == lid:
            noout=1
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
        if not arclakeid < bsid and arclakeid > blid:
            continue
        arcatid = np.unique(cat[lrowcol[:,0],lrowcol[:,1]]) ###### all catchment id containing this lake
        tarid = 0
        ### if there are more than 1 catchment in cat1, determine if they need to be combined
        ### check if these catchment flow into the lake if it is true, change catchment id into lake catchment id
        if len(arcatid)>1:  #
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
                    if float(len(nlake))/float(len(lrowcol)) > 0.1 and cat[catorow,catocol] > bsid:
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
    return cat
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
            cat[lrowcol[:,0],lrowcol[:,1]] = arclakeid
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
        if Poups[rowcol[0,0],rowcol[0,1]] < 0:
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
        if obs[trow,tcol]  >= 0:
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

def Writervhchanl(ocatinfo,outFolder,lenThres,iscalmanningn):
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
        if (temp >= lenThres):
            catlen = float(temp)/1000 #### in km
            strRlen = str(catlen)
        else:
            catlen = -9999
            strRlen = 'ZERO-'
        if catinfo.iloc[i]['ISLAKE'] >= 0 :
            strRlen = 'ZERO-'
        #####################################################3
        Strcat = str(catid)
        if catid == catinfo.iloc[i]['DOWSUBID']:
            StrDid = str(-1)
        else:
            StrDid = str(int(catinfo.iloc[i]['DOWSUBID']))
        pronam = 'Chn_'+ Strcat
        chslope = max(catinfo.iloc[i]['SLOPE3'],0.0001)
        if chslope < 0:
            chslope = catinfo.iloc[i]['BASINSLOPE']
        writechanel(pronam,max(catinfo.iloc[i]['BKFWIDTH'],1),max(catinfo.iloc[i]['BKFDEPTH'],1),
                    chslope,ochn,catinfo.iloc[i]['MEANELEV'],catinfo.iloc[i]['FLOODP_N'],catinfo.iloc[i]['CH_N'],iscalmanningn)
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
    maxcatid = max(catinfo['SUBID'].values)
    for i in range(0,len(catinfo.index)):
        hruid = int(catinfo.iloc[i]['SUBID'])
        catslope = catinfo.iloc[i]['BASINSLOPE']
        if catinfo.iloc[i]['ISLAKE'] > 0:
            if float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LAKEAREA']):
                catarea2 = float(catinfo.iloc[i]['AREA2'])*max((1-float(catinfo.iloc[i]['LAKERATIO'])),0.05)/1000.00/1000.00
            else:
                catarea2 = float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00 - float(catinfo.iloc[i]['LAKEAREA'])
        else:
            catarea2 = float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00
        StrGid =  str(hruid)+tab
        catid = str(int(catinfo.iloc[i]['SUBID']))+tab
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
        orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
        if catinfo.iloc[i]['ISLAKE'] > 0:
            hruid = int(catinfo.iloc[i]['SUBID']) + int(maxcatid)
            catslope = catinfo.iloc[i]['BASINSLOPE']
            if float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LAKEAREA']):
                catarea2 = float(catinfo.iloc[i]['AREA2'])*min((float(catinfo.iloc[i]['LAKERATIO'])),0.95)/1000/1000
            else:
                catarea2 = float(catinfo.iloc[i]['LAKEAREA'])
            StrGid =  str(hruid)+tab
            catid = str(int(catinfo.iloc[i]['SUBID']))+tab
            StrGidarea = str(catarea2)+tab
            StrGidelev = str(catinfo.iloc[i]['MEANELEV'])+tab
            lat = str(catinfo.iloc[i]['INSIDE_Y'])+tab
            lon = str(catinfo.iloc[i]['INSIDE_X'])+tab
            LAND_USE_CLASS = 'WATER'+tab
            VEG_CLASS = 'WATER'+tab
            SOIL_PROFILE ='SOILPROF'+tab
            AQUIFER_PROFILE ='[NONE]'+tab
            TERRAIN_CLASS ='[NONE]'+tab
            SLOPE = str(catslope)+tab
            ASPECT = '200'+tab
            orvh.write("  "+StrGid+tab+StrGidarea+StrGidelev+lat+lon+catid+LAND_USE_CLASS+VEG_CLASS+SOIL_PROFILE+AQUIFER_PROFILE+TERRAIN_CLASS+SLOPE+ASPECT+"\n")
    orvh.write(":EndHRUs"+"\n")
    orvh.write(":RedirectToFile TestLake.rvh")
    orvh.close()
    ochn.close()
    return catinfo
##############################

#########################################################
def writechanel(chname,chwd,chdep,chslope,orchnl,elev,floodn,channeln,iscalmanningn):
    ### Following SWAT instructions, assume a trapezoidal shape channel, with channel sides has depth and width ratio of 2. zch = 2
    zch = 2
    sidwd = zch * chdep ###river side width
    tab = "          "
    botwd = chwd - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*chwd
        sidwd = 0.5*0.5*chwd
        zch = (chwd - botwd)/2/chdep
    if iscalmanningn >= 0:
        mann = str(channeln)
    else:
        mann = str(0.035)
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
    orchnl.write("    0" + tab + str(floodn) +"\n")
    orchnl.write("    " + str(sidwdfp + 2*chwd)+ tab + mann +"\n")
    orchnl.write("    " + str(sidwdfp + 2*chwd + 2*sidwd + botwd)+ tab + str(floodn) +"\n")
    orchnl.write("  :EndRoughnessZones"+"\n")
    orchnl.write(":EndChannelProfile"+"\n")
    orchnl.write("\n")
    orchnl.write("##############new channel ##############################\n")
#########################################################################################################33

def writelake(catinfo,outFolderraven):
    f2 = open(outFolderraven+"TestLake.rvh","w")
    tab = '       '
    maxcatid = max(catinfo['SUBID'].values)
    for i in range(0,len(catinfo.index)):
        if catinfo.iloc[i]['HYLAKEID'] > 0:
            lakeid = int(catinfo.iloc[i]['HYLAKEID'])
            catid = catinfo.iloc[i]['SUBID']
            if float(catinfo.iloc[i]['AREA2'])/1000.00/1000.00 <= float(catinfo.iloc[i]['LAKEAREA']):
                A = float(catinfo.iloc[i]['AREA2'])*min((float(catinfo.iloc[i]['LAKERATIO'])),0.95)
            else:
                A = float(catinfo.iloc[i]['LAKEAREA'])*1000*1000
            A = catinfo.iloc[i]['LAKEAREA']*1000*1000
            h0 = catinfo.iloc[i]['LAKEDEPTH']
            WeirCoe = 0.6
            hruid = int(catinfo.iloc[i]['SUBID']) + int(maxcatid)
            Crewd = catinfo.iloc[i]['BKFWIDTH']
#            if slakeinfo.iloc[0]['Wshd_area'] < 6000 and slakeinfo.iloc[0]['Wshd_area'] > 0:
        ######write lake information to file
            f2.write(":Reservoir"+ "   Lake_"+ str(int(lakeid))+ "   ######## " +"\n")
            f2.write("  :SubBasinID  "+str(int(catid))+ "\n")
            f2.write("  :HRUID   "+str(int(hruid))+ "\n")
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
def Maphru2force(orank,cat,catinfo,fnrows,fncols,outfolder,forcinggrid,outFolderraven):
    arcpy.RasterToPoint_conversion(outfolder + "finalcat.asc", outfolder + "Finalcat_Point.shp", "VALUE")
    ExtractValuesToPoints(outfolder + "Finalcat_Point.shp", forcinggrid, outfolder + "MapForcing.shp",
                      "NONE", "VALUE_ONLY")
    dbftocsv(outfolder + "MapForcing.dbf",outfolder + "MapForcing.csv")
    Mapforcing = pd.read_csv(outfolder + "MapForcing.csv",sep=",",low_memory=False)
    ogridforc = open(outFolderraven+"GriddedForcings2.txt","w")
    ogridforc.write(":GridWeights" +"\n")
    ogridforc.write("   #      " +"\n")
    ogridforc.write("   # [# HRUs]"+"\n")
#    arcpy.AddMessage(catinfo)
    sNhru = len(catinfo) + len(catinfo[catinfo['ISLAKE'] > 0])
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
        maxcatid = max(catinfo['SUBID'].values)
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
        if catinfo.iloc[i]['ISLAKE'] > 0:
            sumwt = 0.0
            for j in range(0,len(ids)):
                StrGid = str(int(catid) + int(maxcatid))+tab
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
#################################################
def Maphru2forceply(forcingply,outfolder,forcinggrid,outFolderraven,Boundaryply,missrow,misscol):
    arcpy.AddMessage(forcingply[:-3])
    dbf1 = Dbf5(forcingply[:-3]+'dbf')
    Forcinfo = dbf1.to_dataframe()
    Focspre = arcpy.Describe(forcingply).spatialReference
    # run the tool
    if Boundaryply != "#":
        arcpy.Project_management(Boundaryply, outfolder+ "Boundary_freferen.shp",Focspre)
        arcpy.Identity_analysis(outfolder+ "Boundary_freferen.shp", forcingply,outfolder+ "Boundary_Frocing.shp")
        dbf3 = Dbf5(outfolder+ "Boundary_Frocing.dbf")
        BounForc = dbf3.to_dataframe()
        Avafgid = BounForc['FGID'].values
        Avafgid = np.unique(Avafgid)
    else:
        Avafgid = Forcinfo['FGID'].values
    arcpy.Project_management(outfolder+"finalcat_info.shp", outfolder+ "finalcat_freferen.shp",Focspre)
    arcpy.Identity_analysis(outfolder+ "finalcat_freferen.shp", forcingply,outfolder+ "finalcat_Frocing.shp")
    arcpy.Dissolve_management(outfolder+ "finalcat_Frocing.shp", outfolder+ "finalcat_Frocing_diso.shp",["SubID", "FGID"
    ,"Row","Col"], "", "", "")
    if SptailRef.type == "Geographic":
        arcpy.env.CoordinateSystem = arcpy.SpatialReference(3573)####wgs84 - north pore canada
    arcpy.AddField_management(outfolder +"finalcat_Frocing_diso.shp","s_area","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management(outfolder +"finalcat_Frocing_diso.shp","s_area","!shape.area@squaremeters!","PYTHON_9.3","#")
    dbf2 = Dbf5(outfolder+ "finalcat_Frocing_diso.dbf")
    Mapforcing = dbf2.to_dataframe()
    dbf2 = Dbf5(outfolder+ "finalcat_info.dbf")
    Catinfofff = dbf2.to_dataframe()
    catids = Mapforcing['SubId'].values
    catids = np.unique(catids)
    Lakeids = Catinfofff['HyLakeId'].values
#    Lakeids = np.unique(Lakeids)
    Lakeids = Lakeids[Lakeids>0]
    arcpy.AddMessage(str(len(Lakeids)) + "......"+ str(len(catids)))
    ogridforc = open(outFolderraven+"GriddedForcings2.txt","w")
    ogridforc.write(":GridWeights" +"\n")
    ogridforc.write("   #      " +"\n")
    ogridforc.write("   # [# HRUs]"+"\n")
    sNhru = len(catids) + len(Lakeids)
    ogridforc.write("   :NumberHRUs       "+ str(sNhru) + "\n")
    sNcell = (max(Forcinfo['Row'].values)+1+missrow) * (max(Forcinfo['Col'].values)+1+misscol)
    ogridforc.write("   :NumberGridCells  "+str(sNcell)+"\n")
    ogridforc.write("   #            "+"\n")
    ogridforc.write("   # [HRU ID] [Cell #] [w_kl]"+"\n")
#    arcpy.AddMessage(Mapforcing)
    maxcatid = max(catids)
    arcpy.AddMessage(" end of mapping ")
    for i in range(len(catids)):
        catid = catids[i]
        cats = Mapforcing.loc[Mapforcing['SubId'] == catid]
#        arcpy.AddMessage(cats['FGID'])
        cats = cats[cats['FGID'].isin(Avafgid)]
        if len(cats) <= 0:
            cats = Mapforcing.loc[Mapforcing['SubId'] == catid]
            arcpy.AddMessage("Following Grid has to be inluded:.......")
            arcpy.AddMessage(cats['FGID'])
#        arcpy.AddMessage(Mapforcing.loc[Mapforcing['SubId'] == catid])
        tarea = sum(cats['s_area'].values)
        fids = cats['FGID'].values
        fids = np.unique(fids)
        sumwt = 0.0
        for j in range(len(fids)):
            scat = cats[cats['FGID'] == fids[j]]
            if j < len(fids) - 1:
                sarea = sum(scat['s_area'].values)
                wt = float(sarea)/float(tarea)
                sumwt = sumwt + wt
            else:
                wt = 1- sumwt
#            arcpy.AddMessage(scat)
            if(len(scat['Row'].values) > 1):
                arcpy.AddMessage(str(catid)+"error: 1 catchement, 1 grid, produce muti sub")
                Strcellid = str(int(scat['Row'].values[0] * (max(Forcinfo['Col'].values) + 1 +misscol) + scat['Col'].values[0])) + "      "
            else:
                Strcellid = str(int(scat['Row'].values * (max(Forcinfo['Col'].values) + 1 +misscol) + scat['Col'].values)) + "      "                    
                ### str((ncrowcol[0,0] * ncncols + ncrowcol[0,1]))
            ogridforc.write("    "+str(int(catid)) + "     "+Strcellid+str(wt) +"\n")
#        arcpy.AddMessage(cats)
        if Catinfofff.loc[Catinfofff['SubId'] == catid]["IsLake"].values[0] > 0:
#        arcpy.AddMessage(Mapforcing.loc[Mapforcing['SubId'] == catid])
            tarea = sum(cats['s_area'].values)
            fids = cats['FGID'].values
            fids = np.unique(fids)
            sumwt = 0.0
            for j in range(len(fids)):
                scat = cats[cats['FGID'] == fids[j]]
                if j < len(fids) - 1:
                    sarea = sum(scat['s_area'].values)
                    wt = float(sarea)/float(tarea)
                    sumwt = sumwt + wt
                else:
                    wt = 1- sumwt
#            arcpy.AddMessage(scat)
                if(len(scat['Row'].values) > 1):
                    arcpy.AddMessage(str(catid)+"error。。。。。。。。。。。。。。")
                    Strcellid = str(int(scat['Row'].values[0] * (max(Forcinfo['Col'].values)+ 1 + misscol) + scat['Col'].values[0])) + "      "
                else:
                    Strcellid = str(int(scat['Row'].values * (max(Forcinfo['Col'].values)+ 1 + misscol) + scat['Col'].values)) + "      " 
                ### str((ncrowcol[0,0] * ncncols + ncrowcol[0,1]))
                ogridforc.write("    "+str(int(catid) + int(maxcatid)) + "     "+Strcellid+str(wt) +"\n")
    ogridforc.write(":EndGridWeights")
    ogridforc.close()
########
#/* example of calcuate grid index
#       	0	1	2	3	4
#       0	0	1	2	3	4
#       1	5	6	7	8	9
#       2	10	11	12	13	14
#       3	15	16	17	18	19
##  we have 4 rows (0-3) and 5 cols (0-4), the index of each cell
#   should be calaulated by row*(max(colnums)+1) + colnum.
#   for example row =2, col=0, index = 2*(4+1)+0 = 10
#   for example row 3, col 3, index = 3*(4+1)+3 = 18
##################################333
######################################################
###################################################################################33
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
Raveinputsfolder = sys.argv[2]
lenThres = int(sys.argv[3])
iscalmanningn = float(sys.argv[4])

arcpy.env.workspace =OutputFolder
dataset = OutputFolder+"/"+"finalcat_info.shp"

SptailRef = arcpy.Describe(dataset).spatialReference
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) 


if not os.path.exists(Raveinputsfolder):
    os.makedirs(Raveinputsfolder)
dbftocsv(OutputFolder +"/"+ "finalcat_info.dbf",OutputFolder +"/"+ "finalcat_info.csv")
ncatinfo = pd.read_csv(OutputFolder +"/"+"finalcat_info.csv",sep=",",low_memory=False)
ncatinfo2 = ncatinfo.drop_duplicates('SUBID', keep='first')
ncatinfo2.to_csv(OutputFolder +"/"+"finalcatcheck.csv",",")
#ncols = int(arcpy.GetRasterProperties_management(dataset, "COLUMNCOUNT").getOutput(0))
#nrows = int(arcpy.GetRasterProperties_management(dataset, "ROWCOUNT").getOutput(0))
Writervhchanl(ncatinfo2,Raveinputsfolder + "/",lenThres,iscalmanningn)
writelake(ncatinfo2,Raveinputsfolder+ "/")

