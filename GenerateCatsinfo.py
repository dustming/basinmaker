
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


##################################################################3
def selectlake(hylake,Str,hylakeinfo,VolThreshold):
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(sl_lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(sl_lake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['HYLAK_ID'] == sl_lid]
        if slakeinfo.iloc[0]['VOL_TOTAL'] > VolThreshold:
            sl_Strinlake = Str[sl_rowcol[:,0],sl_rowcol[:,1]] ### check the str value within lake
            sl_Strid = np.unique(sl_Strinlake[np.argwhere(sl_Strinlake > 0)]).astype(int) ### Get str id
            if len(sl_Strid) <= 0:  ### if there is no stream or number of lake cells smaller than 5
                sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999   ### change lake value into null
        else:
            sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999
    return sl_lake


##################################################################3
def PreparePrograminputs(N60,VolThreshold,OutHyID,OutputFolder,InputsFolder,Databasefolder):
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(SptailRef.factoryCode) ### WGS84
    ###########Get first step input
    if N60 == 0:
        if HyBasinfile == 'Default':
            hyshdply =  Databasefolder + 'Shapefiles/hybas_na_lev01-12_v1c/' + 'hybas_na_'+cle+'_v1c.shp'
        else:
            hyshdply = InputsFolder + HyBasinfile
        hyshddem =  Databasefolder + 'Rasters/na_dem_15s/na_dem_15s'
        hyshdacc =  Databasefolder + 'Rasters/na_acc_15s/na_acc_15s'
        hyshddir =  Databasefolder + 'Rasters/na_dir_15s/na_dir_15s'
        WidDep = Databasefolder + 'Shapefiles/Width_Depth/narivs.shp'
        hyshdinfo = np.genfromtxt(Databasefolder + 'Shapefiles/hybas_na_lev01-12_v1c/' + 'hybas_na_'+cle+'_v1c.csv',delimiter=',')
    if N60 != 0:
        if HyBasinfile == 'Default':
            hyshdply =  Databasefolder + 'Shapefiles/hybas_ar_lev01-12_v1c/' + 'hybas_ar_'+cle+'_v1c.shp'
        else:
            hyshdply = InputsFolder + HyBasinfile
        hyshddem =  Databasefolder + 'Rasters/n60_dem_15s'
        hyshdacc =  Databasefolder + 'Rasters/n60_acc_15s'
        hyshddir =  Databasefolder + 'Rasters/n60_dir_15s'
        WidDep = Databasefolder + 'Shapefiles/Width_Depth/narivs.shp'
#to do### need to transfer dbf to csv here
        hyshdinfo = np.genfromtxt(Databasefolder + 'Shapefiles/hybas_ar_lev01-12_v1c/' + 'hybas_ar_'+cle+'_v1c.csv',delimiter=',')
    ###############Get interested study regions based on basin outlet hydrobasin ID
    if HyBasinfile == 'Default':
        HydroBasins = Defcat(hyshdinfo,OutHyID)
        out_feature_class = OutputFolder +"HyMask.shp"
        where_clause = '"HYBAS_ID" IN'+ " ("
        for i in range(0,len(HydroBasins)):
            if i == 0:
                where_clause = where_clause + str(HydroBasins[i])
            else:
                where_clause = where_clause + "," + str(HydroBasins[i])
        where_clause = where_clause + ")"
        arcpy.Select_analysis(hyshdply, out_feature_class, where_clause)
    else:
        arcpy.CopyFeatures_management(hyshdply, OutputFolder + "HyMask.shp")
    ### Prepare input files
    ##### dir
    arcpy.env.workspace = OutputFolder
#    arcpy.env.outputCoordinateSystem = "GCS_WGS_1984"
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    outExtractByMask = ExtractByMask(hyshddir, OutputFolder +"HyMask.shp")
    outExtractByMask.save(OutputFolder + "dir")
    arcpy.RasterToASCII_conversion(OutputFolder + "dir", OutputFolder + "dir.asc")
    ####### Set envroment variable to dir
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.extent = arcpy.Describe(OutputFolder + "dir").extent
    arcpy.env.snapRaster = OutputFolder + "dir"
    #########################################################
    ##### dem
    outExtractByMask = ExtractByMask(hyshddem, OutputFolder +"HyMask.shp")
    outExtractByMask.save(OutputFolder + "dem")
    arcpy.RasterToASCII_conversion(OutputFolder + "dem", OutputFolder + "dem.asc")
    ##### acc
    outExtractByMask = ExtractByMask(hyshdacc, OutputFolder +"HyMask.shp")
    outExtractByMask.save(OutputFolder + "acc")
    arcpy.RasterToASCII_conversion(OutputFolder + "acc", OutputFolder + "acc.asc")
    ######################################################################################
    #######hydrobasin
    arcpy.PolygonToRaster_conversion(OutputFolder +"HyMask.shp", "FID", OutputFolder + "hybasinfid",
                                 "CELL_CENTER","NONE", cellSize)
    arcpy.RasterToASCII_conversion(OutputFolder + "hybasinfid", OutputFolder + "hybasinfid.asc")
    #######Prepare Lakes
    if Lakefile != 0:
        arcpy.CopyFeatures_management(Lakefile, OutputFolder +"HyLake1.shp")
#        arcpy.Clip_fanalysis(Lakefile, OutputFolder +"HyMask.shp", OutputFolder +"HyLake1.shp", "")
        where_clause = '"Vol_total"> '+ str(VolThreshold)
        arcpy.Select_analysis(OutputFolder +"HyLake1.shp", OutputFolder +"HyLake.shp", where_clause)
        arcpy.PolygonToRaster_conversion(OutputFolder +"HyLake.shp", "Hylak_id", OutputFolder + "hylake",
                                 "MAXIMUM_COMBINED_AREA","Hylak_id", cellSize)
        arcpy.RasterToASCII_conversion(OutputFolder + "hylake", OutputFolder + "hylake.asc")
    else:
        arcpy.Clip_analysis(Databasefolder + 'Shapefiles/HyLake.shp', OutputFolder +"HyMask.shp", OutputFolder +"HyLake.shp", "")
        tdir = np.loadtxt(OutputFolder+ 'dir.asc',dtype = 'i4',skiprows = 6)
        tlake = copy.copy(tdir)
        tlake[:,:] = -9999
        writeraster(OutputFolder + "hylake.asc",tlake,OutputFolder + 'dir')
    #########prepare obspoints
    arcpy.PointToRaster_conversion(obspoint, "FID",
                                OutputFolder + "obs", "MAXIMUM", "", cellSize)
    arcpy.RasterToASCII_conversion( OutputFolder + "obs", OutputFolder + "obs.asc")
    #######converting dbf to csv files
    dbftocsv( OutputFolder +"HyLake.dbf",OutputFolder +"lakeinfo.csv")
    dbftocsv( OutputFolder +"HyMask.dbf",OutputFolder +"hybinfo.csv")
    ##### width and depth
    arcpy.Clip_analysis(WidDep, OutputFolder +"HyMask.shp", OutputFolder + "WidDep.shp")
    arcpy.PolylineToRaster_conversion(OutputFolder + "WidDep.shp", "WIDTH", OutputFolder + "width",
                                  "MAXIMUM_LENGTH", "NONE", cellSize)
    arcpy.PolylineToRaster_conversion(OutputFolder + "WidDep.shp", "DEPTH", OutputFolder + "depth",
                                  "MAXIMUM_LENGTH", "NONE", cellSize)
    arcpy.RasterToASCII_conversion( OutputFolder + "depth", OutputFolder + "depth.asc")
    arcpy.RasterToASCII_conversion( OutputFolder + "width", OutputFolder + "width.asc")
    if ComplRiv > 0:
        rivacc = np.loadtxt(OutputFolder + "acc.asc",dtype = 'i4',skiprows = 6)
        rivtemp = copy.copy(rivacc)
        rowcol = np.argwhere(rivacc > 100)
        rivtemp[:,:] = -9999
        rivtemp[rowcol[:,0],rowcol[:,1]] = 1
        writeraster(OutputFolder + "riv1.asc",rivtemp,OutputFolder + "dir")
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

def GenrateCatchatt(OutputFoldersub,str100):
    nrowcol = np.argwhere(str100 > 0)
    if len(nrowcol) > 0:
#        arcpy.RasterToPolyline_conversion(OutputFoldersub + "str",OutputFoldersub + "riv1.shp" , "ZERO",0, "NO_SIMPLIFY")
#        copyfile( OutputFoldersub + "/"+"HyMask.prj" ,  OutputFoldersub + "/"+"riv1.prj")
        arcpy.Intersect_analysis([OutputFoldersub + "DrainL1.shp", OutputFoldersub + "finalcat.shp"], OutputFoldersub + "riv1_cat.shp", "ALL", 0, "LINE")
        copyfile( OutputFoldersub + "/"+"HyMask.prj" ,  OutputFoldersub + "/"+"riv1_cat.prj")
    #### Get catchment area length
    if SptailRef.type == "Geographic":
        arcpy.env.CoordinateSystem = arcpy.SpatialReference(3573)####wgs84 - north pore canada

    arcpy.AddField_management(OutputFoldersub +"finalcat.shp","area","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management(OutputFoldersub +"finalcat.shp","area","!shape.area@squaremeters!","PYTHON_9.3","#")

    if len(nrowcol) > 0:
        arcpy.AddField_management(OutputFoldersub +"riv1_cat.shp","length","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
        arcpy.CalculateField_management(OutputFoldersub +"riv1_cat.shp","length","!shape.length@meters!","PYTHON_9.3","#")
    if SptailRef.type == "Geographic":    
        arcpy.ProjectRaster_management(OutputFoldersub+ "dem", OutputFoldersub+ "demproj",arcpy.SpatialReference(3573),"NEAREST")
    else:
        arcpy.CopyRaster_management(OutputFoldersub+ "dem", OutputFoldersub+ "demproj")
    
    arcpy.arcpy.env.cellSize = float(arcpy.GetRasterProperties_management(OutputFoldersub+ "demproj", "CELLSIZEX").getOutput(0))
    arcpy.env.extent = arcpy.Describe(OutputFoldersub + "demproj").extent
    arcpy.env.snapRaster = OutputFoldersub + "demproj"
    arcpy.env.CoordinateSystem = arcpy.SpatialReference(arcpy.Describe(OutputFoldersub + "demproj").spatialReference.factoryCode)
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(arcpy.Describe(OutputFoldersub + "demproj").spatialReference.factoryCode)
    
    outSlope = Slope(OutputFoldersub + "demproj", "DEGREE",1)
#    arcpy.AddMessage("After Slope")
    outDivide = Divide(outSlope, 180.00)
    outTimes = Times(outDivide, 3.1415926) 

    finalslope = Tan(outTimes)
#    arcpy.AddMessage("5       " + str(arcpy.Describe(outSlope).spatialReference.factoryCode) + "      " + str(arcpy.Describe(OutputFoldersub + "demproj").spatialReference.factoryCode))   
    ##############################################################333
    arcpy.env.CoordinateSystem = arcpy.SpatialReference(SptailRef.factoryCode)
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.extent = arcpy.Describe(OutputFoldersub + "dir").extent
    arcpy.env.snapRaster = OutputFoldersub + "dir"
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(SptailRef.factoryCode)
    
    if  SptailRef.type == "Geographic":
        arcpy.ProjectRaster_management(finalslope, OutputFoldersub+ "slope",str(SptailRef.factoryCode),"NEAREST",cellSize)
    else:
        finalslope.save(OutputFoldersub+ "slope")
        
    arcpy.PolygonToRaster_conversion(OutputFoldersub +"finalcat.shp", "area", OutputFoldersub + "area",
                                 "MAXIMUM_COMBINED_AREA","NONE", cellSize)
                                 
    arcpy.RasterToASCII_conversion(OutputFoldersub + "area", OutputFoldersub + "area.asc")
    if len(nrowcol) > 0:
        arcpy.PolylineToRaster_conversion(OutputFoldersub + "riv1_cat.shp", "length", OutputFoldersub + "rivlength",  "MAXIMUM_LENGTH", "NONE", cellSize)
        outExtractByMask = ExtractByMask(OutputFoldersub + "rivlength", OutputFoldersub + "strlnk")
        arcpy.RasterToASCII_conversion(outExtractByMask, OutputFoldersub + "rivlength.asc")
    else:
        writeraster(OutputFoldersub + "rivlength.asc",str100,OutputFoldersub + "dir")
    arcpy.RasterToASCII_conversion(OutputFoldersub + "slope", OutputFoldersub + "slope.asc")

####################################################################3
def Checkcat(prow,pcol,nrows,ncols,lid,lake):
    noout=0
    nearcat = np.full(8,-9)
    ### if the point  (prow,pcol) is at the boundary of catchment
    if prow != 0 and prow != nrows -1 and pcol != 0 and pcol != ncols-1:
        if not lake[prow-1,pcol+1] == lid:
            noout=1
            nearcat[0] = lake[prow-1,pcol+1]
        if not lake[prow-1,pcol-1] == lid:
            noout=1
            nearcat[1] = lake[prow-1,pcol-1]
        if not lake[prow-1,pcol] == lid:
            noout=1
            nearcat[2] = lake[prow-1,pcol]
        if not lake[prow,pcol+1] == lid:
            noout=1
            nearcat[3] = lake[prow,pcol+1]
        if not lake[prow,pcol-1] == lid:
            noout=1
            nearcat[4] = lake[prow,pcol-1]
        if not lake[prow+1,pcol-1] == lid:
            noout=1
            nearcat[5] = lake[prow+1,pcol-1]
        if not lake[prow+1,pcol+1] == lid:
            noout=1
            nearcat[6] = lake[prow+1,pcol+1]
        if not lake[prow+1,pcol] == lid:
            noout=1
            nearcat[7] = lake[prow+1,pcol]
    nearcat = nearcat[nearcat > 0]
    return noout,nearcat

###################################################################3
def Getcatrivlenslope(catrow,catcol,rivlen,dem,fac,hydir,finalcat,trow,tcol,nrows,ncols,slope,rivpath):
    rivs = rivlen[catrow,catcol]
    rivs = np.unique(rivs)
    rivs = rivs[rivs > 0]
    rivcinfo = np.full((len(catrow),4),-999999999999)
    rivcinfo[:,0] = rivlen[catrow,catcol]
    rivcinfo[:,1] = fac[catrow,catcol]
    rivcinfo[:,2] = catrow
    rivcinfo[:,3] = catcol
    rivout =  np.full((len(rivs),4),-999999999999.000000) ### store riv info for each river path started at boundary of catchment 
    rivout2 =  np.full((len(rivs),4),-999999999999.000000) ### store riv info for each river path started at inside of catchment 
    for i in range(0,len(rivs)):
        rivsid = rivs[i]
        rivcinfo2 = rivcinfo[rivcinfo[:,0]==rivsid,]
        rivcinfo2 = rivcinfo2[rivcinfo2[:,1].argsort()]
        prow = rivcinfo2[0,2].astype(int)   ## stream cell with lowerst flow accumulation 
        pcol = rivcinfo2[0,3].astype(int)
        lid = finalcat[prow,pcol]     #### catid of start point of stream cell 
        nout, nearcat = Checkcat(prow,pcol,nrows,ncols,lid,finalcat)  #### check if the point (prow,pcol) is close to the catchment boundary, most upstream cell 
        rivtemp = np.full((len(catrow),4),-9999999999.999999)
        icell = 0
        iscase1 = 0
        if len(nearcat) > 0:  ### check if one of the near cat is the upstream catchment
            for incat in range(0,len(nearcat)):
                inearcat = nearcat[incat]
                incat_trow,incat_tcol = Getbasinoutlet(inearcat,finalcat,fac)
                incat_nrow,incat_ncol = Nextcell(hydir,incat_trow,incat_tcol)### get the downstream catchment id
                if finalcat[incat_nrow,incat_ncol] == lid:
                    iscase1 = 1
        if nout > 0 and iscase1 == 1:  #### this means the strem connect at other catchments. or  this catchment only have one stream some head stream was included
#            arcpy.AddMessage("in riv  case 1   " + str(rivsid) + "     " + str(iscase1))
            nrow = prow
            ncol = pcol
            rivpath[nrow,ncol] = 1
            while finalcat[nrow,ncol] == finalcat[trow,tcol]:
                flen_orow,flen_ocol = nrow,ncol
                if flen_orow < 0 or flen_ocol<0 or icell >= len(rivtemp):
                    break
                rivpath[nrow,ncol] = 1
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
                    print "warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol
                    nrow,ncol = 0,0
                if nrow >= nrows or ncol >= ncols:
#                    arcpy.AddMessage("out of the boundary")
                    break                

            rivtemp = rivtemp[rivtemp[:,0]>0,]
            if icell > 0:
                icell = min(icell,len(rivtemp))
                rivout[i,0] = rivtemp[icell-1,2]
                rivout[i,2] = (max(rivtemp[:,1]) - min(rivtemp[:,1]))/rivtemp[icell-1,2]
                rivtemp = rivtemp[rivtemp[:,3]>=0,]
                if len(rivtemp) > 0:
                    rivout[i,1] = np.mean(rivtemp[:,3])
                else:
                    rivout[i,1] = -9999
        rivtemp2 = np.full((len(catrow),4),-9999999999.999999)
        icell = 0
        if len(rivs) > 0 and iscase1 == 0: # for river  not start inside the catchment 
#            arcpy.AddMessage("in riv  case 2   " + str(rivsid))
            nrow = prow
            ncol = pcol
            rivpath[nrow,ncol] = 1
            while finalcat[nrow,ncol] == finalcat[trow,tcol]:
                flen_orow,flen_ocol = nrow,ncol
                if flen_orow < 0 or flen_ocol<0 or icell >= len(rivtemp2):
                    break
                rivpath[nrow,ncol] = 1
                rivtemp2[icell,0] = rivlen[nrow,ncol]
                rivtemp2[icell,1] = dem[nrow,ncol]
                rivtemp2[icell,3] = slope[nrow,ncol]
                if icell > 0:
                    if rivtemp2[icell,0] != rivtemp2[icell - 1,0]:
                        rivtemp2[icell,2] = rivtemp2[icell,0] + rivtemp2[icell - 1,2]
                    else:
                        rivtemp2[icell,2] = rivtemp2[icell-1,2]
                else:
                    rivtemp2[icell,2] = rivtemp2[icell,0]
                icell = icell + 1
                nrow,ncol = Nextcell(hydir,int(flen_orow),int(flen_ocol))
                if nrow < 0 or ncol < 0:
                    nrow,ncol = Nextcell(hydir,int(trow),int(tcol))
                    print "warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol
                    nrow,ncol = 0,0
                if nrow >= nrows or ncol >= ncols:
#                    arcpy.AddMessage("out of the boundary")
                    break 
            rivtemp2 = rivtemp2[rivtemp2[:,0]>0,]
            if icell > 0:
                icell = min(icell,len(rivtemp2))
                rivout2[i,0] = rivtemp2[icell-1,2]
                rivout2[i,2] = (max(rivtemp2[:,1]) - min(rivtemp2[:,1]))/rivtemp2[icell-1,2]
                rivtemp2 = rivtemp2[rivtemp2[:,3]>=0,]
                if len(rivtemp2) > 0:
                    rivout2[i,1] = np.mean(rivtemp2[:,3])
                else:
                    rivout2[i,1] = -9999
    rivout = rivout[rivout[:,0]>0,]
    rivout2 = rivout2[rivout2[:,0]>0,]
    if len(rivout) > 0:
        rivout = rivout[rivout[:,0].argsort()]
        outrivlen = rivout[len(rivout)-1,0]
        outrivslp = rivout[len(rivout)-1,1]
        outrivslp2 = np.mean(rivout[:,2])
    elif len(rivout2) > 0:
        rivout2 = rivout2[rivout2[:,0].argsort()]
        outrivlen = rivout2[len(rivout2)-1,0]
        outrivslp = rivout2[len(rivout2)-1,1]
        outrivslp2 = np.mean(rivout2[:,2])
    else:
        outrivlen = -9999.00
        outrivslp = -9999.00
        outrivslp2 = -9999.00
    return outrivlen, outrivslp,outrivslp2,rivpath
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
                    rivlen,obs,nrows,ncols,slope,landuse,landuseinfo,Q_Mean,wdlen):
    finalcat = copy.copy(Watseds)
    rivpath = copy.copy(Watseds)
    rivpath[:,:] = -9999
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
            catinfo[i,29] = min(float(len(lake[lake == lakeid]))/float(len(finalcat[finalcat == catid])),1.0)
########Check if it is observation points
        if obs[trow,tcol]  >= 0:
#            arcpy.AddMessage(str(catid)+"      "+str(obs[trow,tcol]))
            catinfo[i,23] =  obs[trow,tcol]
########Got basin width and depth
        catrivlen,catrivslp,catrivslp2,rivpath = Getcatrivlenslope(rowcol[:,0],rowcol[:,1],rivlen,dem,fac,fdir,finalcat,
                                                trow,tcol,nrows,ncols,slope,rivpath)
        catwidth,catdepth,catQ = Getcatwd(catid,finalcat,width,depth,Q_Mean,-1,rivlen,wdlen) ### width depth in m
#        arcpy.AddMessage("catid is    " + str(catid) + "    " + str(catwidth))
        catinfo[i,12] = float(sum(dem[rowcol[:,0],rowcol[:,1]])/float(len(rowcol))) ### average elevation
#        catinfo[i,13] = float(sum(area[rowcol[:,0],rowcol[:,1]]))/1000/1000  #### maximum area in km^2
        catinfo[i,14] = max(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
        catinfo[i,15] = min(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
        catinfo[i,16] = dem[trow,tcol] #### outlet elevation
        catinfo[i,17] = catwidth
        catinfo[i,18] = catdepth
        catinfo[i,19] = 0.035
#######Got basin area and rivlen
        catinfo[i,11] = np.mean(area[rowcol[:,0],rowcol[:,1]])
        catinfo[i,20] = catrivlen
        catinfo[i,21] = catrivslp
        slopet = slope[rowcol[:,0],rowcol[:,1]]
        slopet = slopet[slopet>0,]
        catinfo[i,22] = np.mean(slopet)
        catinfo[i,27] = catrivslp2
        if len(slopet) < 1:
            catinfo[i,22] = 0.001
            catinfo[i,21] = 0.001
            catinfo[i,27] = 0.001
        catinfo[i,25] = Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo)
        catinfo[i,26] = catQ
    writeraster(OutputFolder + "/" + "rivpath.asc",rivpath,OutputFolder + "/" + "dir")
    return catinfo

def Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo):
    catidx = finalcat == catid
    rivids = rivlen > 0
    rivincat = np.logical_and(catidx, rivids)
    Landtypes = landuse[rivincat]
    Landtypeid = np.unique(Landtypes)
    Landtypeid1 = Landtypeid[Landtypeid >= 0]
    Landtypeid2 = Landtypeid1[Landtypeid1 > 0]
    Landtypes = Landtypes[Landtypes > 0]
    if len(Landtypes) > 0 and float(len(Landtypeid2))/float(len(Landtypeid1)) >= 0.1:
        sum = 0.0
        for i in range(0,len(Landtypeid2)):
            iid = Landtypeid2[i]
            sum = sum + landuseinfo[landuseinfo['RasterV'] == iid]['MannV'].values*len(np.argwhere(Landtypes == iid))
        floodn = sum/(len(Landtypes))
    else:
        Landtypes2 = landuse[catidx]
        Landtypes2 = Landtypes2[Landtypes2>0]
        Landtypeiicat = np.unique(Landtypes2)
        Landtypeiicat = Landtypeiicat[Landtypeiicat > 0]
        if len(Landtypes2) > 0:
            sum = 0.0
            for i in range(0,len(Landtypeiicat)):
                iid = Landtypeiicat[i]
                sum = sum + landuseinfo[landuseinfo['RasterV'] == iid]['MannV'].values*len(np.argwhere(Landtypes2 == iid))
            floodn = sum/(len(Landtypes2))
        else:
            floodn = 0.035
#    arcpy.AddMessage(floodn)
    return float(floodn)
########################################################3
def Getcatwd(catid,finalcat,width,depth,Q_Mean,DA,rivpath,wdlen):
    catregs = finalcat == catid
    riverp = rivpath > 0
    rivincat = np.logical_and(catregs, riverp)
    wd = width[rivincat]
    dp = depth[rivincat]
    Q = Q_Mean[rivincat]
    wd = wd[wd > 0]
    dp = dp[dp > 0]
    Q =  Q[Q > 0]
    if len(wd) > 0:
        unique, counts = np.unique(wd, return_counts=True)
        catwd = np.average(unique, weights=counts)
        unique, counts = np.unique(dp, return_counts=True)
        catdps = np.average(unique, weights=counts)
        unique, counts =  np.unique(Q, return_counts=True)
        catQ = np.average(unique, weights=counts)
    else:
        catwd = -9
        catdps = -9
        catQ = -9
    return catwd,catdps,catQ
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
    dbfile = OutputFoldersub+ 'finalcat_info.shp'
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
    arcpy.AddField_management(dbfile, "FloodP_n", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Q_Mean", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Ch_n", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Slope2", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "Slope3", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    arcpy.AddField_management(dbfile, "LakeRatio", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    rows = arcpy.UpdateCursor(dbfile)
    for row in rows:
        gridcode = row.gridcode
        catrow = np.argwhere(catinfo[:,0] == gridcode)
        sinfo = catinfo[catinfo[:,0] == gridcode,]
        if len(sinfo) > 0:
            row.SubId = sinfo[0,0]
            row.DowSubId = sinfo[0,1]
            row.Area2 = sinfo[0,11]
            row.Rivlen = sinfo[0,20]
            isupper = 0
            isdown = 0
            if sinfo[0,17] < 0:
                twidth = sinfo[0,17]
                ccurid = sinfo[0,0]
                while(twidth < 0):
                    downid = catinfo[catinfo[:,0] == ccurid,][0,1]
                    if downid < 0:
                        isdown = 0
                        upcats = catinfo[catinfo[:,1] == ccurid,]
                        upcats = upcats[upcats[:,17]>0,]
                        if len(upcats) <= 0:
                            twidth = 1.2345
                            ccurid = -1
                        else:
                            isupper = 1
                            twidth = np.average(upcats[:,17])
                    else:
                        isdown = 1
                        dowcatinfo = catinfo[catinfo[:,0] == downid,]
                        twidth = dowcatinfo[0,17]
                        ccurid = dowcatinfo[0,0]
#                arcpy.AddMessage(str(sinfo[0,0]) +"       "+ str(twidth))
                if twidth == 1.2345 and ccurid == -1:
                    row.BkfWidth = 1.23456
                    row.BkfDepth = 1.23456
                    row.Q_Mean = 1.2345
                    catinfo[catinfo[:,0] == gridcode,17] = 1.2345
                    catinfo[catinfo[:,0] == gridcode,18] = 1.2345
                    catinfo[catinfo[:,0] == gridcode,26] = 1.2345
                elif isdown == 1:
                    row.BkfWidth = dowcatinfo[0,17]
                    row.BkfDepth = dowcatinfo[0,18]
                    row.Q_Mean = dowcatinfo[0,26]
                    catinfo[catinfo[:,0] == gridcode,17] = dowcatinfo[0,17]
                    catinfo[catinfo[:,0] == gridcode,18] = dowcatinfo[0,18]
                    catinfo[catinfo[:,0] == gridcode,26] = dowcatinfo[0,26]
                elif isupper == 1:
                    row.BkfWidth = np.average(upcats[:,17])
                    row.BkfDepth = np.average(upcats[:,18])
                    row.Q_Mean = np.average(upcats[:,26])
                    catinfo[catinfo[:,0] == gridcode,17] = np.average(upcats[:,17])
                    catinfo[catinfo[:,0] == gridcode,18] = np.average(upcats[:,18])
                    catinfo[catinfo[:,0] == gridcode,26] = np.average(upcats[:,26])
            else:
                row.BkfWidth = sinfo[0,17]
                row.BkfDepth = sinfo[0,18]
                row.Q_Mean = sinfo[0,26]
                catinfo[catinfo[:,0] == gridcode,17] = sinfo[0,17]
                catinfo[catinfo[:,0] == gridcode,18] = sinfo[0,18]
            if(sinfo[0,6] > 0):
                row.IsLake = 1
                row.HyLakeId = sinfo[0,4]
                row.LakeVol = sinfo[0,5]
                row.LakeDepth = sinfo[0,7]
                row.LakeArea = sinfo[0,6]
                row.Laketype = sinfo[0,10]
                row.LakeRatio = sinfo[0,29]
            else:
                row.IsLake = -9999.99
                row.HyLakeId =  -9999.99
                row.LakeVol =  -9999.99
                row.LakeDepth =  -9999.99
                row.LakeArea =  -9999.99
                row.Laketype = -9999.99
                row.LakeRatio = -9999.99
            if sinfo[0,23] >= 0:
                row.IsObs = sinfo[0,23]
            else:
                row.IsObs = -9999.99
            row.FloodP_n = sinfo[0,25]
            row.MeanElev = sinfo[0,12]
            twidth = catinfo[catinfo[:,0] == gridcode,17]
            tdepth = catinfo[catinfo[:,0] == gridcode,18]
            tqmean = catinfo[catinfo[:,0] == gridcode,26]
            gropcat1 = catinfo[catinfo[:,17] == twidth,]
            gropcat2 = gropcat1[gropcat1[:,18] == tdepth,]
            gropcat2 = gropcat2[gropcat2[:,26] == tqmean,]
            gropcat2 = gropcat2[gropcat2[:,27] > 0,]
            if len(gropcat2) <=0:
                tslope = float(sinfo[0,22])
            else:
                tslope = np.average(gropcat2[:,27],weights=gropcat2[:,20])
#            if sinfo[0,20] < 0:
#                arcpy.AddMessage("catchment id " + str(sinfo[0,0]) +"     "+str(twidth)+"   "+str(tdepth))
            if tslope >0:
                n = calculateChannaln(twidth,tdepth,tqmean,tslope)
            else:
                n = calculateChannaln(float(sinfo[0,17]),float(sinfo[0,18]),float(sinfo[0,26]),float(sinfo[0,22]))
#            arcpy.AddMessage(n)
#            arcpy.AddMessage(sinfo[0,21])
            if sinfo[0,21] > 0:
                row.RivSlope = max(sinfo[0,21],0.0)
            elif tslope > 0:
                row.RivSlope  = tslope
            else:
                row.RivSlope  = sinfo[0,22]
            row.BasinSlope = sinfo[0,22]
            row.Slope2 = sinfo[0,27]
            row.Slope3 = tslope
            row.Ch_n = n
            catinfo[catinfo[:,0] == gridcode,28] = n
        rows.updateRow(row)
    del row
    del rows
    return catinfo

def calculateChannaln(width,depth,Q,slope):
    zch = 2
    sidwd = zch * depth ###river side width
    tab = "          "
    botwd = width - 2*sidwd ### river
    if (botwd < 0):
        botwd = 0.5*width
        sidwd = 0.5*0.5*width
        zch = (width - botwd)/2/depth
    Ach = botwd*depth + 2*zch*depth*depth/2
#    arcpy.AddMessage(depth)
#    arcpy.AddMessage(zch)
#    arcpy.AddMessage(botwd)
#    arcpy.AddMessage(width)
#    arcpy.AddMessage(slope)

    Pch = botwd + 2*depth*(1+zch**2)**0.5
    Rch = float(Ach)/float(Pch)  ### in meter
    V = float(Q)/float(Ach)
    n = (Rch**(2.0/3.0))*(slope**(1.0/2.0))/V
    return n
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

cellSize = float(arcpy.GetRasterProperties_management(OutputFolder + "/" + "dir", "CELLSIZEX").getOutput(0))
SptailRef = arcpy.Describe(OutputFolder + "/" + "dir").spatialReference
arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84

arcpy.env.workspace =OutputFolder
os.chdir(OutputFolder)
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
####
hylake =  np.loadtxt(OutputFolder + "/"+'hylake.asc',dtype = 'i4',skiprows = 6) # raster of hydro lake
cat = np.loadtxt(OutputFolder + "/"+'hybasinfid.asc',dtype = 'i4',skiprows = 6)   #### raster of hydroshed basin fid
hylakeinfo = pd.read_csv(OutputFolder + "/"+"lakeinfo.csv",sep=",",low_memory=False)       # dataframe of hydrolake database
fac = np.loadtxt(OutputFolder + "/"+'acc.asc',dtype = 'i4',skiprows = 6)   # raster of hydrolakes
hydir = np.loadtxt(OutputFolder + "/"+'dir.asc',dtype = 'i4',skiprows = 6)   #### raster of hydroshed basin fid
hydem = np.loadtxt(OutputFolder + "/"+"dem.asc",dtype = 'i4',skiprows = 6) #### raster of hydroshed dem
obs = np.loadtxt(OutputFolder + "/"+"obs.asc",dtype = 'i4',skiprows = 6)
width = np.loadtxt(OutputFolder + "/"+"width.asc",skiprows = 6)
depth = np.loadtxt(OutputFolder + "/"+"depth.asc",skiprows = 6)
wdlen = np.loadtxt(OutputFolder + "/"+"WD_Len.asc",skiprows = 6)
Q_mean = np.loadtxt(OutputFolder + "/"+"Q_Mean.asc",skiprows = 6)
landuse = np.loadtxt(OutputFolder + "/"+"landuse.asc",dtype = 'i4',skiprows = 6)
landuseinfo = pd.read_csv(OutputFolder + "/"+'landuseinfo.csv',sep=",",low_memory=False)
allsubinfo = pd.read_csv(OutputFolder + "/"+'hybinfo.csv',sep=",",low_memory=False)
allsubinfo['FID'] = pd.Series(allsubinfo['HYBAS_ID'], index=allsubinfo.index)
allLakinfo = pd.read_csv(OutputFolder + "/"+'lakeinfo.csv',sep=",",low_memory=False)
dataset = "dir"
Lake1 = np.loadtxt(OutputFolder + "/"+'Lake1.asc',dtype = 'i4',skiprows = 6)
Str100 = np.loadtxt(OutputFolder + "/"+'strlink.asc',dtype = 'i4',skiprows = 6)
ncols = int(arcpy.GetRasterProperties_management(dataset, "COLUMNCOUNT").getOutput(0))
nrows = int(arcpy.GetRasterProperties_management(dataset, "ROWCOUNT").getOutput(0))
#######################################3
GenrateCatchatt(OutputFolder + "/",Str100)
rivlen = np.loadtxt(OutputFolder+ "/"+ 'rivlength.asc',skiprows = 6)   #### raster of hydroshed basin fid
area = np.loadtxt(OutputFolder+ "/"+"area.asc",skiprows = 6)
slope = np.loadtxt(OutputFolder+ "/"+"slope.asc",skiprows = 6)
finalcat = np.loadtxt(OutputFolder+ "/"+"finalcat.asc",dtype = 'i4',skiprows = 6)
allcatid = np.unique(finalcat)
allcatid = allcatid[allcatid >= 0]
catinfo = np.full((len(allcatid),40),-99999999999.00000000)
catinfo = Generatecatinfo(finalcat,fac,hydir,Lake1,hydem,area,cat,allsubinfo,catinfo,allcatid,allLakinfo,width,depth,rivlen,obs,nrows,ncols,slope,landuse,landuseinfo,Q_mean,wdlen)
arcpy.CopyFeatures_management(OutputFolder + "/" + "finalcat.shp", OutputFolder + "/" + "finalcat_info")
copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"finalcat_info.prj")
catinfo2 = Writecatinfotodbf(OutputFolder + "/",catinfo)
arcpy.AddGeometryAttributes_management(OutputFolder + "/" + "finalcat_info.shp","CENTROID_INSIDE", "","","")
np.savetxt(OutputFolder+"/"+"catinfo.csv",catinfo2,delimiter=",")
arcpy.AddMessage("The generated catchment with lake and calculated parameters is located at OutputFolder with name finalcat_info.shp")
