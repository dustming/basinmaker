
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
            print(" true basin outlet not found for ID...."+ str(ID))
        return crow,ccol        


##################################################################3

##################################################################3
def selectlake(hylake,Str,hylakeinfo,VolThreshold):
    sl_lake = copy.copy(hylake)
    arlakeid = np.unique(sl_lake)
    arlakeid = arlakeid[arlakeid>=0]
    for i in range(0,len(arlakeid)):
        sl_lid = arlakeid[i] ### get lake id
        sl_rowcol = np.argwhere(sl_lake==sl_lid).astype(int) ### get the row and col of lake
        sl_nrow = sl_rowcol.shape[0]
        slakeinfo = hylakeinfo.loc[hylakeinfo['Hylak_id'] == sl_lid]
        if slakeinfo.iloc[0]['Vol_total'] > VolThreshold:
            sl_Strinlake = Str[sl_rowcol[:,0],sl_rowcol[:,1]] ### check the str value within lake
            sl_Strid = np.unique(sl_Strinlake[np.argwhere(sl_Strinlake > 0)]).astype(int) ### Get str id
            if len(sl_Strid) <= 0:  ### if there is no stream or number of lake cells smaller than 5
                sl_lake[sl_rowcol[:,0],sl_rowcol[:,1]] = -9999   ### change lake value into null
        else:
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
    rivs = np.unique(rivs)  ### get all river length within the catchment 
    rivs = rivs[rivs > 0]
    rivcinfo = np.full((len(catrow),4),-999999999999.99999)
    rivcinfo[:,0] = rivlen[catrow,catcol]   #### store river length 
    rivcinfo[:,1] = fac[catrow,catcol]   ### flow accumulation 
    rivcinfo[:,2] = catrow    #### row 
    rivcinfo[:,3] = catcol    ### col 
    rivout =  np.full((len(rivs),4),-999999999999.000000) ### store riv info for each river path started at boundary of catchment 
    rivout2 =  np.full((len(rivs),4),-999999999999.000000) ### store riv info for each river path started at inside of catchment 
    for i in range(0,len(rivs)):  ## loop for each river segments within the catchment  each river segment was labled with their river legth 
        rivsid = rivs[i] 
        rivcinfo2 = rivcinfo[rivcinfo[:,0]==rivsid,]   ### get information of the river segment i
        rivcinfo2 = rivcinfo2[rivcinfo2[:,1].argsort()]
        prow = rivcinfo2[0,2].astype(int)   ## find the first grids of the river segment 
        pcol = rivcinfo2[0,3].astype(int)   ## find the first grids of the river segment  
        lid = finalcat[prow,pcol]     #### catid of start point of stream cell 
        nout, nearcat = Checkcat(prow,pcol,nrows,ncols,lid,finalcat)  #### check if the point (prow,pcol) is close to the catchment boundary, most upstream cell 
        rivtemp = np.full((len(catrow),4),-9999999999.999999)
        icell = 0
        iscase1 = 0
        if len(nearcat) > 0:  ### check if one of the near cat is the upstream catchment
            for incat in range(0,len(nearcat)):
                inearcat = nearcat[incat]
                incat_trow,incat_tcol = Getbasinoutlet(inearcat,finalcat,fac,hydir,nrows,ncols)
                incat_nrow,incat_ncol = Nextcell(hydir,incat_trow,incat_tcol)### get the downstream catchment id
                if incat_nrow >= nrows or incat_nrow < 0 or incat_ncol >= ncols or incat_ncol < 0:
                    continue
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
                    print("warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol)
                    nrow,ncol = 0,0
                if nrow >= nrows or ncol >= ncols:
#                    arcpy.AddMessage("out of the boundary")
                    break                

            rivtemp = rivtemp[rivtemp[:,0]>0,]
            if icell > 0:
                icell = min(icell,len(rivtemp))
                rivout[i,0] = rivtemp[icell-1,2]
                rivout[i,2] = float(max(abs(rivtemp[0,1] - rivtemp[icell-1,1]),0.5))/float(rivtemp[icell-1,2])
#                arcpy.AddMessage("rivsid     case1  " + str(rivtemp[0,1]) + "       end dem " + str(rivtemp[icell-1,1]) + '     river len    ' + str(rivtemp[icell-1,2]))
#                rivout[i,2] = (max(rivtemp[:,1]) - min(rivtemp[:,1]))/rivtemp[icell-1,2]
                rivtemp = rivtemp[rivtemp[:,3]>=0,]
                if len(rivtemp) > 0:
                    rivout[i,1] = np.mean(rivtemp[:,3])
                else:
                    rivout[i,1] = -9999
        rivtemp2 = np.full((len(catrow),4),-9999999999.999999)
        icell = 0
        if len(rivs) > 0 and iscase1 == 0: # for river segment start inside the catchemnt and do not rieve water from upstream catchment 
#            arcpy.AddMessage("in riv  case 2   " + str(rivsid))
            nrow = prow
            ncol = pcol
            rivpath[nrow,ncol] = 1
            while finalcat[nrow,ncol] == finalcat[trow,tcol]:
                flen_orow,flen_ocol = nrow,ncol
                if flen_orow < 0 or flen_ocol<0 or icell >= len(rivtemp2):
                    break
                rivpath[nrow,ncol] = 1
                rivtemp2[icell,0] = rivlen[nrow,ncol]    #### store riv length of each channel cell 
                rivtemp2[icell,1] = dem[nrow,ncol]       #### store riv dem of each channel cell 
                rivtemp2[icell,3] = slope[nrow,ncol]     #### store slope  of each channel cell 
                if icell > 0: ### start from the second cell 
                    if rivtemp2[icell,0] != rivtemp2[icell - 1,0]:   ### come to a new river segment 
                        rivtemp2[icell,2] = rivtemp2[icell,0] + rivtemp2[icell - 1,2]    ### store cumulated river length of different river segment
                    else:
                        rivtemp2[icell,2] = rivtemp2[icell-1,2]   ### still in old river segment, do not update cumulated river length 
                else:
                    rivtemp2[icell,2] = rivtemp2[icell,0]  ### stroe the river length of the first river segment 
                icell = icell + 1   ## move to next cell 
                nrow,ncol = Nextcell(hydir,int(flen_orow),int(flen_ocol))
                if nrow < 0 or ncol < 0:  ### if the next cell move out of the domain 
                    nrow,ncol = Nextcell(hydir,int(trow),int(tcol))
                    print("warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol)
                    nrow,ncol = 0,0
                if nrow >= nrows or ncol >= ncols: ### if the next cell move out of the domain 
#                    arcpy.AddMessage("out of the boundary")
                    break 
            rivtemp2 = rivtemp2[rivtemp2[:,0]>0,]    ### remove the emmpty positions in the array 
            if icell > 0:   ## means there is river grids  icell should equal with number of river grids 
                icell = min(icell,len(rivtemp2))   
                rivout2[i,0] = rivtemp2[icell-1,2]   ####  get the acuumulated river length for river i
                rivout2[i,2] = float(max(abs(rivtemp2[0,1] - rivtemp2[icell-1,1]),0.5))/float(rivtemp2[icell-1,2])  ### (dem_begin - dem_end)/river length, the difference minimum value is 0.5 m
#                arcpy.AddMessage("rivsid   case 2   " + str(rivtemp2[0,1]) + "       end dem " + str(rivtemp2[icell-1,1]) + '     river len    ' + str(rivtemp2[icell-1,2]))
#                rivout2[i,2] = (max(rivtemp2[:,1]) - min(rivtemp2[:,1]))/rivtemp2[icell-1,2] ### ####  get the river slope  for river i 
                rivtemp2 = rivtemp2[rivtemp2[:,3]>=0,]    
                if len(rivtemp2) > 0:
                    rivout2[i,1] = np.mean(rivtemp2[:,3])  ####  get the river slope  for river average 
                else:
                    rivout2[i,1] = -9999
    rivout = rivout[rivout[:,0]>0,]
    rivout2 = rivout2[rivout2[:,0]>0,]
    if len(rivout) > 0:
        rivout = rivout[rivout[:,0].argsort()]   ### sort with river length of each river segment, the river slopes from longest river segment was used
        outrivlen = rivout[len(rivout)-1,0]
        outrivslp =  rivout[len(rivout)-1,2] ### slope np.mean(slope along the river channel )
        outrivslp2 = rivout[len(rivout)-1,1]   ### slope max (den_b - dem end , 0.5)/ river length 
#        arcpy.AddMessage('final slope ' + str(outrivslp) + "     " +  str(outrivslp2) + "    " + str(outrivlen))
    elif len(rivout2) > 0:
        rivout2 = rivout2[rivout2[:,0].argsort()]
        outrivlen = rivout2[len(rivout2)-1,0]
        outrivslp =  rivout2[len(rivout2)-1,2]
        outrivslp2 = rivout2[len(rivout2)-1,1]
#        arcpy.AddMessage('final slope ' + str(outrivslp) + "     " +  str(outrivslp2))
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
        trow,tcol = Getbasinoutlet(GWatids[i],GWat,fac,fdir,nrows,ncols)
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


def Generatecatinfo(Watseds,fac,fdir,lake,dem,area,catinfo,allcatid,lakeinfo,width,depth,
                    rivlen,obs,nrows,ncols,slope,landuse,landuseinfo,Q_Mean):
    finalcat = copy.copy(Watseds)
    rivpath = copy.copy(Watseds)
    rivpath[:,:] = -9999
    for i in range(0,len(allcatid)):
#        print("subid      " + str(allcatid[i].astype(int)))
        catid = allcatid[i].astype(int)
        catinfo.loc[i,'SubId'] = catid
        rowcol = np.argwhere(finalcat==catid).astype(int)
        trow,tcol = Getbasinoutlet(catid,finalcat,fac,fdir,nrows,ncols)
        nrow,ncol = Nextcell(fdir,trow,tcol)### get the downstream catchment id
        if nrow < 0 or ncol < 0:
            catinfo.loc[i,'DowSubId'] = -1
        elif nrow >= nrows or ncol >= ncols:
            catinfo.loc[i,'DowSubId'] = -1
        elif finalcat[nrow,ncol] <= 0:
            catinfo.loc[i,'DowSubId'] = -1
        else:
            catinfo.loc[i,'DowSubId'] = finalcat[nrow,ncol]
#        catinfo[i,2] = trow
#        catinfo[i,3] = tcol
################################## Get lake information
        lakeid = lake[trow,tcol]
        if lakeid > 0:
            slakeinfo = lakeinfo.loc[lakeinfo['Hylak_id'] == lakeid]
            catinfo.loc[i,'IsLake'] = 1
            catinfo.loc[i,'HyLakeId'] = lakeid
            catinfo.loc[i,'LakeVol'] = slakeinfo.iloc[0]['Vol_total']
            catinfo.loc[i,'LakeArea']= slakeinfo.iloc[0]['Lake_area']
            catinfo.loc[i,'LakeDepth']= slakeinfo.iloc[0]['Depth_avg']
            catinfo.loc[i,'Laketype'] = slakeinfo.iloc[0]['Lake_type']
#            catinfo[i,29] = min(float(len(lake[lake == lakeid]))/float(len(finalcat[finalcat == catid])),1.0)
########Check if it is observation points
        if obs[trow,tcol]  >= 0:
#            arcpy.AddMessage(str(catid)+"      "+str(obs[trow,tcol]))
            catinfo.loc[i,'IsObs'] =  obs[trow,tcol]
########Got basin width and depth
        catrivlen,catrivslp,catrivslp2,rivpath = Getcatrivlenslope(rowcol[:,0],rowcol[:,1],rivlen,dem,fac,fdir,finalcat,
                                                trow,tcol,nrows,ncols,slope,rivpath)
        catwidth,catdepth,catQ = Getcatwd(catid,finalcat,width,depth,Q_Mean,-1,rivlen) ### width depth in m
#        arcpy.AddMessage("catid is    " + str(catid) + "    " + str(catwidth))
        catinfo.loc[i,'MeanElev'] = float(sum(dem[rowcol[:,0],rowcol[:,1]])/float(len(rowcol))) ### average elevation
#        catinfo[i,13] = float(sum(area[rowcol[:,0],rowcol[:,1]]))/1000/1000  #### maximum area in km^2
#        catinfo[i,14] = max(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
#        catinfo[i,15] = min(dem[rowcol[:,0],rowcol[:,1]])  #### maximum dem
#        catinfo[i,16] = dem[trow,tcol] #### outlet elevation
        catinfo.loc[i,'BkfWidth'] = catwidth
        catinfo.loc[i,'BkfDepth'] = catdepth
#        catinfo[i,19] = 0.035
#######Got basin area and rivlen
#        catinfo[i,11] = np.mean(area[rowcol[:,0],rowcol[:,1]])
        catinfo.loc[i,'Rivlen'] = catrivlen
        slopet = slope[rowcol[:,0],rowcol[:,1]]
        slopet = slopet[slopet>0,]
        catinfo.loc[i,'BasinSlope'] = np.mean(slopet)
        if catrivslp < 0:
            catinfo.loc[i,'RivSlope'] = np.mean(slopet)
        else：
            catinfo.loc[i,'RivSlope'] = catrivslp
#        catinfo[i,27] = catrivslp2
        if len(slopet) < 1:
            catinfo.loc[i,'RivSlope'] = 0.001
            catinfo.loc[i,'BasinSlope'] = 0.001
        catinfo.loc[i,'FloodP_n'] = Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo)
        catinfo.loc[i,'Q_Mean'] = catQ
#    writeraster(outputFolder + "/" + "rivpath.asc",rivpath,OutputFolder + "/" + "dir")
    return catinfo,rivpath

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
def Getcatwd(catid,finalcat,width,depth,Q_Mean,DA,rivpath):
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
    arcpy.AddField_management(dbfile, "LakeRatio", "FLOAT", fieldPrecision,field_scale,"", "", "NULLABLE","","")
    rows = arcpy.UpdateCursor(dbfile)
    for row in rows:
        gridcode = row.gridcode
        catrow = np.argwhere(catinfo[:,0] == gridcode)
        sinfo = catinfo[catinfo[:,0] == gridcode,]
        if len(sinfo) > 0:
            row.SubId = sinfo[0,0]
            if sinfo[0,0] == sinfo[0,1]:
                row.DowSubId = -1
            else:
                row.DowSubId = sinfo[0,1]
            row.Area2 = sinfo[0,11]
            row.Rivlen = sinfo[0,20]
            isupper = 0
            isdown = 0
            if sinfo[0,17] < 0:           #### if no bankfulll width data avaiable for this catchment
                twidth = sinfo[0,17]      
                ccurid = sinfo[0,0]       ### ccurid  is the current catchment id 
                while(twidth < 0 and ccurid > 0):        
                    downid = catinfo[catinfo[:,0] == ccurid,][0,1]   ### get the downstream if of current catchment 
                    if downid < 0:                                   ### if no donwstream catchment exist 
                        isdown = 0                                   
                        upcats = catinfo[catinfo[:,1] == ccurid,]         #### get upstream catchment id 
                        upcats = upcats[upcats[:,17]>0,]                  #### check if upstream catchment has bankfull width data
                        if len(upcats) <= 0:                              #### if upstream catchment do not have bankfull width data 
                            twidth = 1.2345                               #### define a default value if not downstream exist and upstream do not have bankfull width data
                            ccurid = -1
                        else:
                            isupper = 1                                   ### if upstream has bankfull width data 
                            twidth = np.average(upcats[:,17])                ###  use the averaged bannkfull width data from upstream 
                    else:                                              ### if down stream exist
                        if ccurid == downid:                           ### if downstream id = current catchment id;; downstream catchemnt is not exist
                            twidth = 1.2345
                            ccurid = -1
                            isdown = -1
                        else:
                            isdown = 1                                    
                            dowcatinfo = catinfo[catinfo[:,0] == downid,]  ### get downstream id  
                            twidth = dowcatinfo[0,17]                      ### update twidth catchment with  the downstream id 
                            ccurid = dowcatinfo[0,0]                       ### set currid with
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
            
            if sinfo[0,21] > 0:
                row.RivSlope = max(sinfo[0,21],0.0)
            else:
                row.RivSlope  = sinfo[0,22]
            
            n = calculateChannaln(row.BkfWidth,row.BkfDepth,row.Q_Mean,row.RivSlope)

            row.BasinSlope = sinfo[0,22]
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
import copy
import sys
import shutil
import os
import csv
from simpledbf import Dbf5
import pandas as pd
from shutil import copyfile

def RoutingNetworkTopologyUpdateToolset(projection = 'default',outputFolder = '#'):
    import os
    import sys
    import tempfile
    import shutil
    import os
    import sys
    import tempfile
    import shutil

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
    
    
    if projection != 'default':

        os.system('gdalwarp ' + "\"" + os.path.join(outputFolder, "dem.tif") +"\""+ "    "+ "\""+os.path.join(outputFolder, "dem_proj.tif")+"\""+ ' -t_srs  ' + "\""+projection+"\"")

        project = Session()
        project.open(gisdb=gisdb, location='mytest_project',create_opts=projection)
        grass.run_command("r.import", input = os.path.join(outputFolder, "dem_proj.tif"), output = 'dem_proj', overwrite = True)
        grass.run_command('g.region', raster='dem_proj')  
        
        grass.run_command('v.proj', location='mytest',mapset = 'PERMANENT', input = 'str_finalcat',overwrite = True)
        grass.run_command('v.proj', location='mytest',mapset = 'PERMANENT', input = 'finalcat_F',overwrite = True)
        
        grass.run_command('v.db.addcolumn', map= 'finalcat_F', columns = "Area double precision") 
        grass.run_command('v.db.addcolumn', map= 'str_finalcat', columns = "Length double precision")
              
        
        grass.run_command('v.to.db', map= 'finalcat_F',option = 'area',columns = "Area", units = 'meters') 
        grass.run_command('v.to.db', map= 'str_finalcat',option = 'length', columns = "Length",units = 'meters')
        
        grass.run_command('r.slope.aspect', elevation= 'dem_proj',slope = 'slope',aspect = 'aspect',precision = 'DCELL',overwrite = True)
        grass.run_command('r.mapcalc',expression = 'tanslopedegree = tan(slope) ',overwrite = True) 
        project.close 
        
        
        PERMANENT = Session()
        PERMANENT.open(gisdb=gisdb, location='mytest')
        grass.run_command('g.region', raster='dem')  
        
        grass.run_command('v.proj', location='mytest_project',mapset = 'PERMANENT', input = 'str_finalcat',overwrite = True)
        grass.run_command('v.proj', location='mytest_project',mapset = 'PERMANENT', input = 'finalcat_F',overwrite = True) 
        grass.run_command('r.proj', location='mytest_project',mapset = 'PERMANENT', input = 'tanslopedegree',overwrite = True) 
        grass.run_command('r.proj', location='mytest_project',mapset = 'PERMANENT', input = 'aspect',overwrite = True) 
             
    else:
        PERMANENT = Session()
        PERMANENT.open(gisdb=gisdb, location='mytest')
        grass.run_command('g.region', raster='dem')
        
        grass.run_command('v.db.addcolumn', map= 'finalcat_F', columns = "Area double precision") 
        grass.run_command('v.db.addcolumn', map= 'str_finalcat', columns = "Length double precision")
              
        
        grass.run_command('v.to.db', map= 'finalcat_F',option = 'area',columns = "Area", units = 'meters') 
        grass.run_command('v.to.db', map= 'str_finalcat',option = 'length', columns = "Length",units = 'meters')
        
        grass.run_command('r.slope.aspect', elevation= 'dem_proj',slope = 'slope',aspect = 'aspect',precision = 'DCELL',overwrite = True)
        grass.run_command('r.mapcalc',expression = 'tanslopedegree = tan(slope) ',overwrite = True)         

        
    grass.run_command('v.to.rast',input = 'finalcat_F',output = 'Area',use = 'attr',attribute_column = 'Area',overwrite = True)  
    grass.run_command('v.to.rast',input = 'str_finalcat',output = 'Length',use = 'attr',attribute_column = 'Length',overwrite = True)
    
    grass.run_command('r.out.gdal', input = 'Length',output = outputFolder  + 'rivlength.tif',format= 'GTiff',overwrite = True)
    grass.run_command('r.out.gdal', input = 'finalcat',output = outputFolder  + 'finalcat.tif',format= 'GTiff',overwrite = True)
    
    grass.run_command('r.out.gdal', input = 'tanslopedegree',output = outputFolder  + 'slope.tif',format= 'GTiff',overwrite = True)
    grass.run_command('r.out.gdal', input = 'aspect',output = outputFolder  + 'aspect.tif',format= 'GTiff',overwrite = True)
    
    grass.run_command('v.out.ogr', input = 'str_finalcat',output = outputFolder  + 'str_finalcat.shp',format= 'ESRI_Shapefile',overwrite = True)


#####
    tempinfo = Dbf5(outputFolder + "/"+'Hylake.dbf')#np.genfromtxt(hyinfocsv,delimiter=',')
    allLakinfo = tempinfo.to_dataframe()
    landuseinfo = pd.read_csv(outputFolder + "/"+'landuseinfo.csv',sep=",",low_memory=False)
    
    finalcat_arr = garray.array(mapname="finalcat")
    acc_array = garray.array(mapname="acc_grass")
    dir_array = garray.array(mapname="dir_Arcgis")#ndir_Arcgis
    Lake1_arr = garray.array(mapname="SelectedLakes")    
    dem_array = garray.array(mapname="dem")
    rivlen_array = garray.array(mapname="Length")
    area_array = garray.array(mapname="Area")
    width_array = garray.array(mapname="width")
    depth_array = garray.array(mapname="depth")
    obs_array = garray.array(mapname="obs")
    Q_mean_array = garray.array(mapname="qmean")
    slope_array = garray.array(mapname="tanslopedegree")
    landuse_array = garray.array(mapname="landuse")
    
    
    temparray = garray.array()
    temparray[:,:] = -9999
    ncols = int(temparray.shape[1])
    nrows = int(temparray.shape[0])
    

    allcatid = np.unique(finalcat_arr)
    allcatid = allcatid[allcatid >= 0]
    catinfo2 = np.full((len(allcatid),18),-9999.00000)
    
    catinfodf = pd.DataFrame(catinfo2, columns = ['SubId', "DowSubId","Rivlen",'RivSlope','BasinSlope',
                            'BkfWidth','BkfDepth','IsLake','HyLakeId','LakeVol','LakeDepth',
                             'LakeArea','Laketype','IsObs','MeanElev','FloodP_n','Q_Mean','Ch_n'])
                                 
    catinfo,rivpath= Generatecatinfo(finalcat_arr,acc_array,dir_array,Lake1_arr,dem_array,
             area_array,catinfodf,allcatid,allLakinfo,width_array,depth_array,rivlen_array,obs_array,nrows,ncols,
             slope_array,landuse_array,landuseinfo,Q_mean_array)
             
    catinfo.to_csv (os.path.join(outputFolder, "catinfo.csv"), index = None, header=True)
    
    
    temparray[:,:] = rivpath[:,:]
    temparray.write(mapname="rivpath", overwrite=True)
    grass.run_command('r.null', map='rivpath',setnull=-9999)
    grass.run_command('db.in.ogr', input=os.path.join(outputFolder, "catinfo.csv"),output = 'result',overwrite = True)
    grass.run_command('v.db.join', map= 'finalcat_F',column = 'Gridcode', other_table = 'result',other_column ='SubId', overwrite = True)
    grass.run_command('v.out.ogr', input = 'finalcat_F',output = outputFolder  + 'finalcat_info.shp',format= 'ESRI_Shapefile',overwrite = True)
    PERMANENT.close
    
    
    
    