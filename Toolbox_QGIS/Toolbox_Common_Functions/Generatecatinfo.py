
import numpy as np
from GetBasinoutlet import Getbasinoutlet,Nextcell
from Calculate_River_Len_Slope import Getcatrivlenslope_hydroshed
import copy


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
    if (V > 0):
        n = (Rch**(2.0/3.0))*(slope**(1.0/2.0))/V
    else:
        n = 0.0012345
    return n



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
        catrivlen,catrivslp,catrivslp2,rivpath = Getcatrivlenslope_hydroshed(rowcol[:,0],rowcol[:,1],rivlen,dem,fac,fdir,finalcat,
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
            catinfo.loc[i,'RivSlope'] = 0.0001234
        else:
            catinfo.loc[i,'RivSlope'] = catrivslp
#        catinfo[i,27] = catrivslp2
        if len(slopet) < 1:
            catinfo.loc[i,'RivSlope'] = 0.0001234
            catinfo.loc[i,'BasinSlope'] = 0.0001234
        catinfo.loc[i,'FloodP_n'] = Getfloodplain_n(catid,finalcat,rivlen,landuse,landuseinfo)
        catinfo.loc[i,'Q_Mean'] = catQ
#    writeraster(outputFolder + "/" + "rivpath.asc",rivpath,OutputFolder + "/" + "dir")
    return catinfo,rivpath


     

def Generatecatinfo_riv(Watseds,fac,fdir,lake,dem,catinfo,allcatid,width,depth,
                    obs,slope,aspect,landuse,slop_deg,Q_Mean,landuseinfo,lakeinfo,
                    nrows,ncols,leninfo,areainfo):
    finalcat = copy.copy(Watseds)
    for i in range(0,len(allcatid)):
        catid = allcatid[i].astype(int)
        print(i,catid)
        catinfo.loc[i,'SubId'] = catid
        catmask = Watseds == catid
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

################################## Get lake information        
        lakeinriv = lake[catmask]
        lakeids = np.unique(lakeinriv)
        lakeids = lakeids[lakeids > 0]
        if len(lakeids) == 1:
            lakeid = lakeids[0]
        elif len(lakeids) > 1:
            print('Warning:  stream    ',catid,'connected with ',len(lakeids),'   Lakes')
            lakeid = lakeids[0]
            for j in range(1,len(lakeids)):
                if len(np.argwhere(lakeinriv == lakeid)) < len(np.argwhere(lakeinriv == lakeids[j])):
                    lakeid = lakeids[j]
        else:
            lakeid = -1 
            
        if lakeid > 0:
            slakeinfo = lakeinfo.loc[lakeinfo['Hylak_id'] == lakeid]
            catinfo.loc[i,'IsLake'] = 1
            catinfo.loc[i,'HyLakeId'] = lakeid
            catinfo.loc[i,'LakeVol'] = slakeinfo.iloc[0]['Vol_total']
            catinfo.loc[i,'LakeArea']= slakeinfo.iloc[0]['Lake_area']
            catinfo.loc[i,'LakeDepth']= slakeinfo.iloc[0]['Depth_avg']
            catinfo.loc[i,'Laketype'] = slakeinfo.iloc[0]['Lake_type']
########Check if it is observation points
        if obs[trow,tcol]  >= 0:
#            arcpy.AddMessage(str(catid)+"      "+str(obs[trow,tcol]))
            catinfo.loc[i,'IsObs'] =  obs[trow,tcol]


########Slopes slope,aspect,landuse,slop_deg
        slopeinriv = slope[catmask]
        aspectinriv = aspect[catmask]
        slop_deginriv = slop_deg[catmask]
        deminriv = dem[catmask]
        
        if(len(slop_deginriv[slop_deginriv > 0])) > 0:
            slop_deginriv[slop_deginriv <=0] = np.NaN  
            catinfo.loc[i,'BasSlope'] = np.nanmean(slop_deginriv)
        else:
            catinfo.loc[i,'BasSlope'] = 1.2345      

        if(len(aspectinriv[aspectinriv > 0])) > 0:
            aspectinriv[aspectinriv <=0] = np.NaN  
            catinfo.loc[i,'BasAspect'] = np.nanmean(aspectinriv)
        else:
            catinfo.loc[i,'BasAspect'] = 1.2345  
        
        if(len(deminriv[deminriv > 0])) > 0:
            deminriv[deminriv <=0] = np.NaN
            maxdem = np.nanmax(deminriv)
            mindem = np.nanmin(deminriv)
            catinfo.loc[i,'MeanElev'] =np.nanmean(deminriv)
        else:
            maxdem = 1.2345
            mindem = 1.2345
            
        rivlen = np.unique(leninfo.loc[leninfo['Gridcode'] == catid]['Length_m'].values)  #'Area_m'
        if len(rivlen) == 1:
            catinfo.loc[i,'RivLength'] = rivlen
            if rivlen > 0:
                if max(0,float((maxdem - mindem))/float(rivlen)) == 0:
                    catinfo.loc[i,'RivSlope'] = 0.0012345
                else:
                    catinfo.loc[i,'RivSlope'] = max(0,float((maxdem - mindem))/float(rivlen))
            else:
                catinfo.loc[i,'RivSlope'] = 0.0012345
        else:
            print("Warning  river length of stream  " , catid, "   need check   ", len(rivlen) )
            catinfo.loc[i,'RivLength'] = 1.2345
            catinfo.loc[i,'RivSlope'] = 0.0012345
            
            
        catarea = np.unique(areainfo.loc[areainfo['Gridcode'] == catid]['Area_m'].values)  #'Area_m'
        if len(catarea) == 1:
            catinfo.loc[i,'BasArea'] = catarea
        else:
            print("Warning  basin area of stream  " , catid, "   need check   ", len(catarea) )
            catinfo.loc[i,'BasArea'] = 1.2345

##########
        Landtypes = landuse[catmask]
        Landtypeid = np.unique(Landtypes)
        Landtypeid1 = Landtypeid[Landtypeid >= 0]
        Landtypeid2 = Landtypeid1[Landtypeid1 > 0]
        Landtypes = Landtypes[Landtypes > 0]
        if len(Landtypes) > 0 and float(len(Landtypeid2))/float(len(Landtypeid1)) >= 0.1:
            sum = 0.0
            for j in range(0,len(Landtypeid2)):
                iid = Landtypeid2[j]
                sum = sum + landuseinfo[landuseinfo['RasterV'] == iid]['MannV'].values*len(np.argwhere(Landtypes == iid))
            floodn = sum/(len(Landtypes))
        else:
            floodn = 0.035
        catinfo.loc[i,'FloodP_n'] = floodn
        
        
            
########Got basin width and depth
        widthinriv = width[catmask]
        depthinriv = depth[catmask]
        Q_Meaninriv = Q_Mean[catmask]
        
        widthids = np.unique(widthinriv)
        widthids = widthids[widthids > 0]
        print(i,catid)
        if(len(widthids)) > 0:
            widthinriv[widthinriv <=0] = np.NaN
            depthinriv[depthinriv <=0] = np.NaN
            Q_Meaninriv[Q_Meaninriv <=0] = np.NaN
            catinfo.loc[i,'BkfWidth'] = np.nanmean(widthinriv)
            catinfo.loc[i,'BkfDepth'] = np.nanmean(depthinriv)
            catinfo.loc[i,'Q_Mean'] = np.nanmean(Q_Meaninriv)
#            print(i,len(widthids),np.nanmean(widthinriv))           
        else:
            catinfo.loc[i,'BkfWidth'] = 1.2345
            catinfo.loc[i,'BkfDepth'] = 1.2345
            catinfo.loc[i,'Q_Mean'] =  1.2345
            print(i,len(widthids),np.nanmean(widthinriv))    
        
        if catinfo['BkfWidth'].values[i] != 1.2345 and catinfo['RivSlope'].values[i] != 0.0012345:
            catinfo.loc[i,'Ch_n']  = calculateChannaln(catinfo['BkfWidth'].values[i],catinfo['BkfDepth'].values[i],
                                                        catinfo['Q_Mean'].values[i],catinfo['RivSlope'].values[i])
        else:
            catinfo.loc[i,'Ch_n']  = 0.0012345

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

