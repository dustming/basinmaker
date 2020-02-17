import numpy as np
from GetBasinoutlet import Getbasinoutlet,Nextcell


def Checkcat2(prow,pcol,nrows,ncols,lid,lake):
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



def Getcatrivlenslope_hydroshed(catrow,catcol,rivlen,dem,fac,hydir,finalcat,trow,tcol,nrows,ncols,slope,rivpath):

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
        nout, nearcat = Checkcat2(prow,pcol,nrows,ncols,lid,finalcat)  #### check if the point (prow,pcol) is close to the catchment boundary, most upstream cell 
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
                    if nrow > 0 and ncol >0:  ##### if they are <0, means the downstream catchment is None, so it is not a error 
                        print("warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol)
                    nrow,ncol = -1,-1 ## to break the loop
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
                    if nrow > 0 and ncol >0:  ##### if they are <0, means the downstream catchment is None, so it is not a error 
                        print("warning : check river system for catchment: ",finalcat[trow,tcol],rivs,icell,len(catrow),nrow,ncol,trow,tcol)
                    nrow,ncol = -1,-1 ## to break the loop
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

def Getcatrivlenslope_other(finalcat,catid,rivlen,dem,fac,hydir,trow,tcol,nrows,ncols,slope,rivpath):
    
    catregs = finalcat == catid
    riverp = rivlen > 0
    rivincat = np.logical_and(catregs, riverp)
        
    
    return outrivlen, outrivslp,outrivslp2,rivpath

######################################################
