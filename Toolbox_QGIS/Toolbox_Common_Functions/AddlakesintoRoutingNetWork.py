import copy
import numpy as np
import os
from GetBasinoutlet import Defcat

def Generate_stats_list_from_grass_raster(grass,mode = 1,input_a = 'P_Connect_Lake',
                                          input_b = 'str_grass_r',method='max'):
    
    list_a = []
    list_b = []
    if mode == 1:
        p = grass.pipe_command('r.category', map = input_a, separator = ' ', quiet = True)
    elif mode == 2:
        p = grass.pipe_command('r.stats', input = [input_a,input_b],flags='n', quiet = True)
    else:
        print('error mode opton did not exist')
    
    for line in p.stdout:
        line_str = line.decode('utf8').rstrip('\r\n').split(' ')
        list_a.append(int(line_str[0]))
        if mode != 1:
            list_b.append(int(line_str[1]))
    p.wait()
    
    return list_a,list_b
    
def Check_If_Lake_Have_Multi_OutLet(CL_Id,Str_Id,Routing_info):
    
    ### create a emppty array to store lake id with multi outlet 
    Lakes_WIth_Multi_Outlet = []
    Remove_Str = []
    #### combine lake id and coresponding str id into two dim array 
    CL_Str = np.column_stack((CL_Id, Str_Id))
    #### obtain unique str id 
    Str_Id_unique = np.unique(np.array(Str_Id))
    #### obtain unique lake ids 
    CL_Id_unique = np.unique(np.array(CL_Id))
    
    #### add maxacc of each str to CL_str 
    Acc_Str_CL = np.full((len(CL_Str),1),np.nan)
    for i in range(0,len(Str_Id_unique)):
        strid  = Str_Id_unique[i]
        maxacc = Routing_info.loc[Routing_info['SubId'] ==strid]['MaxAcc_cat'].values[0]
        Acc_Str_CL[CL_Str[:,1] == strid] = maxacc
    CL_Str = np.append(CL_Str, Acc_Str_CL, axis=1)
    
    ### sort CL_str based on max acc of 
    CL_Str = CL_Str[CL_Str[:,2].argsort()]
    
    
    ### routing info 
    routing_info_only      = Routing_info[['SubId','DowSubId']].astype('int').values
#    print(routing_info_only)
    ### check if lakes have multi outlet 
    for i in range(0,len(CL_Id_unique)):
        ### got lake id 
        lake_id = CL_Id_unique[i]
        ### got all str overlaied by this lake id 
        i_CL_Str = CL_Str[CL_Str[:,0] == lake_id]
        
        ### lake is overlaied by one str, so must only have one lake outlet 
        if len(i_CL_Str) <= 1:
            continue 
        
        ### len(i_CL_Str)>1 
        ### check if all str covered by this lake will drainage to the str with 
        ### maximum acc. if it is ture than lake has only one outelt,
        ### if not the lake have multi outlet
        
        ###obtain str id with max acc among strs covered by the lake 
        str_max_acc = i_CL_Str[len(i_CL_Str) - 1,1]
#        print(str_max_acc,len(routing_info_only),"#################")
        ### obtain all upstream str id of max acc str 
        str_to_str_max_acc  = Defcat(routing_info_only,str_max_acc)
        
        if len(str_to_str_max_acc) <= 1: 
            str_to_str_max_acc  = Defcat(routing_info_only,i_CL_Str[len(i_CL_Str) - 2,1])
        ### create a mask for i_CL_Str[:,1], it will be true for in positon 
        ### where it's value in str_to_str_max_acc
        mask = np.isin(i_CL_Str[:,1],str_to_str_max_acc)
        #### not all element in i_CL_Str[:,1] exist in str_to_str_max_acc
        #### the lake have muli outlet 

        Remove_Str_i = []
        
        ### if not all str drainage to outlet str 
        if len(mask[mask == True]) < len(i_CL_Str[:,1]):
            ### obtain strids of str not drainage to outlet str 
            str_notflowto_lakeoutlet = i_CL_Str[np.logical_not(mask),1]
            ### loop for strids of str not drainage to outlet str
            for istr in range(0,len(str_notflowto_lakeoutlet)):
                ### get i str id 
                strid =  str_notflowto_lakeoutlet[istr]
                ### check  streams drainage to current str
#                upstrs  = Defcat(routing_info_only,strid)
                #### no upstream str, ### remove str instead remove lake 
#                if len(upstrs) == 1:
                Remove_Str_i.append(strid)
            
            if len(Remove_Str_i) == len(str_notflowto_lakeoutlet):
                Remove_Str = Remove_Str + Remove_Str_i
            else:
                Lakes_WIth_Multi_Outlet.append(lake_id)
    return Lakes_WIth_Multi_Outlet,Remove_Str
    
                
def DefineConnected_Non_Connected_Lakes(self,grass,con,garray,Routing_info,str_r = 'str_grass_r',
    Lake_r = 'alllake',Remove_Lake_Multi_OL = True):
    
    #### overlay str and lakes 
    exp = "Connect_Lake_IDs = if(isnull('%s'),null(),%s)" %(str_r,Lake_r)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)    
    
    ### obtain connected lake ids 
    Connect_Lake_Ids,temp = Generate_stats_list_from_grass_raster(grass,mode = 1,input_a = 'Connect_Lake_IDs')
    
    #### create non connected lake raster  
    grass.run_command('g.copy',rast = (Lake_r,'Nonconnect_Lake'),overwrite = True)
    grass.run_command('r.null',map='Nonconnect_Lake', setnull = Connect_Lake_Ids,overwrite = True)
    #### create potential connected lake raster 
    exp = "Connect_Lake = if(isnull('%s'),%s,null())" %('Nonconnect_Lake',Lake_r)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)   
    
    if Remove_Lake_Multi_OL == False:
        return 
    
    
    ##### obtain lake id and coresponding overlaied str id 
    CL_Id,Str_Id = Generate_stats_list_from_grass_raster(grass,mode = 2,input_a = 'Connect_Lake',
                                                         input_b = str_r)
    
    Lakes_WIth_Multi_Outlet,Remove_Str = Check_If_Lake_Have_Multi_OutLet(CL_Id,Str_Id,Routing_info)
    ### no lake has multi outlet 

    if len(Lakes_WIth_Multi_Outlet) > 0: 
#        grass.run_command('r.null',map='Connect_Lake', setnull = Lakes_WIth_Multi_Outlet,overwrite = True)
        print("Following lakes have identified outlet :       ")
        print("     ", Lakes_WIth_Multi_Outlet)
    if len(Remove_Str) > 0:
#        grass.run_command('r.null',map='str_grass_r', setnull = Remove_Str,overwrite = True)
#        grass.run_command('r.null',map='cat1', setnull = Remove_Str,overwrite = True)
        print("Following stream have been identified to make each lake has one lake outlet:       ")
        print("     ", Remove_Str)        
    return Remove_Str
    
def changeflowdirectionofedgegrids(N_dir,p_row,p_col,lake1,lid,ncols,nrows,BD_Out_Lakecat_Nriv_mask,Changed_ndir):
    ndir = copy.copy(N_dir)
    changed_ndir = copy.copy(Changed_ndir)
    IS_Change = 0

    if lake1[p_row + 0,p_col + 1] == lid:  #### it is a broudary grids and did not flow to the lake outlet
        ndir[p_row,p_col] = 1  ### change the flow direction of new point to old points
        changed_ndir[p_row,p_col] = 1
        IS_Change = 1
#dir 2
    if lake1[p_row + 1,p_col + 1] == lid:

        ndir[p_row,p_col] = 2
        changed_ndir[p_row,p_col] = 2
        IS_Change = 1
#dir 3
    if lake1[p_row + 1,p_col + 0] == lid:

        ndir[p_row,p_col] = 4
        changed_ndir[p_row,p_col] = 4
        IS_Change = 1
#dir 4
    if lake1[p_row + 1,p_col - 1] == lid:

        ndir[p_row,p_col] = 8
        changed_ndir[p_row,p_col] = 8
        IS_Change = 1

#dir 5
    if lake1[p_row + 0,p_col - 1] == lid:

        ndir[p_row,p_col] = 16
        changed_ndir[p_row,p_col] = 16
        IS_Change = 1

#dir 6
    if lake1[p_row - 1,p_col - 1] == lid:

        ndir[p_row,p_col] = 32
        changed_ndir[p_row,p_col] = 32
        IS_Change = 1
#dir 7
    if lake1[p_row - 1,p_col + 0] == lid:

        ndir[p_row,p_col] = 64
        changed_ndir[p_row,p_col] = 64
        IS_Change = 1
#dir 8
    if lake1[p_row - 1,p_col + 1] == lid:

        ndir[p_row,p_col] = 128
        changed_ndir[p_row,p_col] = 128
        IS_Change = 1
#    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir,IS_Change,changed_ndir



def Dirpoints_v3(N_dir,p_row,p_col,lake1,lid,goodpoint,k,ncols,nrows,BD_Out_Lakecat_Nriv_mask,Changed_ndir):
    ### this function change flow direction of some grids around p_row, p_col, flow to p_row, p_col

    ndir = copy.copy(N_dir)
    changed_ndir = copy.copy(Changed_ndir)
    ip = copy.copy(k) + 1
#dir 1
    if BD_Out_Lakecat_Nriv_mask[p_row + 0,p_col + 1] == lid:  #### it is a broudary grids and did not flow to the lake outlet
        tt = goodpoint[goodpoint[:,0] == p_row + 0,]
        if len(tt[tt[:,1] == p_col + 1,]) < 1:#### the point not exist in good points which store all boundary points
            ndir[p_row + 0,p_col + 1] = 16   ### change the flow direction of new point to old points
            changed_ndir[p_row + 0,p_col + 1] = 16
            goodpoint[ip,0] = p_row + 0
            goodpoint[ip,1] = p_col + 1
            ip = ip + 1
#dir 2
    if BD_Out_Lakecat_Nriv_mask[p_row + 1,p_col + 1] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row + 1,]
        if len(tt[tt[:,1] == p_col + 1,]) < 1:
            ndir[p_row + 1,p_col + 1] = 32
            changed_ndir[p_row + 1,p_col + 1] = 32
            goodpoint[ip,0] = p_row + 1
            goodpoint[ip,1] = p_col + 1
            ip = ip + 1
#dir 3
    if BD_Out_Lakecat_Nriv_mask[p_row + 1,p_col + 0] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row + 1,]
        if len(tt[tt[:,1] == p_col + 0,]) < 1:
            ndir[p_row + 1,p_col + 0] = 64
            changed_ndir[p_row + 1,p_col + 0] = 64
            goodpoint[ip,0] = p_row + 1
            goodpoint[ip,1] = p_col + 0
            ip = ip + 1
#dir 4
    if BD_Out_Lakecat_Nriv_mask[p_row + 1,p_col - 1] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row + 1,]
        if len(tt[tt[:,1] == p_col - 1,]) < 1:
            ndir[p_row + 1,p_col - 1] = 128
            changed_ndir[p_row + 1,p_col - 1] = 128
            goodpoint[ip,0] = p_row + 1
            goodpoint[ip,1] = p_col - 1
            ip = ip + 1
#dir 5
    if BD_Out_Lakecat_Nriv_mask[p_row + 0,p_col - 1] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row + 0,]
        if len(tt[tt[:,1] == p_col - 1,]) < 1:
            ndir[p_row + 0,p_col - 1] = 1
            changed_ndir[p_row + 0,p_col - 1] = 1
            goodpoint[ip,0] = p_row + 0
            goodpoint[ip,1] = p_col - 1
            ip = ip + 1
#dir 6
    if BD_Out_Lakecat_Nriv_mask[p_row - 1,p_col - 1] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row - 1,]
        if len(tt[tt[:,1] == p_col - 1,]) < 1:
            ndir[p_row - 1,p_col - 1] = 2
            changed_ndir[p_row - 1,p_col - 1] = 2
            goodpoint[ip,0] = p_row - 1
            goodpoint[ip,1] = p_col - 1
            ip = ip + 1
#dir 7
    if BD_Out_Lakecat_Nriv_mask[p_row - 1,p_col + 0] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row - 1,]
        if len(tt[tt[:,1] == p_col + 0,]) < 1:
            ndir[p_row - 1,p_col + 0] = 4
            changed_ndir[p_row - 1,p_col + 0] = 4
            goodpoint[ip,0] = p_row - 1
            goodpoint[ip,1] = p_col + 0
            ip = ip + 1
#dir 8
    if BD_Out_Lakecat_Nriv_mask[p_row - 1,p_col + 1] == lid:

        tt = goodpoint[goodpoint[:,0] == p_row - 1,]
        if len(tt[tt[:,1] == p_col + 1,]) < 1:
            ndir[p_row - 1,p_col + 1] = 8
            changed_ndir[p_row - 1,p_col + 1] = 8
            goodpoint[ip,0] = p_row - 1
            goodpoint[ip,1] = p_col + 1
            ip = ip + 1
#    arcpy.AddMessage(goodpoint[goodpoint[:,0]>0,])
    return ndir,goodpoint,ip,changed_ndir


def check_lakecatchment(cat3,lake,fac,fdir,bsid,nrows,ncols,LakeBD_array,nlakegrids,str_array,dir,Pec_Grid_outlier,MaximumLakegrids,Lakemorestream):
    cat = copy.copy(cat3)
    ndir = copy.copy(dir)
    changed_ndir = copy.copy(dir)
    changed_ndir[:,:] = -9999
    BD_problem = copy.copy(dir)
    BD_problem[:,:] = -9999
    arlakeid = np.unique(lake)
    arlakeid = arlakeid[arlakeid>0]
    outlakeids = np.full((1000000,2),-99999.999)
    stream_mask = str_array > 0
    for i in range(0,len(arlakeid)):
        lakeid = arlakeid[i]
        if len(Lakemorestream) >0:
            testmaks = Lakemorestream == lakeid
            if np.sum(testmaks)>=1:
                continue

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
#        print("########################################################################3")
#        print(lakeid,arclakeid,len(nlake),len(lrowcol),float(len(nlake)/len(lrowcol)),Pec_Grid_outlier)
        if outlakeids[i,1] > Pec_Grid_outlier:
            continue

#        if outlakeids[i,1] > Pec_Grid_outlier: ### smaller than 0.97
#            continue

#        print(lakeid,arclakeid,len(nlake),len(lrowcol),float(len(nlake)/len(lrowcol)))
        BD_mask        = LakeBD_array == lakeid
        Lake_mask      = lake == lakeid
        cat_lake_mask  = cat  == arclakeid

        Lakeincat_mask = np.logical_and(Lake_mask,cat_lake_mask)
        nlake2         = np.sum(Lakeincat_mask)

        Lakeoutcat_mask = np.logical_and(Lake_mask,np.logical_not(cat_lake_mask))

        BD_Out_Lakecat_mask = np.logical_and(Lakeoutcat_mask,BD_mask)

#        BD_Out_Lakecat_Nriv_mask = np.logical_and(BD_Out_Lakecat_mask,np.logical_not(stream_mask)) do not allow modify stream grids
        BD_Out_Lakecat_Nriv_mask = BD_Out_Lakecat_mask
        BD_problem[BD_Out_Lakecat_Nriv_mask] = 1

        if np.sum(BD_Out_Lakecat_mask)  > MaximumLakegrids:   ### smaller than nlakegrids or smaller than 0.9
            continue

        print("Total # of lakes: ", len(arlakeid), " Processing    ",i,"th lake" )
        print("Lake ID : ",lakeid,"Lake Cat ID   ",arclakeid, "Total numer of Lake grids   ",len(lrowcol), "Numer of Lake grids in Lake Cat:  ",nlake2,len(nlake))
        print("# of Lake boundary grids:   ", np.sum(BD_mask),"# of grids do not flow to lake catchments",np.sum(Lakeoutcat_mask) ,"# of lake boundary grids not flow to lake catchment   ", np.sum(BD_Out_Lakecat_mask),"  # of lake boundary grids not flow to lake catchment not a river gird  ", np.sum(BD_Out_Lakecat_Nriv_mask))


        #####  Locate the grids that at the target grids ege
        Grid_Nee_Mo = np.argwhere(BD_Out_Lakecat_Nriv_mask == 1)
#        print(print(Grid_Nee_Mo[:,0], Grid_Nee_Mo[:,1]))

        goodpoint = copy.copy(Grid_Nee_Mo)
        goodpoint_2 = copy.copy(Grid_Nee_Mo)
        goodpoint_2[:,:] =  -9999
        idx = 0
        for i in range(0,len(goodpoint)):
            trow = goodpoint[i,0]
            tcol = goodpoint[i,1]
#                print(i,ipo,trow,tcol)
            if trow >= nrows - 1 or tcol == ncols  - 1:
                continue

            ndir,IS_Change,changed_ndir=changeflowdirectionofedgegrids(ndir,trow,tcol,Lakeincat_mask,1,ncols,nrows,BD_Out_Lakecat_Nriv_mask,changed_ndir)
            if IS_Change > 0:
                goodpoint_2[idx] = goodpoint[i,:]
                idx = idx + 1
        goodpoint = goodpoint_2
#        print(goodpoint[goodpoint[:,0]>0,])
        ### add just these flow direction of edge grids to lake domian


        ###3 make BD_Out_Lakecat_Nriv_mask flow to Lakeoutcat_mask
#        print(goodpoint[goodpoint[:,0]>0,])
        ip = len(goodpoint[goodpoint[:,0]>0,])
        k = len(goodpoint[goodpoint[:,0]>0,]) - 1
        ipo = 0
        iter = 0
#        print("############################333")
        while ip > ipo and iter < np.sum(BD_Out_Lakecat_mask) + 10: ###maixmum loop for N points
            iter = iter + 1
#            print(ip,ipo,iter,np.sum(BD_Out_Lakecat_mask),"################3")
            for i in range(ipo,len(goodpoint[goodpoint[:,0]>0,])):
                trow = goodpoint[i,0]
                tcol = goodpoint[i,1]
#                print(i,ipo,trow,tcol)
                if trow >= nrows - 1 or tcol == ncols  - 1:
                    continue
                ndir,goodpoint,k1,changed_ndir= Dirpoints_v3(ndir,trow,tcol,Lakeincat_mask,1,goodpoint,k,ncols,nrows,BD_Out_Lakecat_Nriv_mask,changed_ndir)
                k = k1 - 1

            ipo = ip
            ip = len(goodpoint[goodpoint[:,0]>0,]) - 1
            k = ip

#        print(len(goodpoint[goodpoint[:,0]>0,]),np.sum(BD_Out_Lakecat_mask))

#        print(lakeid,arclakeid,len(nlake),len(lrowcol),float(len(nlake)/len(lrowcol)))
    outlakeids= outlakeids[outlakeids[:,0] > 0]
    return outlakeids,changed_ndir,ndir,BD_problem


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
