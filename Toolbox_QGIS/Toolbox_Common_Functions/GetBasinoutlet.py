import numpy as np
import copy
import pandas as pd 

def Generate_Routing_structure(grass,con,cat = 'Net_cat',acc = 'acc_grass',Name = 'a1',str = '#'):
    
    Raster_res = grass.core.region()['nsres']
    
    ###find outlet point of each cat 
    grass.run_command('r.stats.zonal', base=cat,cover=acc, method='max', output=Name+'_maxacc', overwrite = True)
    ### Find the grid that equal to the max acc, thus this is the outlet grids  
    exp = "'%s' =if(%s == int('%s'),%s,null())" %(Name+'_OL',acc,Name+'_maxacc',cat)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)
    ### change outlet points to point vector 
    grass.run_command('r.to.vect',input = Name+'_OL', output = Name+'_OL_v',type = 'point',flags = 'v',overwrite = True)
    grass.run_command('v.db.addcolumn', map=  Name+'_OL_v', columns = "SubId int")
    grass.run_command('v.db.update', map=  Name+'_OL_v', column = "SubId",qcol = 'cat')
    grass.run_command('v.what.rast', map= Name+'_OL_v', raster=acc,column='OL_acc') 
    grass.run_command('v.buffer',input =  Name+'_OL_v', output = Name+'_OL_v_bf',distance = 1.5 * Raster_res,flags = 't', overwrite = True)
    ### 
    ###find inlet  point of each cat 
    exp = "'%s' =if(isnull(%s),null(),%s)" %(Name+'_acc_riv',str,acc)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)
    grass.run_command('r.null', map=Name+'_acc_riv',setnull=[-9999,0])
    grass.run_command('r.stats.zonal', base=str,cover=Name+'_acc_riv', method='min', output=Name+'_minacc', overwrite = True)
    ### Find the grid that equal to the max acc, thus this is the outlet grids  
    exp = "'%s' =if(%s == int('%s'),%s,null())" %(Name+'_IL',Name+'_acc_riv',Name+'_minacc',str)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)
    ### change outlet points to point vector 
    grass.run_command('r.to.vect',input = Name+'_IL', output = Name+'_IL_v',type = 'point',flags = 'v',overwrite = True)
    grass.run_command('v.db.addcolumn', map=  Name+'_IL_v', columns = "IL_SubId int")
#    grass.run_command('v.db.addcolumn', map=  Name+'_IL_v', columns = "IL_acc int")
    grass.run_command('v.db.update', map=  Name+'_IL_v', column = "IL_SubId",qcol = 'cat')
    grass.run_command('v.what.rast', map=Name+'_IL_v', raster=acc,column='IL_acc') 
#    grass.run_command('v.db.update', map=  Name+'_IL_v', column = "IL_acc",qcol = 'IL_acc_db - 1')
    
    #### add max acc inlet point acc value to Name+'_OL_v_bf'
    grass.run_command('v.vect.stats', point =Name+'_IL_v',areas =Name+'_OL_v_bf',method='max_cat',points_column = 'IL_acc',count_column='IL_Pt_CT',stats_column = 'max_cat')
    grass.run_command('v.vect.stats', point =Name+'_IL_v',areas =Name+'_OL_v_bf',method='min_cat',points_column = 'IL_acc',count_column='IL_Pt_CT2',stats_column = 'min_cat')
    grass.run_command('v.vect.stats', point =Name+'_IL_v',areas =Name+'_OL_v_bf',method='maximum',points_column = 'IL_acc',count_column='IL_Pt_CT3',stats_column = 'ILmaxacc')
    grass.run_command('v.vect.stats', point =Name+'_IL_v',areas =Name+'_OL_v_bf',method='minimum',points_column = 'IL_acc',count_column='IL_Pt_CT4',stats_column = 'ILminacc')
    
    #### add outlet's next inlet subid into Name+'_OL_v_bf'
    grass.run_command('v.db.join', map= Name+'_OL_v_bf',column = 'max_cat', other_table = Name+'_IL_v',other_column ='cat', subset_columns = 'IL_SubId',overwrite = True)
    grass.run_command('v.db.addcolumn', map= Name+'_OL_v_bf', columns = "ILSubIdmax int")
    grass.run_command('v.db.addcolumn', map= Name+'_OL_v_bf', columns = "ILAccstrmin int")
    
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "ILSubIdmax",qcol = 'IL_SubId')
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "ILAccstrmin",qcol = 'ILmaxacc')
    ### add min next inlet subid to Name+'_OL_v_bf'
    grass.run_command('v.db.join', map= Name+'_OL_v_bf',column = 'min_cat', other_table = Name+'_IL_v',other_column ='cat', subset_columns = 'IL_SubId',overwrite = True)
    grass.run_command('v.db.addcolumn', map= Name+'_OL_v_bf', columns = "ILSubIdmin int")
    grass.run_command('v.db.addcolumn', map= Name+'_OL_v_bf', columns = "ILAccstrmax int")
    
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "ILSubIdmin",qcol = 'IL_SubId')
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "ILAccstrmin",qcol = 'ILminacc')
    
    grass.run_command('v.db.addcolumn', map= Name+'_OL_v_bf', columns = "DSubId_str int")
    grass.run_command('v.db.addcolumn', map= Name+'_OL_v_bf', columns = "DSubAccstr int")
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "DSubId_str",qcol = 'ILSubIdmax')
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "DSubAccstr",qcol = 'ILmaxacc')
    
    ### in case cat drainage to itself, drainage to another close subbasin in inlet 
    exp = " '%s' = '%s' " %('SubId','ILSubIdmax')
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "DSubId_str",qcol = 'ILSubIdmin',where = 'SubId = ILSubIdmax')
    grass.run_command('v.db.update', map=  Name+'_OL_v_bf', column = "DSubAccstr",qcol = 'ILminacc',where = 'SubId = ILSubIdmax')
    
    ### add downsubid_Str to point vector 
    grass.run_command('v.db.join', map=  Name+'_OL_v',column = 'SubId', other_table = Name+'_OL_v_bf',other_column ='SubId', subset_columns = ['DSubId_str','IL_Pt_CT','ILSubIdmin','ILSubIdmax','DSubAccstr'],overwrite = True)
 
    ### find downn subid via cat, the reason needs to do this, is that for some cat 
    ### without connected river, above approach can not find downstream cat id for 
    ### these catchment 
    
    ###add max acc in of neighbor grids as each growed grids value  
    exp = "'%s'=if(isnull('%s'),null(),1)" %(Name+'_OL1',Name+'_OL')
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)
    grass.run_command('r.grow', input=Name+'_OL1', output=Name+'_OL1_G', radius=1.5,overwrite = True)  
    grass.run_command('r.clump', input=Name+'_OL1_G', output=Name+'_OL1_G_Clu',overwrite = True)   
        
    ###find inlets of each subbasin 
    grass.run_command('r.stats.zonal', base=Name+'_OL1_G_Clu',cover=acc, method='max', output=Name+'_OL1_G_Clu_maxacc', overwrite = True)
    exp = "'%s'=if(%s == int('%s'),%s,null())" %(Name+'_IL',acc,Name+'_OL1_G_Clu_maxacc',cat)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True) 
    grass.run_command('r.stats.zonal', base=Name+'_OL1_G_Clu',cover=Name+'_IL', method='max', output=Name+'_OL1_G_Clu_IL_SubId', overwrite = True)
        
    grass.run_command('v.what.rast', map=Name+'_OL_v', raster=Name+'_OL1_G_Clu_IL_SubId',column='DSubId_cat')
    grass.run_command('v.what.rast', map=Name+'_OL_v', raster=Name+'_maxacc',column='MaxAcc_cat')
#    grass.run_command('v.db.join', map=  Name+'_OL_v',column = 'SubId', other_table = Name+'_OL_v_bf',other_column ='SubId', subset_columns = ['DSubId_str','OL_acc','IL_Pt_CT','max_cat','IL_SubId_min'],overwrite = True)


    #### create the final downsubid use downsubid from str if it is exist 
    grass.run_command('v.db.addcolumn', map=  Name+'_OL_v', columns = "DowSubId int")
    grass.run_command('v.db.update', map=  Name+'_OL_v', column = "DowSubId", where = 'IL_Pt_CT > 0', qcol = 'DSubId_str')
    grass.run_command('v.db.update', map=  Name+'_OL_v', column = "DowSubId", where = 'IL_Pt_CT <= 0', qcol = 'DSubId_cat')
    ### if it move to itself using dowsubid from cat 
    grass.run_command('v.db.update', map=  Name+'_OL_v', column = "DowSubId", where = 'SubId = DSubId_str', qcol = 'DSubId_cat')
    ### if it move to a cat with smaller acc, use down subid from cat 
    grass.run_command('v.db.update', map=  Name+'_OL_v', column = "DowSubId", where = 'OL_acc > DSubAccstr', qcol = 'DSubId_cat')
    sqlstat="SELECT SubId, DowSubId, MaxAcc_cat FROM %s" % (Name+'_OL_v')
    Routing_info = pd.read_sql_query(sqlstat, con)
    return Routing_info
    
def Defcat(out,outletid):
    otsheds = np.full((1,1),outletid)
    Shedid = np.full((len(out)+1,1),-99999999999)
    psid = 0
    rout = copy.copy(out)
    while len(otsheds) > 0:
        noutshd = np.full((len(out)+1,1),-99999999999)
        poshdid = 0
#        print("################################################a")
        for i in range(0,len(otsheds)):
#            print(otsheds)
#            print(psid,outletid)
            Shedid[psid] = otsheds[i]
#            print(Shedid[psid],otsheds[i])
#            print("##################################################b")
            psid = psid + 1
            irow = np.argwhere(rout[:,1]==otsheds[i]).astype(int)
#            print(len(irow))
            for j in range(0,len(irow)):
                #### if the catchment id already processed skip
                if rout[irow[j],0] in Shedid:
                    continue
                noutshd[poshdid] = rout[irow[j],0]
                poshdid = poshdid + 1
        noutshd = np.unique(noutshd)
        otsheds = noutshd[noutshd>=0]
    Shedid = np.unique(Shedid)
    Shedid = Shedid[Shedid>=0]
    return Shedid
###########



def Getbasinoutlet(ID,basin,fac,dir,nrows,ncols):
    catrowcol = np.argwhere(basin==ID).astype(int)
    catacc = np.full((len(catrowcol),3),-9999)
    catacc[:,0] = catrowcol[:,0]
    catacc[:,1] = catrowcol[:,1]
    catacc[:,2] = fac[catrowcol[:,0],catrowcol[:,1]]
    catacc = catacc[catacc[:,2].argsort()]
#    print(ID,len(catrowcol))
    ### check if it is a real basin outlet 
    crow = catacc[len(catrowcol)-1,0]
    ccol = catacc[len(catrowcol)-1,1]
    
    nrow,ncol =  Nextcell(dir,crow,ccol)
#    print(ID,basin[nrow,ncol],basin[crow,ccol],fac[nrow,ncol],fac[crow,ccol],crow,ccol)
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

