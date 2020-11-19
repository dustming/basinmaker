import numpy as np
import copy
import pandas as pd 

def Generate_Routing_structure(grass,con,cat = 'str_grass_r',acc = 'acc_grass',Name = 'a1'):
    grass.run_command('r.stats.zonal', base=cat,cover=acc, method='max', output=Name+'_maxacc', overwrite = True)
    exp = "'%s' =if(%s == int('%s'),%s,null())" %(Name+'_OL',acc,Name+'_maxacc',cat)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)
    exp = "'%s'=if(isnull('%s'),null(),1)" %(Name+'_OL1',Name+'_OL')
    grass.run_command('r.mapcalc',expression = exp,overwrite = True)
    grass.run_command('r.grow', input=Name+'_OL1', output=Name+'_OL1_G', radius=1.5,overwrite = True)  
    grass.run_command('r.clump', input=Name+'_OL1_G', output=Name+'_OL1_G_Clu',overwrite = True)   
        
    ###find inlets of each subbasin 
    grass.run_command('r.stats.zonal', base=Name+'_OL1_G_Clu',cover=acc, method='max', output=Name+'_OL1_G_Clu_maxacc', overwrite = True)
    exp = "'%s'=if(%s == int('%s'),%s,null())" %(Name+'_IL',acc,Name+'_OL1_G_Clu_maxacc',cat)
    grass.run_command('r.mapcalc',expression = exp,overwrite = True) 
    grass.run_command('r.stats.zonal', base=Name+'_OL1_G_Clu',cover=Name+'_IL', method='max', output=Name+'_OL1_G_Clu_IL_SubId', overwrite = True)
        
    grass.run_command('r.to.vect', input=Name+'_OL', output=Name+'_outlet', type='point',overwrite = True)        
    grass.run_command('v.db.renamecolumn', map=Name+'_outlet', column='value,SubId')
    grass.run_command('v.db.dropcolumn', map=Name+'_outlet', column='label')
    grass.run_command('v.what.rast', map=Name+'_outlet', raster=Name+'_OL1_G_Clu_IL_SubId',column='DowSubId')
    grass.run_command('v.what.rast', map=Name+'_outlet', raster=Name+'_maxacc',column='MaxAcc')
    sqlstat="SELECT SubId, DowSubId, MaxAcc FROM %s" % (Name+'_outlet')
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

