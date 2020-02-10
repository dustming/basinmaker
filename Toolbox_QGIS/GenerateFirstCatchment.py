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
    import numpy as np
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

def WatershedDiscretizationToolset(OutputFolder,Thresholdacc,GRASS_BIN):
    import os
    import sys
    import tempfile
    import shutil
    import numpy as np
    
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
    cellSize = 0.0041666667
    
    ###### run script
    grass.run_command('r.accumulate', direction='dir_Grass',format = '45degree',accumulation ='acc_grass',
                  stream = 'str_grass_v',threshold = 500, overwrite = True)

    grass.run_command('v.to.rast',input = 'str_grass_v',output = 'str_grass_r2',use = 'cat',overwrite = True)
    
    strtemp_array = garray.array(mapname="str_grass_r2")
    acc_array = garray.array(mapname="acc_grass")
    dirarc_array = garray.array(mapname="dir_Arcgis")
    
    strids = np.unique(strtemp_array)
    strids = strids[strids >= 0] 
    
    ncols = int(strtemp_array.shape[1])
    nrows = int(strtemp_array.shape[0])
 
    for i in range(0,len(strids)):
        strid = strids[i]

        trow,tcol = Getbasinoutlet(strid,strtemp_array,acc_array,dirarc_array,nrows,ncols)
        nrow,ncol = Nextcell(dirarc_array,trow,tcol)### get the downstream catchment id
        nstrid = strtemp_array[nrow,ncol]
        
        rowcol = np.argwhere(strtemp_array==nstrid).astype(int)
        catacc = np.full((len(rowcol),3),-9999)
        catacc[:,0] = rowcol[:,0]
        catacc[:,1] = rowcol[:,1]
        catacc[:,2] = acc_array[rowcol[:,0],rowcol[:,1]]
        catacc = catacc[catacc[:,2].argsort()]
    
        if nrow != catacc[0,0] or ncol != catacc[0,1]:
             nnrow,nncol = Nextcell(dirarc_array,nrow,ncol)
             print(nstrid)
             if nnrow <= 0 or nncol <=0 or nnrow >=nrows or nncol >= ncols:
                 continue
             nnstrid = strtemp_array[nnrow,nncol]
             strtemp_array[nrow,ncol] = nnstrid
        
    temparray = garray.array()
    temparray[:,:] = 0
    temparray[:,:] = strtemp_array[:,:]
    temparray.write(mapname="str_grass_rf", overwrite=True)
    
    grass.run_command('r.null', map='str_grass_rf',setnull=0)
    grass.run_command('r.mapcalc',expression = 'str_grass_rfn = int(str_grass_rf)',overwrite = True)
    grass.run_command('r.thin',input = 'str_grass_rfn', output = 'str_grass_r',overwrite = True)
    
    grass.run_command('r.stream.basins',direction = 'dir_Grass', stream = 'str_grass_r', basins = 'cat1',overwrite = True)
    
    grass.run_command('r.out.gdal', input = 'cat1',output = OutputFolder  + 'cat1.tif',format= 'GTiff',overwrite = True)
    grass.run_command('r.out.gdal', input = 'str_grass_r',output = OutputFolder  + 'str.tif',format= 'GTiff',overwrite = True)
    
    grass.run_command('r.to.vect',  input = 'str_grass_r',output = 'str', type ='line' ,overwrite = True)
################ check connected lakes  and non connected lakes 
#    grass.run_command('v.overlay',ainput = 'str_grass_v',binput = 'Hylake',operator = 'and',output = 'str_lake_and',overwrite = True)
#    grass.run_command('v.db.join',map = 'Hylake',column = 'Hylak_id',other_table = 'str_lake_and',other_column ='b_Hylak_id',overwrite = True)    
    grass.run_command('v.select',ainput = 'Hylake',binput = 'str_grass_v',output = 'lake_str',overwrite = True)

    grass.run_command('v.out.ogr', input = 'lake_str',output = os.path.join(OutputFolder, "Connect_lake.shp"),format= 'ESRI_Shapefile',overwrite = True)
    
    os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(cellSize) + "  " +str(cellSize) +"   " + "\"" +  os.path.join(OutputFolder, "Connect_lake.shp") +"\""+ "    "+ "\""+os.path.join(OutputFolder, "cnhylakegdal.tif")+"\"")
    grass.run_command("r.in.gdal", input = os.path.join(OutputFolder, "cnhylakegdal.tif"), output = 'cnlakeraster_in', overwrite = True)
    grass.run_command('r.mapcalc',expression = 'Connect_Lake = int(cnlakeraster_in)',overwrite = True)
#    grass.run_command('v.to.rast',input = 'lake_str',output = 'Connect_Lake',use = 'attr',attribute_column = 'Hylak_id',overwrite = True)
    
    grass.run_command('r.mapcalc',expression = 'Nonconnect_Lake = if(isnull(Connect_Lake),alllake,-9)',overwrite = True)
#    grass.run_command('v.overlay',ainput = 'str_grass_v',binput = 'Hylake',operator = 'or',output = 'str_lake',overwrite = True)    
    
    
    PERMANENT.close()
#    shutil.rmtree(gisdb,ignore_errors=True)
    
    
#    Watershed = Module("r.stream.basins")
#    Watershed(direction='dir.tif', points = 'obspoint.shp', basins = 'cat12s.tif')

#    
#    p.kill()
#    sys.exit(0)
#    print("asdfasdfsadfadsfdsaf")
    
    
    
        
