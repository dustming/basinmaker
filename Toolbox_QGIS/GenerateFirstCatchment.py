def WatershedDiscretizationToolset(OutputFolder,Thresholdacc,GRASS_BIN):
    import os
    import sys
    import tempfile
    import shutil
    
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

    grass.run_command('v.to.rast',input = 'str_grass_v',output = 'str_grass_r',use = 'cat',overwrite = True)
    
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
    
    
    
        
