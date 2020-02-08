
###################################################################3
def Defcat(out,outletid):
    import numpy as np
    import copy
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


def Generatmaskregion(hyshdply,OutHyID,OutHyID2,OutputFolder):
    from simpledbf import Dbf5
    from qgis.core import (
        QgsRasterLayer,
        QgsVectorLayer,
        QgsApplication,
        QgsFeatureRequest,
        QgsVectorFileWriter,
        QgsProcessingFeedback
    )
    hyinfocsv = hyshdply[:-3] + "dbf"
    VolThreshold = 0
    tempinfo = Dbf5(hyinfocsv)#np.genfromtxt(hyinfocsv,delimiter=',')
    hyshdinfo = tempinfo.to_dataframe().values

### obtain the HydroShed polygons that belongs to the target drainage area
    HydroBasins1 = Defcat(hyshdinfo,OutHyID) ### return fid of polygons that needs to be select 

### 
    if OutHyID2 > 0:
        HydroBasins2 = Defcat(hyshdinfo,OutHyID2)
        
###  exculde the Ids in HydroBasins2 from HydroBasins1
        for i in range(len(HydroBasins2)):
            rows =np.argwhere(HydroBasins1 == HydroBasins2[i])
            HydroBasins1 = np.delete(HydroBasins1, rows)
        
        HydroBasins = HydroBasins1
    else:
        HydroBasins = HydroBasins1
    
### Load HydroSHED Layers 
    hyshedl12 = QgsVectorLayer(hyshdply, "")
    
### Build qgis selection expression
    where_clause = '"HYBAS_ID" IN'+ " ("
    for i in range(0,len(HydroBasins)):
        if i == 0:
            where_clause = where_clause + str(HydroBasins[i])
        else:
            where_clause = where_clause + "," + str(HydroBasins[i])
    where_clause = where_clause + ")"

    req = QgsFeatureRequest().setFlags( QgsFeatureRequest.NoGeometry )
    req.setFilterExpression(where_clause)
    it = hyshedl12.getFeatures( req )
    
    ### obtain all feature id of selected polygons
    selectedFeatureID = []

    for feature in it:
#        print(feature.id())
        selectedFeatureID.append(feature.id())
        
    hyshedl12.select(selectedFeatureID)   ### select with polygon id
    
    # Save selected polygons to output 
    _writer = QgsVectorFileWriter.writeAsVectorFormat(hyshedl12, OutputFolder +"HyMask.shp", "UTF-8", hyshedl12.crs(), "ESRI Shapefile", onlySelected=True)


def Generateinputdata_hydrosheds(dem_in = '#',dir_in = '#',hyshdply = '#',WidDep = '#',Lakefile = '#'
                                 ,Landuse = '#',Landuseinfo = '#',obspoint = '#',OutHyID = -1 ,OutHyID2 = -1, 
                                 OutputFolder = '#'):

    from qgis.core import (
        QgsRasterLayer,
        QgsVectorLayer,
        QgsApplication,
        QgsFeatureRequest,
        QgsVectorFileWriter,
        QgsProcessingFeedback
    )
    import qgis
    from qgis.analysis import QgsNativeAlgorithms
    from simpledbf import Dbf5
    import os
    import sys
    import numpy as np
    import shutil
    from shutil import copyfile
    import tempfile
#    sys.path.append('C:/Users/dustm/Documents/GitHub/RoutingTool/Toolbox_QGIS/')
        
    #### set up qgis applicaiton 
    qgisPP = os.environ['QGISPrefixPath']
    RoutingToolPath = os.environ['RoutingToolFolder']
    
    QgsApplication.setPrefixPath(qgisPP, True)
    Qgs = QgsApplication([],False)
    Qgs.initQgis()
    from qgis import processing
    from processing.core.Processing import Processing   
    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
#    Processing.updateAlgsList()
    #######################################3
    
    r_dem_layer = QgsRasterLayer(dem_in, "") ### load DEM raster as a  QGIS raster object to obtain attribute        
    cellSize = float(r_dem_layer.rasterUnitsPerPixelX())  ### Get Raster cell size
    SptailRefid = r_dem_layer.crs().authid()   ### get Raster spatialReference id
    print("Working with a  sptail reference  :   " , r_dem_layer.crs().description(), "      ", SptailRefid)
    print("The cell cize is   ",cellSize)
    
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

    if  hyshdply != '#':   #### input is using hydroshed DEM and hydroshed polygons 
         Generatmaskregion(hyshdply,OutHyID,OutHyID2,OutputFolder) ### creat OutputFolder +"HyMask.shp"
         processing.run('gdal:dissolve', {'INPUT':OutputFolder +"HyMask.shp",'FIELD':'MAIN_BAS','OUTPUT':OutputFolder +"HyMask2.shp"})
         params = {'INPUT': dem_in,'MASK': OutputFolder +"HyMask2.shp",'NODATA': -9999,'ALPHA_BAND': False,'CROP_TO_CUTLINE': True,
                                                                    'KEEP_RESOLUTION': True,
                                                                    'OPTIONS': 'COMPRESS=LZW',
                                                                    'DATA_TYPE': 0,  # Byte
                                                                    'OUTPUT': OutputFolder +"dem.tif"}
         dem = processing.run('gdal:cliprasterbymasklayer',params)  #### extract dem
    
    else:
         params = {'INPUT': dem_in, 'format': 'GTiff', 'OUTPUT': OutputFolder +"dem.tif"}
         processing.run('gdal:translate',params)
#### clip vector fiels  
    processing.run("native:clip", {'INPUT':Lakefile,'OVERLAY':OutputFolder +"HyMask2.shp",'OUTPUT':OutputFolder +"Hylake.shp"})
    processing.run("native:clip", {'INPUT':WidDep,'OVERLAY':OutputFolder +"HyMask2.shp",'OUTPUT':OutputFolder +"WidDep.shp"})
    processing.run("native:clip", {'INPUT':obspoint,'OVERLAY':OutputFolder +"HyMask2.shp",'OUTPUT':OutputFolder + "obspoint.shp"})
    copyfile( OutputFolder + "/"+"HyMask.prj" ,  OutputFolder + "/"+"obspoint.prj")
    
###### set up GRASS environment for translate vector to rasters and clip rasters
    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.pygrass.modules import Module
    from grass_session import Session
    
    gisdb =os.path.join(tempfile.gettempdir(), 'grassdata_toolbox')# "C:/Users/dustm/Documents/ubuntu/share/OneDrive/OneDrive - University of Waterloo/Documents/RoutingTool/Samples/Examples/"

    shutil.rmtree(gisdb,ignore_errors=True)

    try:
        os.stat(gisdb)
    except:
        os.mkdir(gisdb)
        
    os.environ['GISDBASE'] = gisdb    
    os.environ.update(dict(GRASS_COMPRESS_NULLS='1',GRASS_COMPRESSOR='ZSTD'))
    PERMANENT = Session()
    PERMANENT.open(gisdb=gisdb, location='mytest',create_opts=SptailRefid)
    
    grass.run_command("r.import", input = OutputFolder +"dem.tif", output = 'dem', overwrite = True)
    grass.run_command('g.region', raster='dem')
    grass.run_command('r.mask'  , raster='dem', maskcats = '*',overwrite = True)
    
    if dir_in != '#':
        grass.run_command("r.external", input = dir_in, output = 'dir_in', overwrite = True)
        grass.run_command("r.clip", input = 'dir_in', output = 'dir_Arcgis', overwrite = True, flags = 'r')
        grass.run_command('r.reclass', input='dir_Arcgis',output = 'dir_Grass',rules = RoutingToolPath + '/'+ 'Arcgis2GrassDIR.txt',overwrite = True)
    else:
        grass.run_command('r.watershed',elevation = 'dem', drainage = 'dir_Grass', overwrite = True)
        grass.run_command('r.reclass', input='dir_Grass',output = 'dir_Arcgis',rules = RoutingToolPath + '/'+ 'Grass2ArcgisDIR.txt',overwrite = True)
        
    grass.run_command("r.external", input = Landuse, output = 'landuse_in', overwrite = True)
    grass.run_command("r.clip", input = 'landuse_in', output = 'landuse', overwrite = True)        
    
    grass.run_command("v.import", input = OutputFolder +"WidDep.shp", output = 'WidDep', overwrite = True)
    grass.run_command("v.import", input = OutputFolder +"obspoint.shp", output = 'obspoint', overwrite = True)
    grass.run_command("v.import", input = OutputFolder +"Hylake.shp", output = 'Hylake', overwrite = True)

    grass.run_command('v.to.rast',input = 'WidDep',output = 'width',use = 'attr',attribute_column = 'WIDTH',overwrite = True)
    grass.run_command('v.to.rast',input = 'WidDep',output = 'depth',use = 'attr',attribute_column = 'DEPTH',overwrite = True)
    grass.run_command('v.to.rast',input = 'WidDep',output = 'qmean',use = 'attr',attribute_column = 'Q_Mean2',overwrite = True)
    grass.run_command('v.to.rast',input = 'obspoint',output = 'obs',use = 'attr',attribute_column = 'Obs_ID',overwrite = True)
    grass.run_command('v.to.rast',input = 'Hylake',output = 'alllake',use = 'attr',attribute_column = 'Hylak_id',overwrite = True)

    PERMANENT.close()
    Qgs.exit()