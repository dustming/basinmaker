
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



def Generateinputdata_hydrosheds(hyshddem,hyshddir,hyshdacc,hyshdply,WidDep,Lakefile,Landuse,Landuseinfo,obspoint,OutHyID,OutHyID2,OutputFolder):

    from qgis.core import (
        QgsRasterLayer,
        QgsVectorLayer,
        QgsApplication,
        QgsFeatureRequest,
        QgsVectorFileWriter
    )
    from simpledbf import Dbf5
    import os
    import sys
    import numpy as np
    
#    sys.path.append('C:/Users/dustm/Documents/GitHub/RoutingTool/Toolbox_QGIS/')
    
    OutHyID = int(OutHyID)
    OutHyID2 = int(OutHyID2)
    
    QgsApplication.setPrefixPath("C:/QGIS310/apps/qgis", True)
    Qgs = QgsApplication([],False)
    Qgs.initQgis()
    
    r_dem_layer = QgsRasterLayer(hyshddem, "") ### load DEM raster as a  QGIS raster object to obtain attribute
    if not r_dem_layer.isValid():
        print("Layer failed to load!")
            
    cellSize = float(r_dem_layer.rasterUnitsPerPixelX())  ### Get Raster cell size

    SptailRefid = r_dem_layer.crs().authid()   ### get Raster spatialReference id
    
    print("Working with a  sptail reference  :   " , r_dem_layer.crs().description(), "      ", SptailRefid)
    print("The cell cize is   ",cellSize)
    
    if not os.path.exists(OutputFolder):
        os.makedirs(OutputFolder)

        
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
        selectedFeatureID.append(feature.id())
        
    hyshedl12.select(selectedFeatureID)   ### select with polygon id
    
    # Save selected polygons to output 
    _writer = QgsVectorFileWriter.writeAsVectorFormat(hyshedl12, OutputFolder +"HyMask.shp", "UTF-8", hyshedl12.crs(), "ESRI Shapefile", onlySelected=True)

#arcpy.AddMessage(HydroBasins)

    Qgs.exitQgis()
    