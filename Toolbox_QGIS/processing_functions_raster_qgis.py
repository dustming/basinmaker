from qgis.core import *
import qgis
from qgis.analysis import QgsNativeAlgorithms
from qgis.PyQt.QtCore import *
import pandas as pd 
import numpy as np 


def qgis_raster_gdal_warpreproject(processing,Input,TARGET_CRS, Output):
    """ Functions reproject raster layer 
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """

    out = processing.run("gdal:warpreproject", {'INPUT':Input,'SOURCE_CRS':None,'TARGET_CRS':QgsCoordinateReferenceSystem(TARGET_CRS),
                                                'RESAMPLING':0,'NODATA':None,'TARGET_RESOLUTION':None,'OPTIONS':'','DATA_TYPE':0,
                                                'TARGET_EXTENT':None,'TARGET_EXTENT_CRS':None,'MULTITHREADING':False,
                                                'EXTRA':'','OUTPUT':Output})
    return out 
############## 



def qgis_raster_clip_raster_by_mask(processing,Input,MASK,TARGET_CRS, Output):
    """ Functions clip raster by mask
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """

    out = processing.run("gdal:cliprasterbymasklayer", {'INPUT':Input,'MASK':MASK,'SOURCE_CRS':None,
                               'TARGET_CRS':QgsCoordinateReferenceSystem(TARGET_CRS),'NODATA':None,'ALPHA_BAND':False,
                                'CROP_TO_CUTLINE':True,'KEEP_RESOLUTION':False,
                                'SET_RESOLUTION':False,'X_RESOLUTION':None,
                                 'Y_RESOLUTION':None,'MULTITHREADING':False,'OPTIONS':'',
                                 'DATA_TYPE':0,'EXTRA':'',
                                'OUTPUT':Output})
    return out
    

def qgis_raster_slope(processing,Input, Output):
    """ Functions calculate slope from DEM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """

    out = processing.run("qgis:slope", {'INPUT':Input,'Z_FACTOR':1,'OUTPUT':Output})
    return out
    
def qgis_raster_aspect(processing,Input, Output):
    """ Functions calculate aspect from DEM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """

    out = processing.run("qgis:aspect", {'INPUT':Input,'Z_FACTOR':1,'OUTPUT':Output})
    return out


def qgis_raster_zonal_statistics(processing,INPUT_RASTER,INPUT_VECTOR,COLUMN_PREFIX):
    """ Functions calculate zonal statistics between polygon and raster 
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated 
    """

    processing.run("qgis:zonalstatistics", {'INPUT_RASTER':INPUT_RASTER,'RASTER_BAND':1,'INPUT_VECTOR':INPUT_VECTOR,'COLUMN_PREFIX':COLUMN_PREFIX,'STATS':[2]})
    return
    
def qgis_raster_read_raster(processing,INPUT):
    """ qgis function read raster into memeory
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """  
    out = QgsRasterLayer(INPUT, "")
    return out 


def qgis_raster_return_raster_properties(processing,INPUT):
    """ qgis function return raster projection and cellsize
    ----------

    Notes
    -------

    Returns:
    -------
        None, 
    """  
    cellSize = float(INPUT.rasterUnitsPerPixelX())  ### Get Raster cell size
    SpRef_in = INPUT.crs().authid()   ### get Raster spatialReference id
    return cellSize,SpRef_in

        