import os

import numpy as np
import pandas as pd
import qgis
from qgis.analysis import QgsNativeAlgorithms
from qgis.core import *
from qgis.PyQt.QtCore import *


def qgis_raster_gdal_warpreproject(processing, Input, TARGET_CRS, Output):
    """Functions reproject raster layer
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "gdal:warpreproject",
        {
            "INPUT": Input,
            "SOURCE_CRS": None,
            "TARGET_CRS": QgsCoordinateReferenceSystem(TARGET_CRS),
            "RESAMPLING": 0,
            "NODATA": None,
            "TARGET_RESOLUTION": None,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "TARGET_EXTENT": None,
            "TARGET_EXTENT_CRS": None,
            "MULTITHREADING": False,
            "EXTRA": "",
            "OUTPUT": Output,
        },
    )
    return out


##############


def qgis_raster_clip_raster_by_mask(processing, Input, MASK, TARGET_CRS, Output):
    """Functions clip raster by mask
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "gdal:cliprasterbymasklayer",
        {
            "INPUT": Input,
            "MASK": MASK,
            "SOURCE_CRS": None,
            "TARGET_CRS": QgsCoordinateReferenceSystem(TARGET_CRS),
            "NODATA": None,
            "ALPHA_BAND": False,
            "CROP_TO_CUTLINE": True,
            "KEEP_RESOLUTION": False,
            "SET_RESOLUTION": False,
            "X_RESOLUTION": None,
            "Y_RESOLUTION": None,
            "MULTITHREADING": False,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "EXTRA": "",
            "OUTPUT": Output,
        },
    )
    return out


def qgis_raster_clip_raster_by_extent(processing, Input, PROJWIN, Output):
    """Functions clip raster by extent
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """
    out = processing.run(
        "gdal:cliprasterbyextent",
        {
            "INPUT": Input,
            "PROJWIN": PROJWIN,
            "NODATA": None,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "EXTRA": "",
            "OUTPUT": Output,
        },
    )

    return out


def qgis_raster_saga_clip_raster_with_polygon(processing, context, Input, MASK, Output):
    """Functions clip raster by mask,
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "saga:cliprasterwithpolygon",
        {"INPUT": Input, "POLYGONS": MASK, "OUTPUT": Output},
        context=context,
    )

    return out


def qgis_raster_slope(processing, Input, Output):
    """Functions calculate slope from DEM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "qgis:slope", {"INPUT": Input, "Z_FACTOR": 1, "OUTPUT": Output}
    )
    return out


def qgis_raster_aspect(processing, Input, Output):
    """Functions calculate aspect from DEM
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    out = processing.run(
        "qgis:aspect", {"INPUT": Input, "Z_FACTOR": 1, "OUTPUT": Output}
    )
    return out


def qgis_raster_zonal_statistics(processing, INPUT_RASTER, INPUT_VECTOR, COLUMN_PREFIX):
    """Functions calculate zonal statistics between polygon and raster
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    processing.run(
        "qgis:zonalstatistics",
        {
            "INPUT_RASTER": INPUT_RASTER,
            "RASTER_BAND": 1,
            "INPUT_VECTOR": INPUT_VECTOR,
            "COLUMN_PREFIX": COLUMN_PREFIX,
            "STATS": [2],
        },
    )
    return


def qgis_raster_read_raster(processing, INPUT):
    """qgis function read raster into memeory
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = QgsRasterLayer(INPUT, "")
    return out


def qgis_raster_gdal_translate(processing, INPUT, OUTPUT, format="GTiff"):
    """qgis function read raster into memeory
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    params = {"INPUT": INPUT, "format": format, "OUTPUT": OUTPUT}
    out = processing.run("gdal:translate", params)
    return out


def qgis_raster_return_raster_properties(processing, INPUT):
    """qgis function return raster projection and cellsize
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    cellSize = float(INPUT.rasterUnitsPerPixelX())  ### Get Raster cell size
    SpRef_in = INPUT.crs().authid()  ### get Raster spatialReference id
    return cellSize, SpRef_in


def qgis_raster_gdal_polygonize(processing, context, INPUT, OUTPUT):
    """qgis function polygonize rater to polygon
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    out = processing.run(
        "gdal:polygonize",
        {
            "INPUT": INPUT,
            "BAND": 1,
            "FIELD": "DN",
            "EIGHT_CONNECTEDNESS": False,
            "EXTRA": "",
            "OUTPUT": OUTPUT,
        },
        context=context,
    )

    return out


def qgis_raster_gdal_rasterize(
    processing, context, INPUT, Column_nm, cellsize, w, s, e, n, OUTPUT
):
    """qgis gdal function polygonize rater to polygon
    ----------

    Notes
    -------

    Returns:
    -------
        None,
    """
    cmd = (
        "gdal_rasterize -at -of GTiff -a_nodata -9999 -a "
        + Column_nm
        + " -tr  "
        + str(cellsize)
        + "  "
        + str(cellsize)
        + "  -te   "
        + str(w)
        + "   "
        + str(s)
        + "   "
        + str(e)
        + "   "
        + str(n)
        + "   "
        + '"'
        + INPUT
        + '"'
        + "    "
        + '"'
        + OUTPUT
        + '"'
    )
    os.system(cmd)
    return


# os.system('gdal_rasterize -at -of GTiff -a_nodata -9999 -a Hylak_id -tr  '+ str(self.cellSize) + "  " +str(self.cellSize)+'  -te   '+ str(grsregion['w'])+"   " +str(grsregion['s'])+"   " +str(grsregion['e'])+"   " +str(grsregion['n'])+"   " + "\"" +  self.Path_allLakeply +"\""+ "    "+ "\""+self.Path_allLakeRas+"\"")
