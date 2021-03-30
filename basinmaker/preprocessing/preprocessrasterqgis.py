from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import os


def preprocess_raster(
    grassdb,
    grass_location,
    qgis_prefix_path,
    mask,
    raster_path,
    raster_name,
):

    QgsApplication.setPrefixPath(qgis_prefix_path, True)
    Qgs = QgsApplication([], False)
    Qgs.initQgis()
    from processing.core.Processing import Processing
    from processing.tools import dataobjects
    from qgis import processing

    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    mask_layer = qgis_raster_read_raster(
        processing, os.path.join(grassdb, mask + ".tif")
    )  ### load DEM raster as a  QGIS raster object to obtain attribute

    raster_layer = qgis_raster_read_raster(processing, raster_path)

    cellSize, SpRef_in = qgis_raster_return_raster_properties(
        processing, mask_layer
    )  ### Get Raster cell size

    cellSize, SpRef_src = qgis_raster_return_raster_properties(
        processing, raster_layer
    )  ### Get Raster cell size

    processing.run(
        "gdal:warpreproject",
        {
            "INPUT": raster_path,
            "SOURCE_CRS": QgsCoordinateReferenceSystem(SpRef_src),
            "TARGET_CRS": QgsCoordinateReferenceSystem(SpRef_in),
            "RESAMPLING": 0,
            "NODATA": None,
            "TARGET_RESOLUTION": None,
            "OPTIONS": "",
            "DATA_TYPE": 0,
            "TARGET_EXTENT": None,
            "TARGET_EXTENT_CRS": None,
            "MULTITHREADING": False,
            "EXTRA": "",
            "OUTPUT": os.path.join(grassdb, raster_name + "_proj.tif"),
        },
    )

    processing.run(
        "gdal:cliprasterbymasklayer",
        {
            "INPUT": os.path.join(grassdb, raster_name + "_proj.tif"),
            "MASK": os.path.join(grassdb, mask + ".shp"),
            "SOURCE_CRS": None,
            "TARGET_CRS": None,
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
            "OUTPUT": os.path.join(grassdb, raster_name + ".tif"),
        },
    )

    Qgs.exit()

    return
