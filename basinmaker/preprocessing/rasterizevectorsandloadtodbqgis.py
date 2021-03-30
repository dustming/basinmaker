from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import os


def rasterize_vectors_and_load_to_db(
    grassdb,
    grass_location,
    qgis_prefix_path,
    mask,
    vector_path,
    attribue_name,
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
    cellSize, SpRef_in = qgis_raster_return_raster_properties(
        processing, mask_layer
    )  ### Get Raster cell size

    # load grass working location
    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    # get dem array and get nrows and ncols of the domain
    strtemp_array = Return_Raster_As_Array_With_garray(garray, mask)
    ncols = int(strtemp_array.shape[1])
    nrows = int(strtemp_array.shape[0])
    grsregion = gcore.region()

    qgis_raster_gdal_rasterize(
        processing,
        context,
        INPUT=vector_path,
        Column_nm=attribue_name,
        cellsize=cellSize,
        w=grsregion["w"],
        s=grsregion["s"],
        e=grsregion["e"],
        n=grsregion["n"],
        OUTPUT=os.path.join(grassdb, raster_name + ".tif"),
    )

    grass_raster_r_in_gdal(
        grass,
        raster_path=os.path.join(grassdb, raster_name + ".tif"),
        output_nm=raster_name,
    )

    grass_raster_setnull(
        grass, raster_nm=raster_name, null_values=[-9999], create_new_raster=False
    )

    grass_raster_v_import(grass, input_path=vector_path, output_vector_nm=raster_name)
    Qgs.exit()
    PERMANENT.close()


def obtain_lake_vectors_larger_than_threstholds(
    grassdb,
    qgis_prefix_path,
    all_lakes,
    lake_attributes,
    threshold_con_lake,
    threshold_non_con_lake,
    lakes_lg_cl_thres="lakes_lg_cl_thres",
    lakes_lg_ncl_thres="lakes_lg_ncl_thres",
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

    print(threshold_con_lake, threshold_non_con_lake, lake_attributes[2])
    qgis_vector_extract_by_attribute(
        processing,
        context,
        INPUT_Layer=os.path.join(grassdb, all_lakes + ".shp"),
        FIELD=lake_attributes[2],
        OPERATOR=2,
        VALUE=threshold_con_lake,
        OUTPUT=os.path.join(grassdb, lakes_lg_cl_thres + ".shp"),
    )

    qgis_vector_extract_by_attribute(
        processing,
        context,
        INPUT_Layer=os.path.join(grassdb, all_lakes + ".shp"),
        FIELD=lake_attributes[2],
        OPERATOR=2,
        VALUE=threshold_non_con_lake,
        OUTPUT=os.path.join(grassdb, lakes_lg_ncl_thres + ".shp"),
    )
    Qgs.exit()
