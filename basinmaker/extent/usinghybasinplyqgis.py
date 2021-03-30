from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def define_project_extent_using_hybasin_ply(
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_dem_in,
    buffer_distance,
    hybasin_ply,
    down_hybasin_id,
    up_hybasin_id,
    mask="MASK",
    dem="dem",
):

    """Define processing extent

    Function that used to define project processing spatial extent (PSE).
    The processing spatial extent is a region where Toolbox will work in. Toolbox
    will not process grids or features outside the processing spatial extent.
    Several options is available here. The PSE can be defined using Hybasin
    product and a hydrobasin ID. All subbasin drainage to that hydrobasin ID
    will be extracted. And the extent of the extracted polygon will be used as PSE

    Parameters
    ----------
    grassdb                           : path (required)
        It is a path to project grass database folder
    grass_location                    : string (required)
        It is a string of grass location name
    qgis_prefix_path                  : string (required)
        It is a string of qgis prefix path
    path_dem_in                      : string (required)
        It is the path to input dem
    buffer_distance                  : float (optional)
        It is a float number to increase the extent of the PSE
        obtained from Hydrobasins. It is needed when input DEM is not from
        HydroSHEDS. Then the extent of the watershed will be different
        with PSE defined by HydroBASINS.
    hybasin_ply                      : string (optional)
        It is a path to hydrobasin routing product, If it is provided, the
        PSE will be based on the OutHyID and OutHyID2 and
        this HydroBASINS routing product.
    down_hybasin_id                  : int (optional)
        It is a HydroBASINS subbasin ID, which should be the ID of the most
        downstream subbasin in the region of interest.
    up_hybasin_id                    : int (optional)
        It is a HydroBASINS subbasin ID, which should be the ID of the most
        upstream subbasin in the region of interest, normally do not needed.
    mask                             : string (optional)
        It is a output mask name, which will stored in grass_location in both
        vector and raster format
    dem                              : string (optional)
        It is a output dem raster name, which will be stored in grass_location

    Notes
    -------
    Outputs are following files

    MASK                   : raster
        it is a mask raster stored in grass database, which indicate
        the PSE. The grass database is located at
        os.path.join(grassdb, grass_location)
    dem                   : raster
        it is a dem raster stored in grass database, which is
        has the same extent with MASK. The grass database is located at
        os.path.join(grassdb, grass_location)

    Returns:
    -------
       None

    Examples
    -------
    """

    print("mask region:   using hybasin polygon ")

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

    r_dem_layer = qgis_raster_read_raster(
        processing, path_dem_in
    )  ### load DEM raster as a  QGIS raster object to obtain attribute
    cellSize, SpRef_in = qgis_raster_return_raster_properties(
        processing, r_dem_layer
    )  ### Get Raster cell size

    hyshdinfo = Dbf_To_Dataframe(hybasin_ply)
    routing_info = hyshdinfo[["HYBAS_ID", "NEXT_DOWN"]].astype("float").values

    # obtain sub id of subbasins between OutHyID and OutHyID2 in the routing
    # network
    HydroBasins = Return_SubIds_Between_Two_Subbasins_In_Rouing_Network(
        routing_info, down_hybasin_id, up_hybasin_id
    )

    # extract subbasins from hydrobasin product
    Selectfeatureattributes(
        processing,
        Input=hybasin_ply,
        Output=os.path.join(grassdb, mask + "_hy.shp"),
        Attri_NM="HYBAS_ID",
        Values=HydroBasins,
    )

    print("Mask Region:   Using buffered hydroBasin product polygons ")

    # dissolve, buffer and reproject the extracted hydrobasin product
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_hy.shp"),
        FIELD="MAIN_BAS",
        OUTPUT=os.path.join(grassdb, mask + "_1hy.shp"),
    )
    qgis_vector_buffer(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_1hy.shp"),
        Buffer_Distance=buffer_distance,
        OUTPUT=os.path.join(grassdb, mask + "_2hy.shp"),
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_2hy.shp"),
        FIELD="MAIN_BAS",
        OUTPUT=os.path.join(grassdb, mask + "_3hy.shp"),
    )
    qgis_vector_reproject_layers(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_3hy.shp"),
        TARGET_CRS=SpRef_in,
        OUTPUT=os.path.join(grassdb, mask + ".shp"),
    )

    # clip raster layer with this mask
    qgis_raster_clip_raster_by_mask(
        processing,
        Input=path_dem_in,
        MASK=os.path.join(grassdb, mask + ".shp"),
        TARGET_CRS=SpRef_in,
        Output=os.path.join(grassdb, mask + ".tif"),
    )

    # use clipped DEM to great a grass work enviroment
    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    # open/create a grass location
    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT_temp = Session()
    PERMANENT_temp.open(
        gisdb=grassdb,
        location=grass_location + "t1",
        create_opts="EPSG:4326",
    )

    # import clipped dem to target location
    grass_raster_r_in_gdal(
        grass=grass,
        raster_path=os.path.join(grassdb, mask + ".tif"),
        output_nm=dem,
        location=grass_location,
    )
    PERMANENT_temp.close()

    # Define mask and processing region for grass working enviroments
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    # define mask of current working enviroments
    grass_raster_r_mask(grass, dem)
    # define processing extent of the current working enviroment
    grass_raster_g_region(grass, dem)

    PERMANENT.close()
    Qgs.exit()
