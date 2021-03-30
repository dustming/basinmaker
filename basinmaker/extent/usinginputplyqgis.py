from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def define_project_extent_using_input_polygon(
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_dem_in,
    path_extent_ply,
    mask="MASK",
    dem="dem",
):

    """Define processing extent

    Function that used to define project processing spatial extent (PSE).
    The processing spatial extent is a region where Toolbox will work in. Toolbox
    will not process grids or features outside the processing spatial extent.
    Several options is available here. The PSE can be defined using input polygon

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
    path_extent_ply                  : string (optional)
        It is the path of a subregion polygon. It is only used when the Region
        of interest is very large and the resolution of the dem is very high.
        toolbox will first divide the whole region into several small subregions.
        And then using devided subregion polygon as PSE.
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

    print("mask region:   using input polygon ")

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

    # load subregion mask polygon to target grass location
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=grassdb,
        location=grass_location + "t1",
        create_opts="EPSG:4326",
    )
    grass_raster_v_in_org(
        grass,
        input_path=path_extent_ply,
        output_vector_nm=mask,
        location=grass_location,
    )
    PERMANENT.close()

    # unpack dem and generate mask with loaded vector

    PERMANENT = Session()
    # grass unpack function report an error if gisdb is provided using relative
    # path, so absolute path is provided here
    PERMANENT.open(
        gisdb=os.path.abspath(grassdb), location=grass_location, create_opts=""
    )

    # unpack dem dataset
    root, ext = os.path.splitext(path_dem_in)
    if ext == ".pack":
        grass_raster_r_unpack(
            grass, input=os.path.abspath(path_dem_in), output=dem + "_big"
        )
    else:
        grass_raster_r_in_gdal(
            grass=grass, raster_path=path_dem_in, output_nm=dem + "_big"
        )
    ###Use polygon to define a smaller region
    # use big dem define region first and then change it by vector mask
    grass_raster_g_region(grass, dem + "_big")
    # define mask of current working enviroments
    grass_raster_r_mask(grass, raster_nm="#", vector_nm=mask)
    # define new region with this mask
    grass_raster_g_region(grass, raster_nm="#", zoom="MASK")

    # using vector define mask and clip big raster dem
    grass_raster_r_mask(grass, raster_nm="#", vector_nm=mask)
    exp = "%s = %s" % (dem, dem + "_big")
    grass_raster_r_mapcalc(grass, expression=exp)

    # using clipped dem to defing mask
    grass_raster_r_mask(grass, raster_nm=dem)

    # export dem and mask and generate mask polygon
    grass_raster_r_out_gdal(
        grass,
        input_nm=mask,
        output=os.path.join(grassdb, mask + ".tif"),
        format="GTiff",
    )
    grass_raster_r_out_gdal(
        grass,
        input_nm=dem,
        output=os.path.join(grassdb, dem + ".tif"),
        format="GTiff",
    )
    # polygonize exported mask raster, and polygonize and dissove it.
    qgis_raster_gdal_polygonize(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + ".tif"),
        OUTPUT=os.path.join(grassdb, mask + "1.shp"),
    )
    # qgis dissolve function for some reason not working here....
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "1.shp"),
        FIELD="DN",
        OUTPUT=os.path.join(grassdb, mask + ".shp"),
        USING_GDAL_FUNCTION=True,
    )

    PERMANENT.close()

    Qgs.exit()
