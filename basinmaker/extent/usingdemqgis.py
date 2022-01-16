from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def define_project_extent_using_dem(
    grassdb, grass_location, qgis_prefix_path, path_dem_in, mask="MASK", dem="dem"
):

    """Define processing extent

    Function that used to define project processing spatial extent (PSE).
    The processing spatial extent is a region where Toolbox will work in. Toolbox
    will not process grids or features outside the processing spatial extent.
    Several options is available here. The PSE can be defined as the extent of
    input DEM.

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

    print("Mask Region:   Using provided DEM : ")

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
        dict(
            GRASS_COMPRESS_NULLS="1",
            GRASS_COMPRESSOR="ZSTD",
            GRASS_VERBOSE="1",
        )
    )
    PERMANENT = Session()
    PERMANENT.open(
        gisdb=grassdb,
        location=grass_location + "t1",
        create_opts="EPSG:4326",
    )

    # load dem raster to target location
    grass_raster_r_in_gdal(
        grass=grass,
        raster_path=path_dem_in,
        output_nm=dem,
        location=grass_location,
    )
    PERMANENT.close()

    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    exp = "%s = %s/%s" % (
        dem,
        dem,
        str(1),
    )        
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    
    # define mask of current working enviroments
    grass_raster_r_mask(grass, dem)
    # define processing extent of the current working enviroment
    grass_raster_g_region(grass, dem)

    # export generated mask to folder ouside grass work env in
    # tif format
    grass_raster_r_out_gdal(
        grass,
        input_nm="MASK",
        output=os.path.join(grassdb, mask + ".tif"),
        format="GTiff",
    )
    # polygonize exported mask raster, and polygonize and dissove it.
    qgis_raster_gdal_polygonize(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + ".tif"),
        OUTPUT=os.path.join(grassdb, mask + "_t1.shp"),
    )
    # qgis dissolve function for some reason not working here....
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_t1.shp"),
        FIELD="DN",
        OUTPUT=os.path.join(grassdb, mask + ".shp"),
        USING_GDAL_FUNCTION=True,
    )
    PERMANENT.close()

    Qgs.exit()

    return
