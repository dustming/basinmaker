from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def define_project_extent_using_outlet_point(
    grassdb,
    grass_location,
    qgis_prefix_path,
    path_dem_in,
    outlet_pt,
    buffer_distance=0,
    mask="MASK",
    dem="dem",
):

    """Define processing extent

    Function that used to define project processing spatial extent (PSE).
    The processing spatial extent is a region where Toolbox will work in. Toolbox
    will not process grids or features outside the processing spatial extent.
    Several options is available here. The PSE
    can be defined using DEM and an point coordinates. the drainage area
    contribute to that point coordinate will be used as boundary polygon.

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
    outlet_pt                        : list (required)
        It is list that indicate the outlet coordinates of the
        region of interest. If it is provided, the PSE
        will be defined as the drainage area controlled by this point.
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

    print("mask region:   using outlet point ")

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

    PERMANENT_Temp1 = Session()
    PERMANENT_Temp1.open(
        gisdb=grassdb,
        location=grass_location + "t1",
        create_opts="EPSG:4326",
    )

    # import clipped dem to target location
    grass_raster_r_in_gdal(
        grass=grass,
        raster_path=path_dem_in,
        output_nm=dem,
        location=grass_location + "t2",
    )

    PERMANENT_Temp1.close()

    # open location where dem is loaded, generated watershed boundary in this location
    PERMANENT_Temp = Session()
    PERMANENT_Temp.open(
        gisdb=grassdb,
        location=grass_location + "t2",
        create_opts="",
    )

    # define mask of current working enviroments
    grass_raster_r_mask(grass, dem)
    # define processing extent of the current working enviroment
    grass_raster_g_region(grass, dem)

    # generate watrshed boundary for this outlet point

    # define flow direction from dem
    grass_raster_r_watershed(
        grass,
        elevation=dem,
        drainage="dir_grass",
        accumulation="acc_grass2",
        flags="sa",
    )
    # generate watershed boundary
    grass_raster_r_water_outlet(
        grass,
        input_dir_nm="dir_grass",
        output_watshed_nm="wat_mask",
        outlet_coordinates=outlet_pt,
    )
    # define mask with watershed boundary
    grass_raster_r_mask(grass, "wat_mask")

    # export generated mask to folder ouside grass work env in tif format
    grass_raster_r_out_gdal(
        grass,
        input_nm="MASK",
        output=os.path.join(grassdb, mask + "_t1.tif"),
        format="GTiff",
    )

    # polygonize exported mask raster, buffer, and use it clip input dem
    qgis_raster_gdal_polygonize(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_t1.tif"),
        OUTPUT=os.path.join(grassdb, mask + "_t1.shp"),
    )
    qgis_vector_dissolve(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_t1.shp"),
        FIELD="DN",
        OUTPUT=os.path.join(grassdb, mask + "_t2.shp"),
    )
    qgis_vector_buffer(
        processing,
        context,
        INPUT=os.path.join(grassdb, mask + "_t2.shp"),
        Buffer_Distance=buffer_distance,
        OUTPUT=os.path.join(grassdb, mask + ".shp"),
    )
    # using saga function, becasue for some reason the gdal function did
    # not reduce the extent of the raster, the raster still has the same
    # size, but the value out side of mask is set to null
    qgis_raster_saga_clip_raster_with_polygon(
        processing,
        context,
        Input=path_dem_in,
        MASK=os.path.join(grassdb, mask + ".shp"),
        Output=os.path.join(grassdb, dem + ".sdat"),
    )

    # load clipted dem raster to target location
    grass_raster_r_in_gdal(
        grass=grass,
        raster_path=os.path.join(grassdb, dem + ".sdat"),
        output_nm=dem,
        location=grass_location,
    )

    PERMANENT_Temp.close()

    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    # define mask of current working enviroments
    grass_raster_r_mask(grass, dem)
    # define processing extent of the current working enviroment
    grass_raster_g_region(grass, dem)

    grass_raster_r_out_gdal(
        grass,
        input_nm=mask,
        output=os.path.join(grassdb, mask + ".tif"),
        format="GTiff",
    )
    PERMANENT.close()
    Qgs.exit()
