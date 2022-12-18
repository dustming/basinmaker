from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def define_project_extent_using_dem(
    work_folder,
    grass_location,
    qgis_prefix_path,
    path_dem_in,
    mask="MASK",
    dem="dem",
    path_to_snap_raster = "#",
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

    # read and define arcgis work enviroments
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    cellSize = float(arcpy.GetRasterProperties_management(path_dem_in, "CELLSIZEX").getOutput(0))
    extent = arcpy.Describe(path_dem_in).extent
    SptailRef = arcpy.Describe(path_dem_in).spatialReference
    arcgis_env = {
        "extent":extent,
        "georef":arcpy.SpatialReference(int(SptailRef.factoryCode)),
        "cellSize":cellSize,
        }
    if not os.path.exists(work_folder):
        os.makedirs(work_folder)
    arcpy.CreateFileGDB_management(work_folder, "arcgis.gdb")

    
    arcpy.env.XYTolerance = arcgis_env["cellSize"]
    arcpy.arcpy.env.cellSize = arcgis_env["cellSize"]
    arcpy.env.outputCoordinateSystem = arcgis_env["georef"] ### WGS84
    arcpy.env.extent = arcgis_env["extent"]
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")

    if path_to_snap_raster != "#":
        arcpy.env.snapRaster =  path_to_snap_raster

    # copy dem to the work fokder
    arcpy.CopyRaster_management(path_dem_in,dem )
    arcpy.env.extent = arcpy.Describe(dem).extent

    # create mask raster
    mask_raster = Con(dem, 1, 1, "VALUE >= 0")
    mask_raster.save(mask)
    # create mask polygon
    arcpy.RasterToPolygon_conversion(mask_raster, mask+'_ply', "NO_SIMPLIFY", "VALUE")

    return
