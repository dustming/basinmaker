from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def define_project_extent_using_hybasin_ply(
    work_folder,
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

    # read and define arcgis work enviroments 
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    cellSize = float(arcpy.GetRasterProperties_management(path_dem_in, "CELLSIZEX").getOutput(0))
    SptailRef = arcpy.Describe(path_dem_in).spatialReference
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize    
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    if not os.path.exists(work_folder):
        os.makedirs(work_folder)
    arcpy.env.workspace = work_folder


    hyshdinfo = Dbf_To_Dataframe(hybasin_ply)
    hyshdinfo = pd.DataFrame.spatial.from_featureclass(hybasin_ply)
    
    routing_info = hyshdinfo[["HYBAS_ID", "NEXT_DOWN"]].astype("float").values

    # obtain sub id of subbasins between OutHyID and OutHyID2 in the routing
    # network
    HydroBasins = Return_SubIds_Between_Two_Subbasins_In_Rouing_Network(
        routing_info, down_hybasin_id, up_hybasin_id
    )

    hyshdinfo_select = hyshdinfo.loc[hyshdinfo['HYBAS_ID'].isin(HydroBasins)]
    hyshdinfo_select['MASK'] = 1
    
    hyshdinfo_select.spatial.to_featureclass(location=os.path.join(work_folder,mask + "_hy.shp"),overwrite=True,sanitize_columns=False)
    
    arcpy.Dissolve_management(mask + "_hy.shp", mask + "_1hy.shp", ['MASK'])
    
    arcpy.Buffer_analysis(mask + "_1hy.shp",  mask + "_2hy.shp", buffer_distance)
    
    arcpy.Project_management(
        mask + "_2hy.shp",
        mask + ".shp", 
        arcpy.SpatialReference(int(SptailRef.factoryCode)),
        )    

    dem_mask_raster = ExtractByMask(path_dem_in, mask + ".shp")
    dem_mask_raster.save(dem+'.tif')
    dem_mask_raster.save(mask+'.tif')

    return 
