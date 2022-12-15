from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
from basinmaker.addlakeandobs.modifyfdr import modify_lakes_flow_direction

def add_lakes_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_lakefile_in,
    lake_attributes,
    threshold_con_lake,
    threshold_non_con_lake,
    only_included_lake_at_river_interction = False,
    remove_lake_inlets = False,
    path_sub_reg_lake_r="#",
    path_sub_reg_lake_bd_r="#",
    sl_connected_lake="sl_connected_lake",
    sl_non_connected_lake="sl_nonconnect_lake",
    sl_lakes="selected_lakes",
    sl_str_connected_lake="str_sl_connected_lake",
    nfdr_arcgis="narcgis_fdr",
    nfdr_grass="ngrass_fdr",
    cat_add_lake="cat_add_lake",
    pourpoints_with_lakes="pourpoints_with_lakes",
    cat_use_default_acc="cat_use_default_acc",
    lake_outflow_pourpoints="lake_outflow_pourpoints",
    problem_seg="problem_seg",
    max_memroy=1024 * 4,
):
    work_folder = grassdb
    # define required input files names
    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]

    lake_outflow_pourpoints=Internal_Constant_Names["lake_outflow_pourpoints"]
    # define internal file names
    lake_inflow_pourpoints = Internal_Constant_Names["lake_inflow_pourpoints"]
    catchment_pourpoints_outside_lake = Internal_Constant_Names[
        "catchment_pourpoints_outside_lake"
    ]
    cat_add_lake_old_fdr = Internal_Constant_Names["cat_add_lake_old_fdr"]
    str_connected_lake = Internal_Constant_Names["str_connected_lake"]
    alllake = Internal_Constant_Names["all_lakes"]
    lake_boundary = Internal_Constant_Names["lake_boundary"]
    lakes_lg_cl_thres = "lakes_lg_cl_thres"
    lakes_lg_ncl_thres = "lakes_lg_ncl_thres"


    if not os.path.exists(work_folder):
        os.makedirs(work_folder)
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    cellSize = float(arcpy.GetRasterProperties_management(dem, "CELLSIZEX").getOutput(0))
    SptailRef = arcpy.Describe(dem).spatialReference
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(dem).extent
    arcpy.env.snapRaster =  dem
    lowerLeft = arcpy.Point(arcpy.Describe(dem).extent.XMin,arcpy.Describe(dem).extent.YMin)

    pre_process_lake_polygon(path_lakefile_in,alllake,lake_attributes,lake_boundary,mask,cellSize,SptailRef,work_folder,dem)


    define_cl_and_ncl_lakes(str_r,alllake,str_connected_lake,sl_connected_lake,sl_non_connected_lake,sl_lakes,lake_attributes,threshold_con_lake,threshold_non_con_lake,cellSize,SptailRef,work_folder,dem)

    Lakes_WIth_Multi_Outlet, Remove_Str = create_pour_points_with_lakes(str_r,str_v,cat_no_lake,sl_lakes,sl_connected_lake,acc,pourpoints_with_lakes,
                                      lake_inflow_pourpoints,lake_outflow_pourpoints,
                                      catchment_pourpoints_outside_lake,cellSize,SptailRef,work_folder,dem)

    arcpy.PointToRaster_conversion(pourpoints_with_lakes+"_v", "SubId", pourpoints_with_lakes+"_r")
    outWatershed = Watershed(fdr_arcgis, pourpoints_with_lakes+"_r", "VALUE")
    outWatershed.save(cat_add_lake_old_fdr)

    lakeinfo =  pd.DataFrame.spatial.from_featureclass(lake_outflow_pourpoints + "_v")
    lakeinfo["cat"] = lakeinfo["grid_code"]
    lakeinfo["lmax_acc"] = lakeinfo["acc"]

    cat_withlake_array = arcpy.RasterToNumPyArray(cat_add_lake_old_fdr,nodata_to_value=-9999)
    cat_withlake_array = cat_withlake_array.astype(int)
    fdr_arcgis_array = arcpy.RasterToNumPyArray(fdr_arcgis,nodata_to_value=-9999)
    fdr_arcgis_array = fdr_arcgis_array.astype(int)
    str_r_array = arcpy.RasterToNumPyArray(str_r,nodata_to_value=-9999)
    sl_lakes_array = arcpy.RasterToNumPyArray(sl_lakes+"_r",nodata_to_value=-9999)
    sl_lakes_array = sl_lakes_array.astype(int)
    acc_array = arcpy.RasterToNumPyArray(acc,nodata_to_value=-9999)
    ncols = int(cat_withlake_array.shape[1])
    nrows = int(cat_withlake_array.shape[0])
    lake_boundary_array = arcpy.RasterToNumPyArray(lake_boundary+"_r",nodata_to_value=-9999)


    maximumLakegrids = 1000000000
    pec_grid_outlier = 1
    un_modify_fdr_lakeids = []
    outlakeids, chandir, ndir, bd_problem = modify_lakes_flow_direction(
        cat_withlake_array,
        sl_lakes_array,
        acc_array,
        fdr_arcgis_array,
        str_r_array,
        lakeinfo,
        nrows,
        ncols,
        lake_boundary_array,
        pec_grid_outlier,
        maximumLakegrids,
        un_modify_fdr_lakeids,
    )

    ndir_raster = arcpy.NumPyArrayToRaster(ndir,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    ndir_raster.save(nfdr_arcgis)

    chandir_raster = arcpy.NumPyArrayToRaster(chandir,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    chandir_raster.save("chandir")

    bd_problem_raster = arcpy.NumPyArrayToRaster(bd_problem,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    bd_problem_raster.save("bd_problem")

    outWatershed2 = Watershed(nfdr_arcgis, pourpoints_with_lakes+"_r", "VALUE")
    outWatershed2.save(cat_add_lake)

    print("Following lake have multi outlet ")
    print(Lakes_WIth_Multi_Outlet)
    print("following str are corrected to make one lake one outlet")
    print(Remove_Str)


    return
