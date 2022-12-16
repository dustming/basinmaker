from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *



def define_cat_and_riv_without_merge_lake_cats(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_lakefile_in,
    catchment_without_merging_lakes="catchment_without_merging_lakes",
    river_without_merging_lakes="river_without_merging_lakes",
    max_memroy=1024 * 4,
):

    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]
    nfdr_grass = input_geo_names["nfdr_grass"]
    sl_non_connected_lake = input_geo_names["sl_non_connected_lake"]
    pourpoints_add_obs = input_geo_names["pourpoints_add_obs"]
    lake_outflow_pourpoints = input_geo_names["lake_outflow_pourpoints"]

    work_folder = grassdb

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


    cat_outlet =  pd.DataFrame.spatial.from_featureclass("sub_outlet_v")
    lake_inlets =  pd.DataFrame.spatial.from_featureclass("lake_inflow_pourpoints_v")
    lake_outlets =  pd.DataFrame.spatial.from_featureclass("lake_outflow_pourpoints_v")
    obs_points =  pd.DataFrame.spatial.from_featureclass("obs_v")
    headwaters =  pd.DataFrame.spatial.from_featureclass("head_node_v")

    cat_outlet["SubId"] = cat_outlet.index + 1
    cat_outlet["type"]  = "cat"
    cat_outlet = cat_outlet[["SubId","type","SHAPE"]]
    lake_inlets["SubId"] = lake_inlets.index + 1 + len(cat_outlet)
    lake_inlets["type"]  = "lake_inlets"
    lake_inlets = lake_inlets[["SubId","type","SHAPE"]]
    lake_outlets["SubId"] = lake_outlets.index + 1 + len(cat_outlet) + len(lake_inlets)
    lake_outlets["type"]  = "lake_outlets"
    lake_outlets = lake_outlets[["SubId","type","SHAPE"]]
    obs_points["SubId"] = obs_points.index + 1 + len(cat_outlet) + len(lake_inlets) + len(lake_outlets)
    obs_points["type"]  = "obs"
    obs_points = obs_points[["SubId","type","SHAPE"]]
    headwaters["SubId"] = headwaters.index + 1 + len(cat_outlet) + len(lake_inlets) + len(lake_outlets) + len(obs_points)
    headwaters["type"]  = "headwater"
    headwaters = headwaters[["SubId","type","SHAPE"]]

    catchment_pourpoints_outside_lake = pd.concat([cat_outlet, lake_inlets,lake_outlets,obs_points,headwaters], ignore_index=True)
    catchment_pourpoints_outside_lake = catchment_pourpoints_outside_lake.reset_index()
    catchment_pourpoints_outside_lake["SubId"] = catchment_pourpoints_outside_lake.index + 1
    catchment_pourpoints_outside_lake.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","final_pourpoints_v"),overwrite=True,sanitize_columns=False)
    arcpy.PointToRaster_conversion("final_pourpoints_v", 'SubId', "final_pourpoints_r")
    arcpy.RasterToPoint_conversion("final_pourpoints_r", "final_pourpoints_v", "VALUE")

    outWatershed2 = Watershed("nfdr_arcgis", "final_pourpoints_r", "VALUE")
    outWatershed2.save(catchment_without_merging_lakes + "_r")

    riv = Con(IsNull(str_r),-9999,catchment_without_merging_lakes + "_r")
    riv_nonull = SetNull(riv, riv, "Value = -9999")
    riv_nonull.save(river_without_merging_lakes+"_r")

#    arcpy.RasterToPolyline_conversion(in_raster = river_without_merging_lakes+"_r", out_polyline_features = river_without_merging_lakes+"_v",simplify = "NO_SIMPLIFY")
#    StreamToFeature(river_without_merging_lakes+"_r", "nfdr_arcgis", river_without_merging_lakes+"_v","SIMPLIFY")
    arcpy.AddField_management("str_v","fld_div","DOUBLE","#","#","#","#","NULLABLE","NON_REQUIRED","#")
    arcpy.CalculateField_management("str_v", "fld_div", '1', "PYTHON")


    arcpy.RasterToPolygon_conversion(catchment_without_merging_lakes + "_r", os.path.join(work_folder,catchment_without_merging_lakes + "_vt.shp"), "NO_SIMPLIFY", "VALUE")

    # arcpy.AlterField_management(in_table = os.path.join(work_folder,catchment_without_merging_lakes + "_vt.dbf"),
    #                            field =  'gridcode',
    #                            new_field_name = 'gridcode',
    #                            field_type = 'LONG')
    arcpy.Dissolve_management(os.path.join(work_folder,catchment_without_merging_lakes + "_vt.shp"), os.path.join(work_folder,catchment_without_merging_lakes + "_v.shp"), ["gridcode"])

    StreamToFeature(river_without_merging_lakes+"_r", "fdr_arcgis", os.path.join(work_folder,river_without_merging_lakes + "_vt.shp"),"NO_SIMPLIFY")
    arcpy.Dissolve_management(os.path.join(work_folder,river_without_merging_lakes + "_vt.shp"), os.path.join(work_folder,river_without_merging_lakes + "_v.shp"), ["grid_code"])



    return
