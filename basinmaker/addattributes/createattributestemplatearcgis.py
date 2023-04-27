from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *



def create_catchments_attributes_template_table(
    grassdb,
    grass_location,
    input_geo_names,
    columns,
    obs_attributes=[],
    lake_attributes=[],
    path_landuse="#",
    path_landuse_info="#",
    path_k_c_zone_polygon = '#',
    output_folder='#',
):

    catchments = input_geo_names["catchment_without_merging_lakes"]
    work_folder = grassdb
    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]
    work_folder = grassdb
    catchment_without_merging_lakes="catchment_without_merging_lakes"
    river_without_merging_lakes="river_without_merging_lakes"

    attri_table = pd.DataFrame({
                       'SubId'     : pd.Series(dtype='int'),
                       'DowSubId'  : pd.Series(dtype='int'),
                       # 'RivSlope'  : pd.Series(dtype='float'),
                       # 'RivLength' : pd.Series(dtype='float'),
                       # 'BasSlope'  : pd.Series(dtype='float'),
                       # 'BasAspect' : pd.Series(dtype='float'),
                       # 'BasArea'   : pd.Series(dtype='int64'),
                       'BkfWidth'  : pd.Series(dtype='float'),
                       'BkfDepth'  : pd.Series(dtype='float'),
                       'Lake_Cat'  : pd.Series(dtype='int'),
                       'HyLakeId'  : pd.Series(dtype='int'),
                       'LakeVol'   : pd.Series(dtype='float'),
                       'LakeDepth' : pd.Series(dtype='float'),
                       'LakeArea'  : pd.Series(dtype='float'),
                       'Laketype'  : pd.Series(dtype='int'),
                       'Has_POI'   : pd.Series(dtype='int'),
                       # 'MeanElev'  : pd.Series(dtype='int'),
                       # 'FloodP_n'  : pd.Series(dtype='float'),
                       'Q_Mean'    : pd.Series(dtype='float'),
                       'Ch_n'      : pd.Series(dtype='float'),
                       'DrainArea' : pd.Series(dtype='int64'),
                       'Strahler'  : pd.Series(dtype='int'),
                       'Seg_ID'    : pd.Series(dtype='int'),
                       'Seg_order' : pd.Series(dtype='int'),
                       # 'Max_DEM'   : pd.Series(dtype='float'),
                       # 'Min_DEM'   : pd.Series(dtype='float'),
                       'DA_Obs'    : pd.Series(dtype='float'),
                       'DA_error'  : pd.Series(dtype='float'),
                       # 'Obs_NM'    : pd.Series(dtype='float'),
                       # 'SRC_obs'   : pd.Series(dtype='float'),
                       # 'centroid_x': pd.Series(dtype='float'),
                       # 'centroid_y': pd.Series(dtype='float'),
                       'DA_Chn_L'  : pd.Series(dtype='float'),
                       'DA_Slope'  : pd.Series(dtype='float'),
                       'DA_Chn_Slp': pd.Series(dtype='float'),
                       # 'outletLat' : pd.Series(dtype='float'),
                       # 'outletLng' : pd.Series(dtype='float'),
                       # 'k'         : pd.Series(dtype='float'),
                       # 'c'         : pd.Series(dtype='float'),
                       }
                      )


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


    ## GIS spatial calculations

    # rasterize landuse data


    arcpy.RasterToPolygon_conversion(catchment_without_merging_lakes + "_r", os.path.join(work_folder,catchment_without_merging_lakes + "_vt.shp"), "NO_SIMPLIFY", "VALUE")

    # arcpy.AlterField_management(in_table = os.path.join(work_folder,catchment_without_merging_lakes + "_vt.dbf"),
    #                            field =  'gridcode',
    #                            new_field_name = 'gridcode',
    #                            field_type = 'LONG')
    arcpy.Dissolve_management(os.path.join(work_folder,catchment_without_merging_lakes + "_vt.shp"), os.path.join(work_folder,catchment_without_merging_lakes + "_v.shp"), ["gridcode"])

    StreamToFeature(river_without_merging_lakes+"_r", "fdr_arcgis", os.path.join(work_folder,river_without_merging_lakes + "_vt.shp"),"NO_SIMPLIFY")
    arcpy.Dissolve_management(os.path.join(work_folder,river_without_merging_lakes + "_vt.shp"), os.path.join(work_folder,river_without_merging_lakes + "_v.shp"), ["grid_code"])


    landuse_in_region = ExtractByMask(path_landuse, dem)
    landuse_in_region.save("landuse")

    arcpy.Project_management(
        path_k_c_zone_polygon,
        "kc_zone_proj",
        arcpy.SpatialReference(int(SptailRef.factoryCode)),
        )
    arcpy.arcpy.analysis.PairwiseClip("kc_zone_proj", mask + '_ply', "kc_zone_proj_clip", "")

    arcpy.PolygonToRaster_conversion("kc_zone_proj", "k", "k")
    arcpy.PolygonToRaster_conversion("kc_zone_proj", "c", "c")


    arcpy.management.CalculateGeometryAttributes(os.path.join(work_folder,river_without_merging_lakes+"_v"), [["Length_m", "LENGTH"]], "Meters")
    arcpy.management.CalculateGeometryAttributes(os.path.join(work_folder,catchment_without_merging_lakes+"_v"), [["Area_m", "AREA"]], "Meters")
    arcpy.management.CalculateGeometryAttributes(
                                                 in_features = os.path.join(work_folder,catchment_without_merging_lakes+"_v"),
                                                 geometry_property = [["centroid_x", "CENTROID_X"],["centroid_y", "CENTROID_Y"]],
                                                 coordinate_format = "DD",
                                                 )

    arcpy.management.CalculateGeometryAttributes(
                                                 in_features = "final_pourpoints_v",
                                                 geometry_property = [["outletLng", "POINT_X"],["outletLat", "POINT_Y"]],
                                                 coordinate_format = "DD",
                                                 )

    slope_degree = Slope(dem, "DEGREE")
    slope_degree.save("slope_degree")
    slope_pec = Slope(dem, "PERCENT_RISE")
    slope_pec.save("slope_pec")
    aspect = Aspect(dem)
    aspect.save("aspect")

    ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value",slope_degree,
                                 "sub_degree", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value",aspect,
                                 "sub_aspect", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value",dem,
                                 "sub_dem", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value","k",
                                 "cat_k", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value","c",
                                 "cat_c", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(river_without_merging_lakes+"_r", "Value",slope_pec,
                                 "riv_slope", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(river_without_merging_lakes+"_r", "Value",dem,
                                 "riv_dem", "DATA", "MIN_MAX_MEAN")
    ZonalStatisticsAsTable(river_without_merging_lakes+"_r", "Value","landuse",
                                 "riv_landuse", "DATA", "MIN_MAX_MEAN")
    ExtractMultiValuesToPoints("final_pourpoints_v", [["sl_connected_lake_r","cl_lake_id"],["sl_nonconnect_lake_r","ncl_kake_id"],["obs_r","obsid"]], "NONE")

    ## END gis calculations  GIS spatial calculations

    # routing information
    StreamToFeature(river_without_merging_lakes+"_r", "nfdr_arcgis","riv_link_v","SIMPLIFY")
    routing_info_strs = return_routing_info_using_str_v("riv_link_v",cellSize,SptailRef,work_folder,dem)
    routing_info_strs = routing_info_strs.sort_values(by=['SubId','DowSubId'], ascending=True)
    routing_info_strs = routing_info_strs.drop_duplicates(subset=['SubId'], keep='last')

    final_pourpoints = pd.DataFrame.spatial.from_featureclass("final_pourpoints_v")
    final_pourpoints["SubId"] = final_pourpoints["grid_code"].astype(int)
    final_pourpoints = final_pourpoints.merge(routing_info_strs,on='SubId',how='left')
    pourpoints_no_dowsubid = final_pourpoints[final_pourpoints['DowSubId'].isnull()]
    pourpoints_no_dowsubid.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","final_pp_no_downsubid"),overwrite=True,sanitize_columns=False)
    arcpy.analysis.Buffer(in_features = "final_pp_no_downsubid", out_feature_class = "final_pp_no_downsubid_bf",buffer_distance_or_field = str(2.5*cellSize) + "  Meters",method="PLANAR")
    ZonalStatisticsAsTable("final_pp_no_downsubid_bf", "SubId",catchment_without_merging_lakes+"_r" ,
                                 "routing_table_not_on_river", "DATA", "ALL",percentile_values = [90,80,70,60,50,40,30,20,10,5,1,0.1,0.01])

    routing_for_pp_not_on_river = read_table_as_pandas("routing_table_not_on_river",['SubId','VARIETY','PCT90','PCT80','PCT70','PCT60','PCT50','PCT40','PCT30','PCT20','PCT10','PCT5','PCT1','PCT0_1','PCT0_01','MIN','MAX'],work_folder)
    routing_for_pp_not_on_river['DowSubId2'] = -1.2345
    mask = routing_for_pp_not_on_river['VARIETY'] == 1
    routing_for_pp_not_on_river.loc[mask,'DowSubId2'] = -1
    loop_pct = ['PCT90','PCT80','PCT70','PCT60','PCT50','PCT40','PCT30','PCT20','PCT10','PCT5','PCT1','PCT0_1','PCT0_01','MIN','MAX']

    for pct in loop_pct:
        mask_notprocessed    = routing_for_pp_not_on_river['DowSubId2'] == -1.2345
        mask_pck = routing_for_pp_not_on_river[pct] != routing_for_pp_not_on_river['SubId']
        maskfinal = np.logical_and(mask_pck,mask_notprocessed)
        routing_for_pp_not_on_river.loc[maskfinal,'DowSubId2'] = routing_for_pp_not_on_river.loc[maskfinal,pct]

    not_processed = routing_for_pp_not_on_river[routing_for_pp_not_on_river['DowSubId2'] == -1.2345].copy(deep=True)
    if len(not_processed) > 0:
        print(not_processed)
        print("contact Ming for this error")
        quit()

    final_pourpoints = final_pourpoints.merge(routing_for_pp_not_on_river,on='SubId',how='left')
    final_pourpoints.loc[final_pourpoints['DowSubId'].isnull(),'DowSubId'] = final_pourpoints.loc[final_pourpoints['DowSubId'].isnull(),'DowSubId2']

    outlet_mask = final_pourpoints['DowSubId'] == final_pourpoints['SubId']
    final_pourpoints.loc[outlet_mask,'DowSubId'] = -1

    final_pourpoints["HyLakeId"] = final_pourpoints["cl_lake_id"]
    final_pourpoints.loc[final_pourpoints['HyLakeId'].isnull(),'HyLakeId'] = final_pourpoints.loc[final_pourpoints['HyLakeId'].isnull(),'ncl_kake_id']

    alllake123 = pd.DataFrame.spatial.from_featureclass("all_lakes_v")
    alllake123['HyLakeId'] = alllake123[lake_attributes[0]]
    alllake123 = alllake123[lake_attributes + ['HyLakeId']]
    obs_v123 = pd.DataFrame.spatial.from_featureclass("obs_clip")
    obs_v123["obsid"] = obs_v123[obs_attributes[0]]
    obs_v123 = obs_v123[obs_attributes + ['obsid']]

    final_pourpoints = final_pourpoints.merge(alllake123,on='HyLakeId',how='left')
    final_pourpoints = final_pourpoints.merge(obs_v123,on='obsid',how='left')
    final_pourpoints.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","final_pp_with_routing"),overwrite=True,sanitize_columns=False)

    final_pourpoints = final_pourpoints.sort_values(by='SubId', ascending=True)

    attri_table['SubId'] = final_pourpoints['SubId']
    attri_table = attri_table.fillna(-1.2345)
    attri_table['DowSubId'] = final_pourpoints['DowSubId'].fillna(0).astype(int)
    attri_table['HyLakeId'] = final_pourpoints['HyLakeId'].fillna(0).astype(int)
    attri_table['Lake_Cat'] = 0
    attri_table['LakeVol']  = final_pourpoints[lake_attributes[3]].fillna(0)
    attri_table['LakeDepth'] = final_pourpoints[lake_attributes[4]].fillna(0)
    attri_table['LakeArea'] =final_pourpoints[lake_attributes[2]].fillna(0)*1000*1000
    attri_table['Laketype'] = final_pourpoints[lake_attributes[1]].fillna(0)
    attri_table.loc[~final_pourpoints['ncl_kake_id'].isnull(),'Lake_Cat'] = 2
    attri_table.loc[~final_pourpoints['cl_lake_id'].isnull(),'Lake_Cat'] = 1
    attri_table['Has_POI'] = int(0)
    attri_table.loc[~final_pourpoints['obsid'].isnull(),'Has_POI'] = int(1)
    attri_table['DA_Obs'] = final_pourpoints[obs_attributes[2]].fillna(0)
    attri_table['Obs_NM'] = final_pourpoints[obs_attributes[1]].astype('str').fillna(" ")
    attri_table['SRC_obs'] = final_pourpoints[obs_attributes[3]].astype('str').fillna(" ")
    attri_table['outletLng'] = final_pourpoints['outletLng'].fillna(-1.2345)
    attri_table['outletLat'] = final_pourpoints['outletLat'].fillna(-1.2345)


    sub_attri = read_table_as_pandas("sub_degree",['Value','MEAN'],work_folder)
    sub_attri['BasSlope'] = sub_attri['MEAN'].fillna(-1.2345)
    sub_attri['SubId'] = sub_attri['Value']
    sub_attri = sub_attri.drop(columns=['MEAN', 'Value'])
    attri_table = attri_table.merge(sub_attri,on='SubId',how='left')


    sub_asp = read_table_as_pandas("sub_aspect",['Value','MEAN'],work_folder)
    sub_asp['BasAspect'] = sub_asp['MEAN'].fillna(-1.2345)
    sub_asp['SubId'] = sub_asp['Value']
    sub_asp = sub_asp.drop(columns=['MEAN', 'Value'])
    attri_table = attri_table.merge(sub_asp,on='SubId',how='left')

    sub_dem = read_table_as_pandas("sub_dem",['Value','MEAN'],work_folder)
    sub_dem['MeanElev'] = sub_dem['MEAN'].fillna(-1.2345)
    sub_dem['SubId'] = sub_dem['Value']
    sub_dem = sub_dem.drop(columns=['MEAN', 'Value'])
    attri_table = attri_table.merge(sub_dem,on='SubId',how='left')


    sub_c = read_table_as_pandas("cat_c",['Value','MEAN'],work_folder)
    sub_c['c'] = sub_c['MEAN'].fillna(-1.2345)
    sub_c['SubId'] = sub_c['Value']
    sub_c = sub_c.drop(columns=['MEAN', 'Value'])
    attri_table = attri_table.merge(sub_c,on='SubId',how='left')

    sub_k = read_table_as_pandas("cat_k",['Value','MEAN'],work_folder)
    sub_k['k'] = sub_k['MEAN'].fillna(-1.2345)
    sub_k['SubId'] = sub_k['Value']
    sub_k = sub_k.drop(columns=['MEAN', 'Value'])
    attri_table = attri_table.merge(sub_k,on='SubId',how='left')
#    print(len(attri_table),"f")

    cat_ply = pd.DataFrame.spatial.from_featureclass(os.path.join(work_folder,catchment_without_merging_lakes+"_v"))
    cat_ply['SubId'] = cat_ply['gridcode']
    cat_ply['BasArea'] = cat_ply['Area_m'].fillna(-1.2345)
    cat_ply2 = cat_ply[['SubId','BasArea','centroid_x','centroid_y']].copy(deep=True)
    attri_table = attri_table.merge(cat_ply2,on='SubId',how='left')
#    print(len(attri_table),"e")
#    print(len(cat_ply))

    riv_line = pd.DataFrame.spatial.from_featureclass(os.path.join(work_folder,river_without_merging_lakes+"_v"))
    riv_line['SubId'] = riv_line['grid_code']
    riv_line['RivLength'] = riv_line['Length_m'].fillna(-1.2345)
    riv_line2 = riv_line[['SubId','RivLength']].copy(deep=True)
    attri_table = attri_table.merge(riv_line2,on='SubId',how='left')
#    print(len(attri_table),"d")


    riv_attri = read_table_as_pandas("riv_slope",['Value','MEAN'],work_folder)
    riv_attri['SubId'] = riv_attri['Value']
    riv_attri = riv_attri.drop(columns=['MEAN', 'Value'])
    attri_table = attri_table.merge(riv_attri,on='SubId',how='left')
#    print(len(attri_table),"c")
    riv_dem = read_table_as_pandas("riv_dem",['Value','MIN','MAX'],work_folder)
    riv_dem['Max_DEM'] = riv_dem['MAX'].fillna(-1.2345)
    riv_dem['Min_DEM'] = riv_dem['MIN'].fillna(-1.2345)
    riv_dem['SubId'] = riv_dem['Value']
    riv_dem = riv_dem.drop(columns=['MIN', 'MAX','Value'])
    attri_table = attri_table.merge(riv_dem,on='SubId',how='left')
    attri_table['RivSlope'] = (attri_table['Max_DEM'] - attri_table['Min_DEM']) / attri_table['RivLength']
#    print(len(attri_table),"b")
    riv_landuse = read_table_as_pandas("riv_landuse",['Value','MEAN'],work_folder)
    riv_landuse['FloodP_n'] = riv_landuse['MEAN'] / 1000.0
    riv_landuse['FloodP_n'] = riv_landuse['FloodP_n'].fillna(-1.2345)
    riv_landuse['SubId'] = riv_landuse['Value']
    riv_landuse = riv_landuse.drop(columns=['MEAN','Value'])
    attri_table = attri_table.merge(riv_landuse,on='SubId',how='left')

    attri_table = attri_table.drop_duplicates(subset=['SubId'], keep='first')
    attri_table = streamorderanddrainagearea(attri_table)
    attri_table = calculate_bkf_width_depth(attri_table)

    attri_table = update_non_connected_catchment_info(attri_table)
    attri_table.loc[attri_table['DA_Obs'] > 0,'DA_error'] = (attri_table.loc[attri_table['DA_Obs'] > 0,'DrainArea']/1000/1000)/attri_table.loc[attri_table['DA_Obs'] > 0,'DA_Obs']
    attri_table.loc[attri_table['RivLength'] == -1.2345,'RivSlope'] = -1.2345
    attri_table.loc[attri_table['RivLength'] == -1.2345,'FloodP_n'] = -1.2345
    attri_table.loc[attri_table['RivLength'] == -1.2345,'Max_DEM'] = -1.2345
    attri_table.loc[attri_table['RivLength'] == -1.2345,'Min_DEM'] = -1.2345
    attri_table.loc[attri_table['RivLength'] == -1.2345,'Ch_n'] = -1.2345

    cat_ply= cat_ply[['SubId','SHAPE']]
    cat_ply = cat_ply.merge(attri_table,on='SubId',how='left')
    cat_ply = cat_ply.fillna(-1.2345)
    cat_ply['Obs_NM'] = cat_ply['Obs_NM'].astype('str')
    cat_ply['SRC_obs'] = cat_ply['SRC_obs'].astype('str')
    cat_ply = cat_ply.drop_duplicates(subset=['SubId'], keep='first')
    cat_ply.spatial.to_featureclass(location=os.path.join(output_folder,catchment_without_merging_lakes+"_v1-0"),overwrite=True,sanitize_columns=False)

    riv_line= riv_line[['SubId','SHAPE']]
    riv_line = riv_line.merge(attri_table,on='SubId',how='left')
    riv_line = riv_line.fillna(-1.2345)
    riv_line['Obs_NM'] = riv_line['Obs_NM'].astype('str')
    riv_line['SRC_obs'] = riv_line['SRC_obs'].astype('str')
    riv_line = riv_line.drop_duplicates(subset=['SubId'], keep='first')
    riv_line.spatial.to_featureclass(location=os.path.join(output_folder,river_without_merging_lakes+"_v1-0"),overwrite=True,sanitize_columns=False)

    
    arcpy.FeatureClassToFeatureClass_conversion("sl_connected_lake_v", output_folder,"sl_connected_lake_v1-0.shp")
    arcpy.FeatureClassToFeatureClass_conversion("sl_nonconnect_lake_v",output_folder,"sl_non_connected_lake_v1-0.shp")

    arcpy.DeleteField_management(os.path.join(output_folder,"sl_connected_lake_v1-0.shp"),
        ["area_ratio", "Id","ORIG_FID", "BUFF_DIST","newarea0","newarea1","Shape_Leng","Shape_Area"]
    )
    arcpy.DeleteField_management(os.path.join(output_folder,"sl_non_connected_lake_v1-0.shp"),
        ["area_ratio", "Id","ORIG_FID", "BUFF_DIST","newarea0","newarea1","Shape_Leng","Shape_Area"]
    )

    obs_v = pd.DataFrame.spatial.from_featureclass("obs_v")
    obs_v['obsid'] = obs_v['grid_code']
    obs_v2 = pd.DataFrame.spatial.from_featureclass("obs_clip")
    obs_v2['obsid'] = obs_v2['Obs_ID']
    obs_v_missing = obs_v2[~obs_v2['obsid'].isin(final_pourpoints['obsid'].values)].copy(deep=True)

    final_pourpoints_2 = final_pourpoints[['obsid','SubId']]
    obs_v = obs_v.merge(final_pourpoints_2,on='obsid',how='left')
    obs_v = obs_v[['SubId','SHAPE']]
    cat_ply_att = cat_ply[['SubId','DA_Obs','SRC_obs','DrainArea','DA_error','Obs_NM']]
    obs_v = obs_v.merge(cat_ply_att,on='SubId',how='left')
    obs_v = obs_v[['SubId','DA_Obs','SRC_obs','DrainArea','DA_error','Obs_NM','SHAPE']]
    obs_v.spatial.to_featureclass(location=os.path.join(output_folder,"poi_v1-0.shp"),overwrite=True,sanitize_columns=False)
    if len(obs_v_missing) > 0:
        obs_v_missing.spatial.to_featureclass(location=os.path.join(output_folder,"poi_missing.shp"),overwrite=True,sanitize_columns=False)

    arcpy.DeleteField_management(os.path.join(output_folder,"poi_v1-0.shp"),
                             ["Id"])
    arcpy.DeleteField_management(os.path.join(output_folder,catchment_without_merging_lakes+"_v1-0"),
                             ["Id"])
    arcpy.DeleteField_management(os.path.join(output_folder,river_without_merging_lakes+"_v1-0"),
                             ["Id"])

    return
