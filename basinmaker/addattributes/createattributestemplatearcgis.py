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
                       'RivSlope'  : pd.Series(dtype='float'),
                       'RivLength' : pd.Series(dtype='float'),
                       'BasSlope'  : pd.Series(dtype='float'),
                       'BasAspect' : pd.Series(dtype='float'),
                       'BasArea'   : pd.Series(dtype='int64'),
                       'BkfWidth'  : pd.Series(dtype='float'),
                       'BkfDepth'  : pd.Series(dtype='float'),
                       'Lake_Cat'  : pd.Series(dtype='int'),
                       'HyLakeId'  : pd.Series(dtype='int'),
                       'LakeVol'   : pd.Series(dtype='float'),
                       'LakeDepth' : pd.Series(dtype='float'),
                       'LakeArea'  : pd.Series(dtype='float'),
                       'Laketype'  : pd.Series(dtype='int'),
                       'Has_POI'   : pd.Series(dtype='int'),
                       'MeanElev'  : pd.Series(dtype='int'),
                       'FloodP_n'  : pd.Series(dtype='float'),
                       'Q_Mean'    : pd.Series(dtype='float'),
                       'Ch_n'      : pd.Series(dtype='float'),
                       'DrainArea' : pd.Series(dtype='int64'),
                       'Strahler'  : pd.Series(dtype='int'),
                       'Seg_ID'    : pd.Series(dtype='int'),
                       'Seg_order' : pd.Series(dtype='int'),
                       'Max_DEM'   : pd.Series(dtype='float'),
                       'Min_DEM'   : pd.Series(dtype='float'),
                       'DA_Obs'    : pd.Series(dtype='float'),
                       'DA_error'  : pd.Series(dtype='float'),
                       'Obs_NM'    : pd.Series(dtype='float'),
                       'SRC_obs'   : pd.Series(dtype='float'),
                       'centroid_x': pd.Series(dtype='float'),
                       'centroid_y': pd.Series(dtype='float'),
                       'DA_Chn_L'  : pd.Series(dtype='float'),
                       'DA_Slope'  : pd.Series(dtype='float'),
                       'DA_Chn_Slp': pd.Series(dtype='float'),
                       'outletLat' : pd.Series(dtype='float'),
                       'outletLng' : pd.Series(dtype='float'),
                       'k'         : pd.Series(dtype='float'),
                       'c'         : pd.Series(dtype='float'),
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



    StreamToFeature(river_without_merging_lakes+"_r", "fdr_arcgis", river_without_merging_lakes+"_vt","SIMPLIFY")
    arcpy.analysis.PairwiseDissolve(river_without_merging_lakes+"_vt", river_without_merging_lakes+"_v", "grid_code")

    # extract lake information to the table

#    ExtractMultiValuesToPoints("final_pourpoints_v", [["sl_connected_lake_r","cl_lake_id"],["sl_nonconnect_lake_r","ncl_kake_id"],["obs_r","obsid"]], "NONE")


    # routing information
    StreamToFeature(river_without_merging_lakes+"_r", "nfdr_arcgis","riv_link_v","SIMPLIFY")
    routing_info_strs = return_routing_info_using_str_v("riv_link_v",cellSize,SptailRef,work_folder,dem)
    routing_info_strs = routing_info_strs.sort_values(by=['SubId','DowSubId'], ascending=True)
    routing_info_strs = routing_info_strs.drop_duplicates(subset=['SubId'], keep='last')

    final_pourpoints = pd.DataFrame.spatial.from_featureclass("final_pourpoints_v")
    final_pourpoints["SubId"] = final_pourpoints["grid_code"].astype(int)
    # final_pourpoints = pd.concat([final_pourpoints,routing_info_strs],axis=1,ignore_index=True)
    final_pourpoints = final_pourpoints.merge(routing_info_strs,on='SubId',how='left')
    pourpoints_no_dowsubid = final_pourpoints[final_pourpoints['DowSubId'].isnull()]

    pourpoints_no_dowsubid.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","final_pp_no_downsubid"),overwrite=True,sanitize_columns=False)
    arcpy.analysis.Buffer(in_features = "final_pp_no_downsubid", out_feature_class = "final_pp_no_downsubid_bf",buffer_distance_or_field = str(2.5*cellSize) + "  Meters",method="PLANAR")
    ZonalStatisticsAsTable("final_pp_no_downsubid_bf", "SubId",catchment_without_merging_lakes+"_r" ,
                                 "routing_table_not_on_river", "NODATA", "ALL")
    routing_for_pp_not_on_river = read_table_as_pandas("routing_table_not_on_river",['SubId','MIN','MAX'],work_folder)
    mask = routing_for_pp_not_on_river['SubId'] == routing_for_pp_not_on_river['MIN']
    routing_for_pp_not_on_river['DowSubId2'] = routing_for_pp_not_on_river['SubId']
    routing_for_pp_not_on_river['DowSubId2'] = -1
    routing_for_pp_not_on_river.loc[mask,'DowSubId2'] = routing_for_pp_not_on_river.loc[mask,'MAX']
    routing_for_pp_not_on_river.loc[~mask,'DowSubId2'] = routing_for_pp_not_on_river.loc[~mask,'MIN']
    final_pourpoints = final_pourpoints.merge(routing_for_pp_not_on_river,on='SubId',how='left')
    final_pourpoints.loc[final_pourpoints['DowSubId'].isnull(),'DowSubId'] = final_pourpoints.loc[final_pourpoints['DowSubId'].isnull(),'DowSubId2']

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

    attri_table['SubId'] = final_pourpoints['SubId']
    attri_table['DowSubId'] = final_pourpoints['DowSubId'].fillna(0).astype(int)
    attri_table['HyLakeId'] = final_pourpoints['HyLakeId'].fillna(0).astype(int)
    attri_table['Lake_Cat'] = 0
    attri_table['LakeVol']  = final_pourpoints[lake_attributes[3]].fillna(0)
    attri_table['LakeDepth'] = final_pourpoints[lake_attributes[4]].fillna(0)
    attri_table['LakeArea'] =final_pourpoints[lake_attributes[2]].fillna(0)
    attri_table['Laketype'] = final_pourpoints[lake_attributes[1]].fillna(0)
    attri_table.loc[~final_pourpoints['ncl_kake_id'].isnull(),'Lake_Cat'] = 2
    attri_table.loc[~final_pourpoints['cl_lake_id'].isnull(),'Lake_Cat'] = 1
    attri_table['Has_POI'] = int(0)
    attri_table.loc[~final_pourpoints['obsid'].isnull(),'Has_POI'] = int(1)
    attri_table['DA_Obs'] = final_pourpoints[obs_attributes[2]]
    attri_table['Obs_NM'] = final_pourpoints[obs_attributes[1]]
    attri_table['SRC_obs'] = final_pourpoints[obs_attributes[3]]

    # slope_degree = Slope(dem, "DEGREE")
    # slope_degree.save("slope_degree")
    # slope_pec = Slope(dem, "PERCENT_RISE")
    # slope_pec.save("slope_pec")
    # aspect = Aspect(dem)
    # aspect.save("aspect")
    #
    # ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value",slope_degree,
    #                              "sub_degree", "NODATA", "MIN_MAX_MEAN")
    # ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value",aspect,
    #                              "sub_aspect", "NODATA", "MIN_MAX_MEAN")
    # ZonalStatisticsAsTable(catchment_without_merging_lakes+"_r", "Value",aspect,
    #                              "sub_dem", "NODATA", "MIN_MAX_MEAN")
    # ZonalStatisticsAsTable(river_without_merging_lakes+"_r", "Value",slope_pec,
    #                              "riv_slope", "NODATA", "MIN_MAX_MEAN")

    print(attri_table)
    adf

    return
