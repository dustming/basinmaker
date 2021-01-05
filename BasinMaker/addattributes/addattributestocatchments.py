def add_attributes_to_catchments(
    input_geo_names,
    path_bkfwidthdepth="#",
    bkfwd_attributes=[],
    path_landuse="#",
    path_landuse_info="#",
    path_lake_ply="#",
    out_cat_name="catchment_without_merging_lakes",
    out_riv_name="river_without_merging_lakes",
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
    projection="EPSG:3573",
    obs_v="obs_snap_r2v",
    obs_r="obs",
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    outlet_obs_id=1,
    path_sub_reg_outlets_v="#",
    output_folder = '#'
):
    columns = [
        "SubId",
        "DowSubId",
        "RivSlope",
        "RivLength",
        "BasSlope",
        "BasAspect",
        "BasArea",
        "BkfWidth",
        "BkfDepth",
        "IsLake",
        "HyLakeId",
        "LakeVol",
        "LakeDepth",
        "LakeArea",
        "Laketype",
        "IsObs",
        "MeanElev",
        "FloodP_n",
        "Q_Mean",
        "Ch_n",
        "DA",
        "Strahler",
        "Seg_ID",
        "Seg_order",
        "Max_DEM",
        "Min_DEM",
        "DA_Obs",
        "DA_error",
        "Obs_NM",
        "SRC_obs",
    ]

    coltypes = [
        "Integer",
        "Integer",
        "Real",
        "Real",
        "Real",
        "Real",
        "Real",
        "Real",
        "Real",
        "Integer",
        "Integer",
        "Real",
        "Real",
        "Real",
        "Integer",
        "Integer",
        "Real",
        "Real",
        "Real",
        "Real",
        "Real",
        "Integer",
        "Integer",
        "Integer",
        "Real",
        "Real",
        "Real",
        "Real",
        "Character",
        "Character",
    ]

    if path_lake_ply != "#":
        dem = input_geo_names["dem"]
        mask = input_geo_names["mask"]
        fdr = input_geo_names["nfdr_grass"]
        acc = input_geo_names["acc"]
        cat_use_default_acc = input_geo_names["cat_use_default_acc"]
        river_r = input_geo_names["river_without_merging_lakes"]
        river_v = input_geo_names["str_v"]
        catchments = input_geo_names["catchment_without_merging_lakes"]
        sl_connected_lake = input_geo_names["sl_connected_lake"]
        sl_non_connected_lake = input_geo_names["sl_nonconnect_lake"]
    else:
        dem = input_geo_names["dem"]
        mask = input_geo_names["mask"]
        fdr = input_geo_names["fdr_grass"]
        acc = input_geo_names["acc"]
        cat_use_default_acc = input_geo_names["cat_no_lake"]
        river_r = input_geo_names["str_r"]
        river_v = input_geo_names["str_v"]
        catchments = input_geo_names["cat_no_lake"]

    if gis_platform == "qgis":
        assert (
            grassdb != "#"
        ), "grass database folder is needed, when gis_platform = qgis"
        assert (
            grass_location != "#"
        ), "grass location name is needed, when gis_platform = qgis"
        assert (
            qgis_prefix_path != "#"
        ), "qgis prefix path is needed, when gis_platform = qgis"
        from addattributes.createattributestemplateqgis import (
            create_catchments_attributes_template_table,
        )
        from addattributes.calculatebasicattributesqgis import (
            calculate_basic_attributes,
        )
        from addattributes.addlakeattributesqgis import add_lake_attributes
        from addattributes.adddaandstreamorder import streamorderanddrainagearea
        from addattributes.addnclcatchmentsinfo import (
            update_non_connected_catchment_info,
        )
        from addattributes.joinpandastoattributesqgis import (
            join_pandas_table_to_vector_attributes,
        )
        from addattributes.exportoutputsqgis import export_files_to_output_folder
        from addattributes.addgaugeattributesqgis import add_gauge_attributes
        from addattributes.calfloodmanningnqgis import calculate_flood_plain_manning_n
        from addattributes.calbkfwidthdepthqgis import (
            calculate_bankfull_width_depth_from_polyline,
        )
        from addattributes.returninterestcatchmentsinfo import (
            return_interest_catchments_info,
        )

        attr_template = create_catchments_attributes_template_table(
            grassdb=grassdb,
            grass_location=grass_location,
            columns=columns,
            catchments=catchments,
        )
        
        attr_basic = calculate_basic_attributes(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            catchments=catchments,
            river_r=river_r,
            river_v=river_v,
            dem=dem,
            fdr=fdr,
            acc=acc,
            cat_use_default_acc=cat_use_default_acc,
            projection=projection,
            catinfo=attr_template,
        )
        
        if path_lake_ply != "#":
            attr_lake = add_lake_attributes(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                sl_connected_lake=sl_connected_lake,
                sl_non_connected_lake=sl_non_connected_lake,
                catchments=catchments,
                path_lake_ply=path_lake_ply,
                catinfo=attr_basic,
            )
        else:
            attr_lake = attr_basic
        
        if obs_r != "#":
            attr_obs = add_gauge_attributes(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                pourpoints="Final_OL_v",
                obs_v="obs_snap_r2v",
                obs_r="obs",
                obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
                catinfo=attr_lake,
            )
        else:
            attr_obs = attr_lake
        
        if outlet_obs_id > 0:
            attr_select = return_interest_catchments_info(
                catinfo=attr_obs,
                outlet_obs_id=outlet_obs_id,
                path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            )
        else:
            attr_select = attr_obs
        
        if path_landuse != "#":
            attr_landuse = calculate_flood_plain_manning_n(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                catinfo=attr_select,
                path_landuse=path_landuse,
                path_landuse_info=path_landuse_info,
                riv_seg="nstr_nfinalcat_F",
            )
        else:
            attr_landuse = attr_select
        
        attr_da = streamorderanddrainagearea(attr_landuse)
        
        if path_bkfwidthdepth != "#":
            attr_bkf = calculate_bankfull_width_depth_from_polyline(
                grassdb=grassdb,
                grass_location=grass_location,
                qgis_prefix_path=qgis_prefix_path,
                path_bkfwidthdepth=path_bkfwidthdepth,
                bkfwd_attributes=bkfwd_attributes,
                catchments=catchments,
                catinfo=attr_da,
                mask=mask,
            )
        else:
            attr_bkf = attr_da
        
        attr_ncl = update_non_connected_catchment_info(attr_bkf)
        
        join_pandas_table_to_vector_attributes(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            vector_name="nstr_nfinalcat_F",
            pd_table=attr_ncl,
            column_types=coltypes,
            columns_names=columns,
        )
        
        join_pandas_table_to_vector_attributes(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            vector_name="Net_cat_F",
            pd_table=attr_ncl,
            column_types=coltypes,
            columns_names=columns,
        )

        export_files_to_output_folder(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            input_riv="nstr_nfinalcat_F",
            input_cat="Net_cat_F",
            output_riv=out_riv_name,
            output_cat=out_cat_name,
            input_lake_path=path_lake_ply,
            obs_v="obs_snap_r2v",
            output_folder = output_folder,
        )
