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
