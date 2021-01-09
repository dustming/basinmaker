from utilities.utilities import Internal_Constant_Names


def add_lakes_and_obs_into_existing_watershed_delineation(
    input_geo_names,
    path_lakefile_in="#",
    lake_attributes=[],
    path_obsfile_in="#",
    obs_attributes=[],
    path_sub_reg_outlets_v="#",
    threshold_con_lake=0,
    threshold_non_con_lake=0,
    search_radius=100,
    sl_connected_lake="sl_connected_lake",
    sl_non_connected_lake="sl_nonconnect_lake",
    sl_lakes="selected_lakes",
    catchment_without_merging_lakes="catchment_without_merging_lakes",
    river_without_merging_lakes="river_without_merging_lakes",
    snapped_obs_points="snapped_obs_points",
    cat_use_default_acc="cat_use_default_acc",
    max_memroy=1024 * 4,
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):
    # define internal file names
    cat_add_lake = Internal_Constant_Names["cat_add_lake"]
    pourpoints_with_lakes = Internal_Constant_Names["pourpoints_with_lakes"]
    pourpoints_add_obs = Internal_Constant_Names["pourpoints_add_obs"]
    alllake = Internal_Constant_Names["all_lakes"]
    lake_boundary = Internal_Constant_Names["lake_boundary"]
    connected_lake = Internal_Constant_Names["connect_lake"]
    non_connected_lake = Internal_Constant_Names["nonconnect_lake"]
    str_connected_lake = Internal_Constant_Names["str_connected_lake"]
    sl_str_connected_lake = Internal_Constant_Names["str_sl_connected_lake"]
    nfdr_arcgis = Internal_Constant_Names["nfdr_arcgis"]
    nfdr_grass = Internal_Constant_Names["nfdr_grass"]
    lake_outflow_pourpoints = Internal_Constant_Names["lake_outflow_pourpoints"]
    
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
        from addlakeandobs.addlakesqgis import (
            add_lakes_into_existing_watershed_delineation,
        )
        from addlakeandobs.addobsqgis import add_obs_into_existing_watershed_delineation
        from addlakeandobs.definecatrivqgis import (
            define_cat_and_riv_without_merge_lake_cats,
        )
    if path_lakefile_in != "#":
        add_lakes_into_existing_watershed_delineation(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            input_geo_names=input_geo_names,
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            threshold_con_lake=threshold_con_lake,
            threshold_non_con_lake=threshold_non_con_lake,
            alllake=alllake,
            lake_boundary=lake_boundary,
            connected_lake=connected_lake,
            non_connected_lake=non_connected_lake,
            str_connected_lake=str_connected_lake,
            sl_connected_lake=sl_connected_lake,
            sl_non_connected_lake=sl_non_connected_lake,
            sl_lakes=sl_lakes,
            sl_str_connected_lake=sl_str_connected_lake,
            nfdr_arcgis=nfdr_arcgis,
            nfdr_grass=nfdr_grass,
            cat_add_lake=cat_add_lake,
            pourpoints_with_lakes=pourpoints_with_lakes,
            cat_use_default_acc=cat_use_default_acc,
            lake_outflow_pourpoints = lake_outflow_pourpoints,
            max_memroy=max_memroy,
        )
    if path_obsfile_in != "#":
        add_obs_into_existing_watershed_delineation(
            grassdb=grassdb,
            grass_location=grass_location,
            qgis_prefix_path=qgis_prefix_path,
            input_geo_names=input_geo_names,
            path_obsfile_in=path_obsfile_in,
            path_lakefile_in=path_lakefile_in,
            obs_attributes=obs_attributes,
            search_radius=search_radius,
            pourpoints_with_lakes=pourpoints_with_lakes,
            lake_outflow_pourpoints = lake_outflow_pourpoints,
            cat_add_lake = cat_add_lake,
            path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            max_memroy=max_memroy,
            pourpoints_add_obs=pourpoints_add_obs,
            snapped_obs_points=snapped_obs_points,
        )
    else:
        pourpoints_add_obs = pourpoints_with_lakes

    define_cat_and_riv_without_merge_lake_cats(
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        input_geo_names=input_geo_names,
        path_lakefile_in=path_lakefile_in,
        pourpoints_add_obs=pourpoints_add_obs,
        catchment_without_merging_lakes=catchment_without_merging_lakes,
        river_without_merging_lakes=river_without_merging_lakes,
        max_memroy=max_memroy,
    )
