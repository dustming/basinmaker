def add_lakes_and_obs_into_existing_watershed_delineation(
    input_geo_names,
    path_lakefile_in="#",
    lake_attributes = [],
    path_obsfile_in="#",
    obs_attributes = [],
    threshold_con_lake = 0,
    threshold_non_con_lake = 0,
    alllake = 'all_lakes',
    lake_boundary ='lake_boundary',
    connected_lake = 'connect_lake', 
    non_connected_lake = 'nonconnect_lake',
    str_connected_lake = 'str_connected_lake', 
    sl_connected_lake = 'sl_connected_lake',  
    sl_non_connected_lake = 'sl_nonconnect_lake', 
    sl_lakes = 'selected_lakes' ,
    sl_str_connected_lake = 'str_sl_connected_lake',
    max_memroy=1024 * 4,
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
):    
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
        from addlakeandobs.addlakesqgis import add_lakes_into_existing_watershed_delineation
        
    if path_lakefile_in != '#':
        add_lakes_into_existing_watershed_delineation(
            grassdb,
            grass_location,
            qgis_prefix_path,
            input_geo_names,
            path_lakefile_in,
            lake_attributes,
            threshold_con_lake,
            threshold_non_con_lake,
            alllake = alllake,
            lake_boundary = lake_boundary,
            connected_lake = connected_lake, 
            non_connected_lake = non_connected_lake,
            str_connected_lake = str_connected_lake, 
            sl_connected_lake = sl_connected_lake,  
            sl_non_connected_lake = sl_non_connected_lake, 
            sl_lakes = sl_lakes ,
            sl_str_connected_lake = sl_str_connected_lake,
            max_memroy = max_memroy,
        )

