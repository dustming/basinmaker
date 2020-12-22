def add_lakes_and_obs_into_existing_watershed_delineation(
    mode,
    mask="MASK",
    dem="dem",
    path_lakefile_in="#",
    lake_attributes = [],
    path_obsfile_in="#",
    obs_attributes = [],
    fdr_arcgis="fdr_arcgis",
    fdr_grass="fdr_grass",
    str_r="str_r",
    str_v="str_v",
    acc="acc",
    cat_no_lake="cat_no_lake",
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
            mask,
            dem,
            path_lakefile_in,
            lake_attributes,
            fdr_arcgis,
            fdr_grass,
            str_r,
            str_v,
            acc,
            cat_no_lake,
            max_memroy,
        )

