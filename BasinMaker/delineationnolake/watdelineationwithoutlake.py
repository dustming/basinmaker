def watershed_delineation_without_lake(
    mode,
    mask="MASK",
    dem="dem",
    acc_thresold=-1,
    fdr_path="#",
    subreg_fdr_path="#",
    subreg_acc_path="#",
    subreg_str_r_path="#",
    subreg_str_v_path="#",
    fdr_arcgis="",
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
        from delineation.watusingdemqgis import delineate_watershed_no_lake_using_dem
        from delineation.watusingfdrqgis import delineate_watershed_no_lake_using_fdr
        from delineation.watusingsubregionddata import delineate_watershed_no_lake_using_subregion_data
        
    if mode == "dem":
        assert dem != "#", "The name of dem is needed to delineate watershed from dem"
        assert (
            acc_thresold > 0
        ), "the flow accumulation thresthold is needed to delineate watershed from dem"

    if mode == "fdr":
        assert (
            fdr_path != "#"
        ), "The path of the provided flow direction data is needed to delineate watershed from flow direction"

    if mode == "subreg":
        assert (
            subreg_acc_path != "#"
        ), "subregion flow accumulation is neeeded to delineate watershed for current subregion "

        assert (
            subreg_fdr_path != "#"
        ), "subregion flow direction is neeeded to delineate watershed for current subregion "

        assert (
            subreg_str_r_path != "#"
        ), "subregion stream raster is neeeded to delineate watershed for current subregion "

        assert (
            subreg_str_v_path != "#"
        ), "subregion stream vector is neeeded to delineate watershed for current subregion "

    if mode == "usingdem":
        delineate_watershed_no_lake_using_dem(
            grassdb,
            grass_location,
            qgis_prefix_path,
            mask,
            dem,
            acc_thresold,
            fdr_arcgis,
            fdr_grass,
            str_r,
            str_v,
            acc,
            cat_no_lake,
            max_memroy,
        )

    if mode == "usingfdr":
        delineate_watershed_no_lake_using_fdr(
            grassdb,
            grass_location,
            qgis_prefix_path,
            mask,
            fdr_path,
            acc_thresold,
            fdr_arcgis,
            fdr_grass,
            str_r,
            str_v,
            acc,
            cat_no_lake,
            max_memroy,
        )
        
    if mode == "usingssubregiondata":
        delineate_watershed_no_lake_using_subregion_data(
            grassdb,
            grass_location,
            qgis_prefix_path,
            mask,
            subreg_fdr_path,
            subreg_acc_path,
            subreg_str_v_path,
            subreg_str_r_path,
            acc_thresold,
            fdr_arcgis,
            fdr_grass,
            str_r,
            str_v,
            acc,
            cat_no_lake,
            max_memroy,
        )
