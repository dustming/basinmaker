def watershed_delineation_without_lake(
    acc_thresold,
    mode,
    input_geo_names,
    max_memroy=1024 * 4,
    grassdb="#",
    grass_location="#",
    qgis_prefix_path="#",
    gis_platform="qgis",
    fdr_path="#",
    subreg_fdr_path="#",
    subreg_acc_path="#",
    subreg_str_r_path="#",
    subreg_str_v_path="#",
    fdr_arcgis="fdr_arcgis",
    fdr_grass="fdr_grass",
    str_r="str_r",
    str_v="str_v",
    acc="acc",
    cat_no_lake="cat_no_lake",
):
    """Generate a subbasin delineation without considering lake

    Function that used to Generate a subbasin delineation and river
    network using user provied flow accumulation thresthold
    without considering lake.

    Parameters
    ----------
    accthresold       : float
        It is the flow accumulation thresthold, used to determine
        subbsains and river network. Increasing of accthresold will
        increase the size of generated subbasins, reduce the number
        subbasins and reduce the number of generated stream segments
    mode              : string (required)
        It is a string indicate which dataset will be used to delineate
        watershed.
        'usingdem'             : dem is used for delineation
        'usingfdr'             : flow direction data is used for delineation
        'usingsubreg'          : predefined subregion inputs is used for
    input_geo_names    : dict
        it is a dictionary that list the required input file names,should at
        least indicate the name of following items:

        mask                   : raster
            it is a mask raster stored in grass database, which indicate
            the PSE. The grass database is located at
            os.path.join(grassdb, grass_location)
        dem                   : raster
            it is a dem raster stored in grass database, which is
            has the same extent with MASK. The grass database is located at
            os.path.join(grassdb, grass_location)
    max_memroy        : integer
        It is the maximum memeory that allow to be used.
    grassdb           : path (required)
        It is a path to project grass database folder
    grass_location    : string (required)
        It is a string of grass location name
    qgis_prefix_path  : string (required)
        It is a string of qgis prefix path
    gis_platform      : string
        It is a string indicate with gis platform is used:
        'qgis'                : the basinmaker is running within QGIS
        'arcgis'              : the basinmaker is running within ArcGIS
    fdr_path          : string
        It is a string indicate path of flow direction dataset
    subreg_fdr_path   : string
        It is a string indicate path of subregion flow direction dataset
    subreg_acc_path   : string
        It is a string indicate path of subregion flow accumlation dataset
    subreg_str_r_path : string
        It is a string indicate path of subregion stream raster flow
        direction dataset
    subreg_str_v_path : string
        It is a string indicate path of subregion stream vector datasets

    Notes
    -------
    Outputs are following files

    fdr_grass              : raster
        it is a raster represent flow direction dataset, which is
        using 1 - 8 to represent different directions
    fdr_arcgis             : raster
        it is a raster represent flow direction dataset, which is
        using 1,2,4,...64,128 to represent different directions
    str_v                  : vector
        it is a river network in vector format
    str_r                  : raster
        it is a river network in raster format
    cat_no_lake            : raster
         it is the raster represent the delineated subbasins without
         considering lakes
    acc                    : raster
         it is the raster represent the flow accumulation

    Returns:
    -------
       None

    Examples
    -------
    """

    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]

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
        from delineationnolake.watusingdemqgis import (
            delineate_watershed_no_lake_using_dem,
        )
        from delineationnolake.watusingfdrqgis import (
            delineate_watershed_no_lake_using_fdr,
        )
        from delineationnolake.watusingsubregionddata import (
            delineate_watershed_no_lake_using_subregion_data,
        )

    if mode == "usingdem":
        assert dem != "#", "The name of dem is needed to delineate watershed from dem"
        assert (
            acc_thresold > 0
        ), "the flow accumulation thresthold is needed to delineate watershed from dem"

    if mode == "usingfdr":
        assert (
            fdr_path != "#"
        ), "The path of the provided flow direction data is needed to delineate watershed from flow direction"

    if mode == "usingsubreg":
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

    if mode == "usingsubreg":
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
