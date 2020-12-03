from extent.usingdem import define_project_extent_using_dem


def define_project_extent(
    grassdb,
    grass_location,
    qgis_prefix_path,
    mode,
    path_dem_in,
    outlet_pt=[-1, -1],
    path_extent_ply="#",
    buffer_distance=0.0,
    hybasin_ply="#",
    down_hybasin_id=-1,
    up_hybasin_id=-1,
    mask="MASK",
    dem="dem",
):
    """Define processing extent

    Function that used to define project processing spatial extent (PSE).
    The processing spatial extent is a region where Toolbox will work in. Toolbox
    will not process grids or features outside the processing spatial extent.
    Several options is available here. 1) The PSE can be defined as the extent of
    input DEM. 2)The PSE can be defined using Hybasin product and a hydrobasin
    ID. All subbasin drainage to that hydrobasin ID will be extracted. And
    the extent of the extracted polygon will be used as PSE. 3)The PSE
    can be defined using DEM and an point coordinates. the drainage area
    contribute to that point coordinate will be used as boundary polygon. 4)
    The PSE can be defined using

    Parameters
    ----------
    grassdb                           : path (required)
        It is a path to project grass database folder
    grass_location                    : string (required)
        It is a string of grass location name
    qgis_prefix_path                  : string (required)
        It is a string of qgis prefix path
    mode                              : string (required)
        It is a string indicate which method to define project processing
        spatial extent
        'using_dem'            : the extent of input dem will be used
        'using_hybasin'        : the extent will be defined using hydrobasin
                                 product
        'using_outlet_pt'      : the extent will be defined with provided outlet
                                 point
        'using_provided_ply'   : the extent will be defined by provided polygon
    path_dem_in                      : string (required)
        It is the path to input dem
    outlet_pt                        : list (optional)
        It is list that indicate the outlet coordinates of the
        region of interest. If it is provided, the PSE
        will be defined as the drainage area controlled by this point.
    path_extent_ply                  : string (optional)
        It is the path of a subregion polygon. It is only used when the Region
        of interest is very large and the resolution of the dem is very high.
        toolbox will first divide the whole region into several small subregions.
        And then using devided subregion polygon as PSE.
    buffer_distance                  : float (optional)
        It is a float number to increase the extent of the PSE
        obtained from Hydrobasins. It is needed when input DEM is not from
        HydroSHEDS. Then the extent of the watershed will be different
        with PSE defined by HydroBASINS.
    hybasin_ply                      : string (optional)
        It is a path to hydrobasin routing product, If it is provided, the
        PSE will be based on the OutHyID and OutHyID2 and
        this HydroBASINS routing product.
    down_hybasin_id                  : int (optional)
        It is a HydroBASINS subbasin ID, which should be the ID of the most
        downstream subbasin in the region of interest.
    up_hybasin_id                    : int (optional)
        It is a HydroBASINS subbasin ID, which should be the ID of the most
        upstream subbasin in the region of interest, normally do not needed.
    mask                             : string (optional)
        It is a output mask name, which will stored in grass_location in both
        vector and raster format
    dem                              : string (optional)
        It is a output dem raster name, which will be stored in grass_location

    Notes
    -------
    Outputs are following files

    MASK                   : raster
        it is a mask raster stored in grass database, which indicate
        the PSE. The grass database is located at
        os.path.join(grassdb, grass_location)
    dem                   : raster
        it is a dem raster stored in grass database, which is
        has the same extent with MASK. The grass database is located at
        os.path.join(grassdb, grass_location)

    Returns:
    -------
       None

    Examples
    -------
    """

    if mode == "using_hybasin":
        assert (
            hybasin_ply != "#"
        ), "HydroBasin product is needed to define processing extent with mode = using_hybasin"
        assert (
            down_hybasin_id != "#"
        ), "The watershed outlet HydroBasin sub Id is needed to define processing extent with mode = using_hybasin"

    if mode == "using_outlet_pt":
        assert (
            outlet_pt[0] != "-1"
        ), "outlet_pt needs to be define processing extent with mode = using_outlet_pt"

    if mode == "using_provided_ply":
        assert (
            outlet_pt[0] != "-1"
        ), "outlet_pt needs to be define processing extent with mode = using_outlet_pt"

    if mode == "using_dem":
        define_project_extent_using_dem(
            grassdb,
            grass_location,
            qgis_prefix_path,
            path_dem_in,
            mask=mask,
            dem=dem,
        )
