import os
from extent.projectextent import define_project_extent

class BasinMakerQGIS:

    """
    QGIS/GRASSToolsets to delelineate lake river routing structure
    ...

    Attributes
    ----------


    Methods
    -------

    """

    def __init__(
        self,
        path_output_folder,
        path_working_folder,
    ):
        
        # define drived values
        # create folders
        self.path_output_folder = path_output_folder
        self.path_working_folder = path_working_folder
        
        os.makedirs(self.path_output_folder, exist_ok  = True)
        os.makedirs(self.path_working_folder, exist_ok = True)

        # obtain qgis prefix path
        self.qgispp = os.environ["QGISPrefixPath"]
        # obtain basinmaker path
        self.basinmaker_path = os.environ["RoutingToolFolder"]

        # define grass database folder
        self.grassdb = os.path.join(self.path_working_folder, "grassdb")
        if not os.path.exists(self.grassdb):
            os.makedirs(self.grassdb)

        # add grass database folder into enviroment variable
        os.environ["GISDBASE"] = self.grassdb

        # define grass location names
        self.grass_location_geo = "main_working_location"

        # grass sql databse folder
        self.sqlpath = os.path.join(
            self.grassdb, "main_working_location", "PERMANENT", "sqlite", "sqlite.db"
        )

        # define constants

        # field names in the attribute table
        self.field_names = [
            "SubId",
            "HRU_IsLake",
            "Landuse_ID",
            "Soil_ID",
            "Veg_ID",
            "O_ID_1",
            "O_ID_2",
            "HRU_Area",
            "HRU_ID",
            "LAND_USE_C",
            "VEG_C",
            "SOIL_PROF",
            "HRU_CenX",
            "HRU_CenY",
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
            "HRU_S_mean",
            "HRU_A_mean",
            "HRU_E_mean",
            "centroid_x",
            "centroid_y",
            "Rivlen",
            "area",
        ]

        # default channel manning's coefficient
        DEFAULT_CHN = 0.035
        # minimum channel slope
        MIN_RIV_SLP = 0.00001

        # default pre processed and well prepared spatial data name

        # name of dem
        self.dem = "dem"
        # name of projected dem
        self.dem_projected = "dem_proj"
        # name of flow direction dataset for arcgis and grass
        self.fdr_grass = "fdr_grass"
        self.fdr_arcgis = "fdr_arcgis"
        # name of flow accumulation data
        self.acc = "acc"
        # name of lake dataset
        self.all_lakes = "all_lakes"
        self.all_selected_lake = "selected_lakes"
        # name of observation dataset
        self.observation_point = "obspoint"
        # name of flood plain manning's n datset
        self.flood_n = "flood_n"
        # name of landuse dataset
        self.landuse = "landuse"
        # name of bankfull width and depth dataset
        self.bankf_wd = "bk_wd"
        # name of mask
        self.mask = "MASK"
        # name of the result attribute table data type in grass format
        self.attribte_type = "attribute_type.csvt"

        # variable needs to be calcuated later

        # raster cellsize
        self.cellsize = -9.9999
        # input crs projection id
        self.crs_in = "#"
        # numboer of columns in processing extent
        self.ncols = -9999
        # number of rows in processing extent
        self.nrows = -9999
        # a list indicate potential wrong river reaches
        self.remove_str = []

    # first modulized methods
    def define_project_extent_method(
        self,
        mode,
        path_dem_in,
        outlet_pt=[-1, -1],
        path_extent_ply="#",
        buffer_distance=0.0,
        hybasin_ply="#",
        down_hybasin_id=-1,
        up_hybasin_id=-1,
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

        define_project_extent(
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            mode=mode,
            path_dem_in=path_dem_in,
            outlet_pt=outlet_pt,
            path_extent_ply=path_extent_ply,
            buffer_distance=buffer_distance,
            hybasin_ply=hybasin_ply,
            down_hybasin_id=down_hybasin_id,
            up_hybasin_id=up_hybasin_id,
            mask=self.mask,
            dem=self.dem,
        )
