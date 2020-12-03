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
        os.makedirs(self.path_output_folder,exists_ok = False)
        os.makedirs(self.path_working_folder,exists_ok = False)

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
        self.grass_location_geo_temp = "temporary_location1"
        self.grass_location_geo_temp1 = "temporary_location2"
        self.grass_location_pro = "project_crs_location"

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
