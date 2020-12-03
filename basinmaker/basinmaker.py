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
        path_dem_in="#",
        path_wid_dep_in="#",
        path_lake_in="#",
        path_dir_in="#",
        path_landuse_in="#",
        path_landuse_n_in="#",
        path_obspoint="#",
        path_subregion_info_folder="#",
        debug=False,
        mode="using_dem",
    ):

        # initialize various input dataset folder
        self.path_dem_in = path_dem_in
        self.path_wid_dep_in = path_wid_dep_in
        self.path_lake_in = path_lake_in
        self.path_dir_in = path_dir_in
        self.path_landuse_in = path_landuse_in
        self.path_landuse_n_in = path_landuse_n_in
        self.path_output_folder = path_outputfolder
        self.path_working_folder = path_working_folder
        self.path_subregion_info_folder = path_subregion_info_folder
        self.debug = debug
        self.mode = mode

        # validate the inputs
        if self.mode == "using_dem":
            assert (
                self.path_dem_in != "#"
            ), "need provide dem data when mode = using_dem"

        if self.mode == "using_dir":
            assert (
                self.path_dir_in != "#"
            ), "need provide flow direction data when mode = using_dir"

        if self.mode == "using_subregion":
            assert (
                self.path_subregion_info_folder != "#"
            ), "need provide subregion data folder when mode = using_subregion"

        # define drived values
        # create folders
        if not os.path.exists(self.path_output_folder):
            os.makedirs(self.path_output_folder)

        if not os.path.exists(self.path_working_folder):
            os.makedirs(self.path_working_folder)

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

        # define subregion data path if mode is using_subregion
        if self.mode == "using_subregion":

            self.path_subregion_outlets_r = os.path.join(
                self.path_subregion_info_folder, "Sub_Reg_Outlets_point_r.pack"
            )
            self.path_subregion_outlets_v = os.path.join(
                self.path_subregion_info_folder, "Sub_Reg_Outlets_point_v.pack"
            )
            self.path_subregion_grass_dir = os.path.join(
                self.path_subregion_info_folder, "dir_grass.pack"
            )
            self.path_subregion_arcgis_dir = os.path.join(
                self.path_subregion_info_folder, "dir_Arcgis.pack"
            )
            self.path_subregion_grass_acc = os.path.join(
                self.path_subregion_info_folder, "acc_grass.pack"
            )
            self.path_subregion_grass_str_r = os.path.join(
                self.path_subregion_info_folder, "Sub_Reg_str_grass_r.pack"
            )
            self.path_subregion_grass_str_v = os.path.join(
                self.path_subregion_info_folder, "Sub_Reg_str_grass_v.pack"
            )
            self.path_subregion_dem = os.path.join(
                self.path_subregion_info_folder, "Sub_Reg_dem.pack"
            )

        # define constants

        # field names in the attribute table
        self.fieldnamelist = [
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
        self.default_chn = 0.035
        # minimum channel slope
        self.min_riv_slp = 0.00001

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
        self.Remove_Str = []
