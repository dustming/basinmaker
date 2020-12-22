import os


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

        os.makedirs(self.path_output_folder, exist_ok=True)
        os.makedirs(self.path_working_folder, exist_ok=True)

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
        self.mask = "MASK"
        # name of projected dem
        self.dem_projected = "dem_proj"
        # name of flow direction dataset for arcgis and grass
        self.fdr_grass = "fdr_grass"
        self.fdr_arcgis = "fdr_arcgis"
        # name of flow accumulation data
        self.acc = "acc"
        self.str_r = "str_r"
        self.str_v = "str_v"
        self.cat_no_lake = "cat_no_lake"
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
        gis_platform="qgis",
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

        from extent.projectextent import define_project_extent

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
            gis_platform="qgis",
        )

    def preprocessing_inputs_method(
        self,
        path_lakefile_in="#",
        lake_attributes=[],
        gis_platform="qgis",
    ):

        """Preprocessing input dataset

        Function that used to project and clip input dataset such as
        DEM, Land use, Lake polygon etc with defined processing extent
        by function Generatmaskregion. And then it will rasterize these
        vector files. Generated raster files will has the
        same extent and resolution as the "dem" generated by Generatmaskregion
        and will be stored in a grass database located at
        os.path.join(self.tempFolder, 'grassdata_toolbox').

        Parameters
        ----------

        path_widep_in                      : string (optional)
            It is a string to indicate the full path of the
            polyline shapefile that having bankfull width and
            depth data. At least following three columns should
            exist in the shapefile attribute table:  1) WIDTH,
            is the Bankfull width in m; 2) DEPTH, is the Bankfull
            depth in m; 3) Q_Mean, is the annual mean discharge
            in m3/s. If it is not provided, -9999 will be used
            for bankfull width and depth in the generated routing
            structure
        path_lakefile_in                   : string (optional)
            It is a string to indicate the full path of the polygon
            shapefile that include lake data. Follow attributes needs
            to be include in the lake polygon shpfile: 1) Hylak_id,
            the unique Id of each lake within the lake polygon shpfile;
            2) Lake_type, type of the lake should be integer; 3) Vol_total,
            the volume of the lake in km3; 4) Depth_avg, the average depth
            of the lake in m; 5) Lake_area, the area of the lake in km2.
            If it is not provided, -9999 will be used for lake related
            parameters in the generated routing structure
        lake_attributes                    : list
            It is a list of column name in the lake polygon attribute table
            lake_attributes[0] should be the lake Id column
            lake_attributes[1] should be the Lake type column
            lake_attributes[2] should be the Lake area column
            lake_attributes[3] should be the lake volume column
            lake_attributes[4] should be the lake depth column
        path_landuse_in                    : string (optional)
            It is a string to indicate the full path of the landuse data.
            It will be used to estimate the floodplain roughness
            coefficient. Should have the same projection with the DEM data
            in ".tif" format.
        path_landuseinfo_in                : string (optional)
            It is a string to indicate the full path of the table in '.csv'
            format.The table describe the floodplain roughness coefficient
            correspond to a given landuse type. The table should have two
            columns: RasterV and MannV. RasterV is the landuse value in the
            landuse raster for each land use type and the MannV is the
            roughness coefficient value for each landuse type.
        path_flow_direction_in             : string (optional)
            It is a string to indicate the full path of the input flow direction
            dataset in '.tif' format

        Notes
        -------

        Raster and vector files that are generated by this function and
        will be used in next step are list in following. All files are
        stored at a grass database in os.path.join(self.tempFolder,
        'grassdata_toolbox')

        alllake              : raster
            it is a raster represent all lakes within the processing extent.
        Lake_Bound           : raster
            it is a raster represent the lake boundary grids
        acc_grass            : raster
            it is the flow accumulation raster generated by 'r.watershed'
        width                : raster
            it is the rasterized bankfull width depth polyline shapefile
            using column "WIDTH"
        depth                : raster
            it is the rasterized bankfull width depth polyline shapefile
            using column "DEPTH"
        qmean                : raster
            it is the rasterized bankfull width depth polyline shapefile
            using column "Q_Mean"
        up_area              : raster
            it is the rasterized bankfull width depth polyline shapefile
            using column "UP_AREA"
        SubId_WidDep         : raster
            it is the rasterized bankfull width depth polyline shapefile
            using column "HYBAS_ID"
        landuse              : raster
            it is the landuse raster.

        Returns:
        -------
           None

        Examples
        -------
        """

        from preprocessing.preprocessinginput import preprocessing_inputs

        preprocessing_inputs(
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            mask=self.mask,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform="qgis",
        )

    def watershed_delineation_without_lake_method(
        self,
        acc_thresold=-1,
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
    ):
        from delineationnolake.watdelineationwithoutlake import (
            watershed_delineation_without_lake,
        )

        watershed_delineation_without_lake(
            mode,
            mask=self.mask,
            dem=self.dem,
            acc_thresold=acc_thresold,
            fdr_path="#",
            subreg_fdr_path="#",
            subreg_acc_path="#",
            subreg_str_r_path="#",
            subreg_str_v_path="#",
            fdr_arcgis=self.fdr_arcgis,
            fdr_grass=self.fdr_grass,
            str_r=self.str_r,
            str_v=self.str_v,
            acc=self.acc,
            cat_no_lake=self.cat_no_lake,
            max_memroy=max_memroy,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform="qgis",
        )


    def watershed_delineation_with_lakes_method(
        self,
        path_lakefile_in="#",
        lake_attributes=[],
        path_obsfile_in = '#',
        obs_attributes=[],
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
    ):
        from addlakeandobs.addlakeandobsintowatershed import (
            add_lakes_and_obs_into_existing_watershed_delineation,
        )

        add_lakes_and_obs_into_existing_watershed_delineation(
            mode,
            mask=self.mask,
            dem=self.dem,
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            path_obsfile_in = path_lakefile_in,
            obs_attributes=obs_attributes,
            fdr_arcgis=self.fdr_arcgis,
            fdr_grass=self.fdr_grass,
            str_r=self.str_r,
            str_v=self.str_v,
            acc=self.acc,
            cat_no_lake=self.cat_no_lake,
            max_memroy=max_memroy,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform="qgis",
        )
        
        