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
        self.geofilenames = {
            "dem": "dem",
            "mask": "MASK",
            "dem_proj": "dem_proj",
            "fdr_grass": "fdr_grass",
            "fdr_arcgis": "fdr_arcgis",
            "acc": "acc",
            "str_r": "str_r",
            "str_v": "str_v",
            "cat_no_lake": "cat_no_lake",
            "all_lakes": "all_lakes",
            "selected_lakes": "selected_lakes",
            "connect_lake": "connect_lake",
            "nonconnect_lake": "nonconnect_lake",
            "str_connected_lake": "str_connected_lake",
            "sl_nonconnect_lake": "sl_nonconnect_lake",
            "sl_connected_lake": "sl_connected_lake",
            "str_sl_connected_lake": "str_sl_connected_lake",
            "obspoint": "obspoint",
            "flood_n": "flood_n",
            "bk_wd": "bk_wd",
            "lake_boundary": "lake_boundary",
            "nfdr_arcgis": "nfdr_arcgis",
            "nfdr_grass": "nfdr_grass",
            "catchment_without_merging_lakes": "catchment_without_merging_lakes",
            "river_without_merging_lakes": "river_without_merging_lakes",
            "cat_use_default_acc": "cat_use_default_acc",
            "snapped_obs_points": "snapped_obs_points",
        }

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

        mask                   : raster
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
            mask=self.geofilenames["mask"],
            dem=self.geofilenames["dem"],
            gis_platform=gis_platform,
        )

    def watershed_delineation_without_lake_method(
        self,
        acc_thresold=-1,
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
        subreg_fdr_path="#",
        subreg_acc_path="#",
        subreg_str_r_path="#",
        subreg_str_v_path="#",
        fdr_path="#",
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
                                     delineation
        max_memroy        : integer
            It is the maximum memeory that allow to be used.
        gis_platform      : string
            It is a string indicate with gis platform is used:
            'qgis'                : the basinmaker is running within QGIS
            'arcgis'              : the basinmaker is running within ArcGIS
        subreg_fdr_path   : string
            It is a string indicate path of subregion flow direction dataset
        subreg_acc_path   : string
            It is a string indicate path of subregion flow accumlation dataset
        subreg_str_r_path : string
            It is a string indicate path of subregion stream raster flow
            direction dataset
        subreg_str_v_path : string
            It is a string indicate path of subregion stream vector datasets
        fdr_path          : string
            It is a string indicate path of flow direction dataset

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
        from delineationnolake.watdelineationwithoutlake import (
            watershed_delineation_without_lake,
        )

        watershed_delineation_without_lake(
            mode=mode,
            input_geo_names=self.geofilenames,
            acc_thresold=acc_thresold,
            fdr_path=fdr_path,
            subreg_fdr_path=subreg_fdr_path,
            subreg_acc_path=subreg_acc_path,
            subreg_str_r_path=subreg_str_r_path,
            subreg_str_v_path=subreg_str_v_path,
            fdr_arcgis=self.geofilenames["fdr_arcgis"],
            fdr_grass=self.geofilenames["fdr_grass"],
            str_r=self.geofilenames["str_r"],
            str_v=self.geofilenames["str_v"],
            acc=self.geofilenames["acc"],
            cat_no_lake=self.geofilenames["cat_no_lake"],
            max_memroy=max_memroy,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

    def watershed_delineation_add_lake_control_points(
        self,
        path_lakefile_in="#",
        lake_attributes=[],
        threshold_con_lake=0,
        threshold_non_con_lake=0,
        path_obsfile_in="#",
        obs_attributes=[],
        path_sub_reg_outlets_v="#",
        search_radius=100,
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
    ):
        """Add lake inflow and outflow points as new subabsin outlet

        Update the subbasin delineation result generated by
        "watershed_delineation_without_lake_method". Lake's inflow and outflow
        points will be added into subbasin delineation result as a new subbasin
        outlet. Result from this tool is not the final delineation result.
        Because some lake may cover several subbasins. these subbasin
        are not megered into one subbasin yet.

        Parameters
        ----------
        path_lakefile_in                   : string (optional)
            It is a string to indicate the full path of the polygon
            shapefile that include lake data.
        lake_attributes                    ：list
            the columns names that indicate following items has to be included
            1) column name for the unique Id of each lake within the lake polygon shpfile;
            2) column name for type of the lake should be integer;
            3) column name for the volume of the lake in km3;
            4) column name for the average depth of the lake in m;
            5) column name for the area of the lake in km2.
        threshold_con_lake                 : float
            It is a lake area thresthold for connected lakes in km2
            Connected Lake with lake area below this value will be
            not added into the subbasin delineation
        threshold_non_con_lake             : float
            It is a lake area thresthold for non-connected lakes in km2
            Non connected Lake with lake area below this value will be
            not added into the subbasin delineation
        path_obsfile_in
            It is a string to indicate the full path of the point
            shapefile that indicate observation gauges
        obs_attributes                    ：list
            the columns names that indicate following items has to be included
            1) column name for the unique Id of each observation point;
            2) column name for the unique name of each observation point;
            3) column name for the drainage area of each observation point in km3;
            4) column name for the source of the observation point:
                'CA' for observation in canada;
                'US' for observation in US;
        search_radius       : integer
            It is the search ratio in number of grids to snap observation
            point into the river network.
        path_sub_reg_outlets_v  : path
            It is the full path of the subregion outlet vector file

        Notes
        -------
        Raster and vector files that are generated by this function and
        will be used by next step are list as following. All files are
        stored at a grass database

        selected_lakes                    : raster
            it is a raster represent all lakes that are selected by two lake
            area threstholds
        sl_nonconnect_lake       : raster
            it is a raster represent all non connected lakes that are selected
            by lake area threstholds
        sl_connected_lake           : raster
            it is a raster represent allconnected lakes that are selected
            by lake area threstholds
        river_without_merging_lakes                         : raster/vector
            it is the updated river segment for each subbasin
        catchment_without_merging_lakes                     : raster/vector
            it is a raster represent updated subbasins after adding lake inflow
            and outflow points as new subbasin outlet.
        snapped_obs_points                                  : raster/vector
            it is a name of the FIG file represent successfully sanpped
            observation points
        Returns:
        -------
           None

        Examples
        -------

        """
        from addlakeandobs.addlakeandobsintowatershed import (
            add_lakes_and_obs_into_existing_watershed_delineation,
        )

        add_lakes_and_obs_into_existing_watershed_delineation(
            input_geo_names=self.geofilenames,
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            path_obsfile_in=path_obsfile_in,
            obs_attributes=obs_attributes,
            path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            threshold_con_lake=threshold_con_lake,
            threshold_non_con_lake=threshold_non_con_lake,
            search_radius=search_radius,
            sl_connected_lake=self.geofilenames["sl_connected_lake"],
            sl_non_connected_lake=self.geofilenames["sl_nonconnect_lake"],
            sl_lakes=self.geofilenames["selected_lakes"],
            catchment_without_merging_lakes=self.geofilenames[
                "catchment_without_merging_lakes"
            ],
            river_without_merging_lakes=self.geofilenames[
                "river_without_merging_lakes"
            ],
            snapped_obs_points=self.geofilenames["snapped_obs_points"],
            cat_use_default_acc=self.geofilenames["cat_use_default_acc"],
            nfdr_arcgis=self.geofilenames["nfdr_arcgis"],
            nfdr_grass=self.geofilenames["nfdr_grass"],
            max_memroy=max_memroy,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

    def add_attributes_to_catchments_method(
        self,
        path_bkfwidthdepth="#",
        bkfwd_attributes=[],
        path_landuse="#",
        path_landuse_info="#",
        gis_platform="qgis",
        lake_attributes = [],
        obs_attributes=[],
        outlet_obs_id=1,
        path_sub_reg_outlets_v="#",
        output_folder="#",
    ):
        from addattributes.addattributestocatchments import add_attributes_to_catchments

        add_attributes_to_catchments(
            input_geo_names=self.geofilenames,
            path_bkfwidthdepth=path_bkfwidthdepth,
            bkfwd_attributes=bkfwd_attributes,
            path_landuse=path_landuse,
            path_landuse_info=path_landuse_info,
            out_cat_name=self.geofilenames["catchment_without_merging_lakes"],
            out_riv_name=self.geofilenames["river_without_merging_lakes"],
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
            obs_attributes=obs_attributes,
            lake_attributes = lake_attributes,
            outlet_obs_id=outlet_obs_id,
            path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            output_folder=output_folder,
        )

    def combine_catchments_covered_by_the_same_lake_method(
        self,
        OutputFolder="#",
        Path_final_rivply="#",
        Path_final_riv="#",
        gis_platform="qgis",
    ):
        from postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake,
        )

        combine_catchments_covered_by_the_same_lake(
            OutputFolder=OutputFolder,
            Path_final_rivply=Path_final_rivply,
            Path_final_riv=Path_final_riv,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

    def simplify_routing_structure_by_filter_lakes_method(
        self,
        Path_final_riv_ply="#",
        Path_final_riv="#",
        Path_Con_Lake_ply="#",
        Path_NonCon_Lake_ply="#",
        Thres_Area_Conn_Lakes=-1,
        Thres_Area_Non_Conn_Lakes=-1,
        Selection_Method="ByArea",
        Selected_Lake_List_in=[],
        OutputFolder="#",
        gis_platform="qgis",
    ):
        from postprocessing.postprocessingfunctions import (
            simplify_routing_structure_by_filter_lakes,
        )

        simplify_routing_structure_by_filter_lakes(
            Path_final_riv_ply=Path_final_riv_ply,
            Path_final_riv=Path_final_riv,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
            Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
            Selection_Method=Selection_Method,
            Selected_Lake_List_in=Selected_Lake_List_in,
            OutputFolder=OutputFolder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

    def simplify_routing_structure_by_drainage_area_method(
        self,
        Path_final_riv_ply="#",
        Path_final_riv="#",
        Path_Con_Lake_ply="#",
        Path_NonCon_Lake_ply="#",
        Area_Min=-1,
        OutputFolder="#",
        gis_platform="qgis",
        qgis_prefix_path="#",
    ):
        from postprocessing.postprocessingfunctions import (
            simplify_routing_structure_by_drainage_area,
        )

        simplify_routing_structure_by_drainage_area(
            Path_final_riv_ply=Path_final_riv_ply,
            Path_final_riv=Path_final_riv,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            Area_Min=Area_Min,
            OutputFolder=OutputFolder,
            gis_platform=gis_platform,
            qgis_prefix_path=self.qgispp,
        )
