import os


class basinmaker:

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
        path_output_folder = '#',
        path_working_folder = '#',
    ):

        # define drived values
        # create folders
        self.path_output_folder = path_output_folder
        self.path_working_folder = path_working_folder

#        os.makedirs(self.path_output_folder, exist_ok=True)
        os.makedirs(self.path_working_folder, exist_ok=True)

        # obtain qgis prefix path
        self.qgispp = os.environ["QGISPrefixPath"]
        # obtain basinmaker path
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
            "sub_reg_str_r":"sub_reg_str_r",
            "sub_reg_str_v":"sub_reg_str_v",
            "sub_reg_nfdr_grass":"sub_reg_nfdr_grass",
            "sub_reg_nfdr_arcgis":"sub_reg_nfdr_arcgis",
            "sub_reg_acc":"sub_reg_acc",
            "sub_reg_dem":"sub_reg_dem",
            "problem_seg":"problem_seg",
            "lake_outflow_pourpoints":"lake_outflow_pourpoints",
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
            it is a name of the point gis file represent successfully sanpped
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
        projection = "EPSG:3573",
        k_in=-1,
        c_in=-1,
        gis_platform="qgis",
        lake_attributes = [],
        obs_attributes=[],
        outlet_obs_id=-1,
        path_sub_reg_outlets_v="#",
        output_folder="#",
    ):
        """Calculate hydrological paramters

        Calculate hydrological paramters for each subbasin generated by
        "AutomatedWatershedsandLakesFilterToolset". The result generaed
        by this tool can be used as inputs for Define_Final_Catchment
        and other post processing tools

        Parameters
        ----------
        path_bkfwidthdepth             : string
            It is a string to indicate the full path of the
            polyline shapefile that having bankfull width and
            depth data        
        bkfwd_attributes               : 
            the columns names that indicate following items has to be included
            1) column name for the Bankfull width in m;
            2) column name for the Bankfull depth in m;
            3) column name for the annual mean discharge in m3/s; 
        path_landuse                   : string 
            It is a string to indicate the full path of the landuse data.
            It will be used to estimate the floodplain roughness
            coefficient. Should have the same projection with the DEM data
            in ".tif" format.        
        path_landuse_info              : string 
            It is a string to indicate the full path of the table in '.csv'
            format.The table describe the floodplain roughness coefficient
            correspond to a given landuse type. The table should have two
            columns: RasterV and MannV. RasterV is the landuse value in the
            landuse raster for each land use type and the MannV is the
            roughness coefficient value for each landuse type.        
        gis_platform                   : string
            It is a string indicate with gis platform is used:
            'qgis'                : the basinmaker is running within QGIS
            'arcgis'              : the basinmaker is running within ArcGIS
        lake_attributes                : list 
            the columns names that indicate following items has to be included
            1) column name for the unique Id of each lake within the lake polygon shpfile;
            2) column name for type of the lake should be integer;
            3) column name for the volume of the lake in km3;
            4) column name for the average depth of the lake in m;
            5) column name for the area of the lake in km2.        
        obs_attributes                 : list 
            the columns names that indicate following items has to be included
            1) column name for the unique Id of each observation point;
            2) column name for the unique name of each observation point;
            3) column name for the drainage area of each observation point in km3;
            4) column name for the source of the observation point:
                'CA' for observation in canada;
                'US' for observation in US;   
        outlet_obs_id                  : int 
            It is one 'Obs_ID' in the provided observation gauge
            shapefile. If it is larger than zero. Subbasins that
            do not drainage to this gauge will be removed from
            delineation result.                   
        projection                     : string
            It is a string indicate a projected coordinate system,
            which wiil be used to calcuate area, slope and aspect.
        output_folder                  : string
            The path to a folder to save outputs
        path_sub_reg_outlets_v         : string 
        
        Notes
        -------
        Five vector files will be generated by this function. these files
        can be used to define final routing structure by "Define_Final_Catchment"
        or be used as input for other postprocessing tools. All files
        are stored at self.OutputFolder

        catchment_without_merging_lakes.shp             : shapefile
            It is the subbasin polygon before merging lakes catchments and
            need to be processed before used.
        river_without_merging_lakes.shp                 : shapefile
            It is the subbasin river segment before merging lakes catchments and
            need to be processed before used.
        sl_connected_lake.shp                           : shapefile
            It is the connected lake polygon. Connected lakes are lakes that
            are connected by  Path_final_riv.
        sl_non_connected_lake.shp                       : shapefile
            It is the  non connected lake polygon. Connected lakes are lakes
            that are not connected by Path_final_cat_riv or Path_final_riv.
        obs_gauges                                      : shapefile
            It is the point shapefile that represent the observation gauge
            after snap to river network.
            
        Returns:
        -------
           None

        Examples
        -------

        """
        
        from addattributes.addattributestocatchments import add_attributes_to_catchments

        add_attributes_to_catchments(
            input_geo_names=self.geofilenames,
            path_bkfwidthdepth=path_bkfwidthdepth,
            bkfwd_attributes=bkfwd_attributes,
            path_landuse=path_landuse,
            path_landuse_info=path_landuse_info,
            projection = projection,
            k_in=k_in,
            c_in=k_in,
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
        """Define final lake river routing structure

        Generate the final lake river routing structure by merging subbasin
        polygons that are covered by the same lake.
        The input are the catchment polygons and river segements
        before merging for lakes. The input files can be output of
        any of following functions:
        SelectLakes, Select_Routing_product_based_SubId,
        Customize_Routing_Topology,RoutingNetworkTopologyUpdateToolset_riv
        The result is the final catchment polygon that ready to be used for
        hydrological modeling

        Parameters
        ----------
        OutputFolder                   : string
            Folder name that stores generated extracted routing product
        Path_final_riv_ply             : string
            Path to the catchment polygon which is the routing product
            before merging lakes catchments and need to be processed before
            used. It is the input for simplify the routing product based
            on lake area or drianage area.
            routing product and can be directly used.
        Path_final_riv                 : string
            Path to the river polyline which is the routing product
            before merging lakes catchments and need to be processed before
            used. It is the input for simplify the routing product based
            on lake area or drianage area.

        Notes
        -------
        This function has no return values, instead will generate following
        files. They are catchment polygons and river polylines that can be
        used for hydrological modeling.
        os.path.join(OutputFolder,'finalcat_info.shp')
        os.path.join(OutputFolder,'finalcat_info_riv.shp')

        Returns:
        -------
        None

        Examples
        -------

        """
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
        """Simplify the routing product by lake area

        Function that used to simplify the routing product by user
        provided lake area thresthold.
        The input catchment polygons is the routing product before
        merging for lakes. It is provided with the routing product.
        The result is the simplified catchment polygons. But
        result from this fuction still not merging catchment
        covering by the same lake. Thus, The result generated
        from this tools need further processed by
        Define_Final_Catchment, or can be further processed by
        Customize_Routing_Topology

        Parameters
        ----------

        Path_final_riv_ply             : string
            Path to the catchment polygon which is the routing product
            before merging lakes catchments and need to be processed before
            used. It is the input for simplify the routing product based
            on lake area or drianage area.
        Path_final_riv                 : string
            Path to the river polyline which is the routing product
            before merging lakes catchments and need to be processed before
            used. It is the input for simplify the routing product based
            on lake area or drianage area.
        Path_Con_Lake_ply              : string
            Path to a connected lake polygon. Connected lakes are lakes that
            are connected by Path_final_riv.
        Path_NonCon_Lake_ply           : string
            Path to a non connected lake polygon. Connected lakes are lakes
            that are not connected by Path_final_riv.
        Thres_Area_Conn_Lakes          : float (optional)
            It is the lake area threshold for connated lakes, in km2
        Thres_Area_Non_Conn_Lakes      : float (optional)
            It is the lake area threshold for non connated lakes, in km2
        Selection_Method               : string
            It is a string indicate lake selection methods
            "ByArea" means lake in the routing product will be selected based
            on two lake area thresthold Thres_Area_Conn_Lakes and
            Thres_Area_Non_Conn_Lakes
            "ByLakelist" means lake in the routing product will be selected
            based on user provided hydrolake id, in Selected_Lake_List_in
        Selected_Lake_List_in          : list
            A list of lake ids that will be keeped in the routing product.
            Lakes not in the list will be removed from routing product.
        OutputFolder                   : string
            Folder name that stores generated extracted routing product

        Notes
        -------
        This function has no return values, instead will generate following
        files. The output tpye will be the same as inputs, but the routing
        network will be simplified by removing lakes.

        os.path.join(OutputFolder,os.path.basename(Path_final_riv_ply))
        os.path.join(OutputFolder,os.path.basename(Path_final_riv))
        os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply))
        os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply))

        Returns:
        -------
        None

        Examples
        -------

        """
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
        """Simplify the routing product by drainage area

        Function that used to simplify the routing product by
        using user provided minimum subbasin drainage area.
        The input catchment polygons are routing product before
        merging for lakes. It is provided with routing product.
        The result is the simplified catchment polygons. But
        result from this fuction still not merging catchment
        covering by the same lake. Thus, The result generated
        from this tools need further processed by
        Define_Final_Catchment, or can be further processed by
        SelectLakes

        Parameters
        ----------

        Path_final_riv_ply             : string
            Path to the catchment polygon which is the routing product
            before merging lakes catchments and need to be processed before
            used. It is the input for simplify the routing product based
            on lake area or drianage area.
        Path_final_riv                 : string
            Path to the river polyline which is the routing product
            before merging lakes catchments and need to be processed before
            used. It is the input for simplify the routing product based
            on lake area or drianage area.
        Path_Con_Lake_ply              : string
            Path to a connected lake polygon. Connected lakes are lakes that
            are connected by Path_final_riv.
        Path_NonCon_Lake_ply           : string
            Path to a non connected lake polygon. Connected lakes are lakes
            that are not connected by Path_final_riv.
        Area_Min                       : float
            The minimum drainage area of each catchment in km2
        OutputFolder                   : string
            Folder name that stores generated simplified routing product

        Notes
        -------
        This function has no return values, instead will generate following
        files. The output tpye will be the same as inputs, but the routing
        network will be simplified by increase subbasin size, reduce
        number of subbasins and number of river segments.

        os.path.join(OutputFolder,os.path.basename(Path_final_riv_ply))
        os.path.join(OutputFolder,os.path.basename(Path_final_riv))
        os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply))
        os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply))

        Returns:
        -------
        None

        Examples
        -------

        """
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

    def select_part_of_routing_product_methods(
        self,
        Path_Points = '#',
        Gauge_NMS = '#',
        OutputFolder = '#',
        mostdownid = -1,
        mostupid = -1,
        Path_Catchment_Polygon="#",
        Path_River_Polyline="#",
        Path_Con_Lake_ply="#",
        Path_NonCon_Lake_ply="#",
        qgis_prefix_path="#",
        gis_platform = 'qgis',
    ):
        """Extract region of interest based on provided pourpoints or gauge name

        Parameters
        ----------
        Path_Points                    : string (Optional)
            It is the path of the point shapefile. If the point shapefile is
            provided. The function will return subids of those catchment
            polygons that includes these point in the point shapefile
        Gauge_NMS                      : list
            Name of the streamflow gauges, such as ['09PC019'], if the gauge
            name is provided, the subbasin ID that contain this gauge will be
            returned
        OutputFolder                   : string
            Folder path that stores extracted routing product
        Path_Catchment_Polygon         : string
            Path to the catchment polygon
        Path_River_Polyline            : string (optional)
            Path to the river polyline
        Path_Con_Lake_ply              : string (optional)
            Path to a connected lake polygon. Connected lakes are lakes that
            are connected by Path_final_cat_riv or Path_final_riv.
        Path_NonCon_Lake_ply           : string (optional)
            Path to a non connected lake polygon. Connected lakes are lakes
            that are not connected by Path_final_cat_riv or Path_final_riv.


        Notes
        -------
        This function has no return values, instead following fiels will be
        generated. The output files have are same as inputs expect the extent
        are different.

        os.path.join(OutputFolder,os.path.basename(Path_Catchment_Polygon))
        os.path.join(OutputFolder,os.path.basename(Path_River_Polyline))
        os.path.join(OutputFolder,os.path.basename(Path_Con_Lake_ply))
        os.path.join(OutputFolder,os.path.basename(Path_NonCon_Lake_ply))

        Returns:
        -------
        None

        Examples
        -------

        """
        from postprocessing.postprocessingfunctions import (
            select_part_of_routing_product,
        )
        select_part_of_routing_product(
            Path_Points=Path_Points,
            Gauge_NMS=Gauge_NMS,
            OutputFolder=OutputFolder,
            mostdownid = mostdownid,
            mostupid = mostupid,
            Path_Catchment_Polygon=Path_Catchment_Polygon,
            Path_River_Polyline=Path_River_Polyline,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            qgis_prefix_path=qgis_prefix_path,
            gis_platform = gis_platform,
        )
        
    def generate_hrus_methods(
        self,
        Path_Subbasin_Ply,
        Landuse_info,
        Soil_info,
        Veg_info,
        Sub_Lake_ID="HyLakeId",
        Sub_ID="SubId",
        Path_Connect_Lake_ply="#",
        Path_Non_Connect_Lake_ply="#",
        Lake_Id="Hylak_id",
        Path_Landuse_Ply="#",
        Landuse_ID="Landuse_ID",
        Path_Soil_Ply="#",
        Soil_ID="Soil_ID",
        Path_Veg_Ply="#",
        Veg_ID="Veg_ID",
        Path_Other_Ply_1="#",
        Other_Ply_ID_1="O_ID_1",
        Path_Other_Ply_2="#",
        Other_Ply_ID_2="O_ID_2",
        DEM="#",
        Project_crs="EPSG:3573",
        OutputFolder="#",
        qgis_prefix_path='#',
        gis_platform='qgis',
    ):
        """Generate HRU polygons and their attributes needed by hydrological model

        Function that be used to overlay: subbasin polygon, lake polygon (optional)
        , Land use polygon (optional), soil type polygon(optional),
        vegetation polygon (optional), and two other user defined polygons
        (optional).
        Non-lake HRU polygons in a subbasin is defined by an unique
        combination of all user provided datasets.
        A lake HRU polygon is defined the same as the provided lake polygon.
        All value of landuse and Veg polygon covered by lake will
        be changed to 1, indicating it is a covered by lake.
        All value of the soil polygon covered by the lake will be change to
        the soil id of the polygon covered by the lake with largest area.

        Parameters
        ----------
        Path_Subbasin_Ply                 : string
            It is the path of the subbasin polygon, which is generated by
            toolbox. if not generated by toolbox, the attribute table should
            including following attribute.
            ##############Subbasin related attributes###########################
            SubID           - integer, The subbasin Id
            DowSubId        - integer, The downstream subbasin ID of this
                                       subbasin
            IsLake          - integer, If the subbasin is a lake / reservior
                                       subbasin. 1 yes, <0, no
            IsObs           - integer, If the subbasin contains a observation
                                       gauge. 1 yes, < 0 no.
            RivLength       - float,   The length of the river in current
                                       subbasin in m
            RivSlope        - float,   The slope of the river path in
                                       current subbasin, in m/m
            FloodP_n        - float,   Flood plain manning's coefficient, in -
            Ch_n            - float,   main channel manning's coefficient, in -
            BkfWidth        - float,   the bankfull width of the main channel
                                       in m
            BkfDepth        - float,   the bankfull depth of the main channel
                                       in m
            HyLakeId        - integer, the lake id
            LakeVol         - float,   the Volume of the lake in km3
            LakeDepth       - float,   the average depth of the lake m
            LakeArea        - float,   the area of the lake in m2
        Landuse_info                      : string
            Path to a csv file that contains landuse information, including
            following attributes:
            Landuse_ID (can be any string)  - integer, the landuse ID in the
                                                       landuse polygon
            LAND_USE_C                      - string,  the landuse class name
                                                       for each landuse Type
        Soil_info                        : string
            Path to a csv file that contains soil information, including
            following attributes:
            Soil_ID (can be any string)     - integer, the Soil ID in the
                                                       soil polygon
            SOIL_PROF                       - string,  the Soil profile name
                                                       for each soil type
        Veg_info                         : string
            Path to a csv file that contains vegetation information, including
            following attributes:
            Veg_ID (can be any string)      - integer, the vegetation ID in the
                                                       vegetation polygon
            VEG_C                           - string,  the vegetation class name
                                                       for each vegetation Type
        Sub_Lake_ID                      : string (optional)
            The column name of the lake id in the subbasin polygon
        Sub_ID                           : string (optional)
            The column name of the subbasin id in the subbasin polygon
        Path_Connect_Lake_ply            : string (Optional)
            Path to the connected lake's polygon
        Path_Non_Connect_Lake_ply        : string (Optional)
            Path to the non connected lake's polygon
        Lake_Id                          : string (Optional)
            The the column name in lake polygon indicate the lake ID.
        Path_Landuse_Ply                 : string (Optional)
            Path to the landuse polygon. when Path_Landuse_Ply is not
            provided. The Landuse ID in Landuse_info should be
            1: land, -1: lake
        Landuse_ID                       : string (Optional)
            the the column name in landuse polygon and Landuse_info csv
            indicate the landuse ID. when Path_Landuse_Ply is not
            provided. The Landuse ID should be
            1: land, -1: lake.
        Path_Soil_Ply                    : string (Optional)
            Path to the soil polygon. when soil polygon is not
            provided. The Soil ID in Soil_info should be the same
            as Landuse ID.
        Soil_ID                          : string (Optional)
            the the column name in soil polygon and soil_info csv
            indicate the soil ID. when soil polygon is not
            provided. The Soil ID in Soil_info should be the same
            as Landuse ID.
        Path_Veg_Ply                     : string (Optional)
            Path to the vegetation polygon. when Veg polygon is not
            provided. The Veg ID in Veg_info should be the same
            as Landuse ID.
        Veg_ID                           : string (Optional)
            the the column name in vegetation polygon and veg_info csv
            indicate the vegetation ID. when Veg polygon is not
            provided. The Veg ID in Veg_info should be the same
            as Landuse ID.
        Path_Other_Ply_1                 : string (Optional)
            Path to the other polygon that will be used to define HRU,
            such as elevation band, or aspect.
        Other_Ply_ID_1                   : string (Optional)
            the the column name in Other_Ply_1 polygon
            indicate the landuse ID.
        Path_Other_Ply_2                 : string (Optional)
            Path to the other polygon that will be used to define HRU,
            such as elevation band, or aspect.
        Other_Ply_ID_2                   : string (Optional)
            the the column name in Other_Ply_2 polygon
            indicate the landuse ID.
        DEM                              : string (optional)
            the path to a raster elevation dataset, that will be used to
            calcuate average apspect, elevation and slope within each HRU.
            if no data is provided, basin average value will be used for
            each HRU.
        Project_crs                      : string
            the EPSG code of a projected coodinate system that will be used to
            calcuate HRU area and slope.
        OutputFolder                     : string
            The path to the folder that will save output HRU polygon.

        Notes
        -------
        Following ouput files will be generated in "<OutputFolder>/"
        'finalcat_hru_info.shp'              - HRU polygon and it's attributes


        Returns:
        -------
           None

        Examples
        -------
        >>> from ToolboxClass import LRRT
        >>> import pandas as pd
        >>> DataFolder = "C:/Path_to_foldr_of_example_dataset_provided_in_Github_wiki/"
        >>> RTtool=LRRT()
        >>> RTtool.GenerateHRUS(OutputFolder = DataFolder,
                               Path_Subbasin_Ply = os.path.join(DataFolder,"finalcat_info.shp"),
                               Path_Connect_Lake_ply = os.path.join(DataFolder,'Con_Lake_Ply.shp'),
                               Path_Non_Connect_Lake_ply = os.path.join(DataFolder,'Non_Con_Lake_Ply.shp'),
                               Path_Landuse_Ply = os.path.join(DataFolder,'modislanduse_exp_lg_pre.shp'),
                               Landuse_ID = 'gridcode',
                               Path_Soil_Ply = os.path.join(DataFolder,'ca_all_slc_v3r2_exp_lg.shp'),
                               Soil_ID = 'POLY_ID',
                               Landuse_info=os.path.join(DataFolder,'landuse_info.csv'),
                               Soil_info=os.path.join(DataFolder,'soil_info.csv'),
                               Veg_info=os.path.join(DataFolder,'veg_info.csv'),
                               DEM = os.path.join(DataFolder,'na_dem_15s_1.tif')
                               )

        """
        from postprocessing.postprocessingfunctions import (
            generate_hrus,
        )
        generate_hrus(
            Path_Subbasin_Ply=Path_Subbasin_Ply,
            Landuse_info=Landuse_info,
            Soil_info=Soil_info,
            Veg_info=Veg_info,
            Sub_Lake_ID=Sub_Lake_ID,
            Sub_ID=Sub_ID,
            Path_Connect_Lake_ply=Path_Connect_Lake_ply,
            Path_Non_Connect_Lake_ply=Path_Non_Connect_Lake_ply,
            Lake_Id=Lake_Id,
            Path_Landuse_Ply=Path_Landuse_Ply,
            Landuse_ID=Landuse_ID,
            Path_Soil_Ply=Path_Soil_Ply,
            Soil_ID=Soil_ID,
            Path_Veg_Ply=Path_Veg_Ply,
            Veg_ID=Veg_ID,
            Path_Other_Ply_1=Path_Other_Ply_1,
            Other_Ply_ID_1=Other_Ply_ID_1,
            Path_Other_Ply_2=Path_Other_Ply_2,
            Other_Ply_ID_2=Other_Ply_ID_2,
            DEM=DEM,
            Project_crs=Project_crs,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
        )
        
    def divide_domain_into_sub_regions_method(
        self,
        path_lakefile_in,
        lake_attributes,
        path_bkfwidthdepth = '#',
        bkfwd_attributes = '#',
        Min_Num_Domain=9,
        Max_Num_Domain=13,
        Initaial_Acc=5000,
        Delta_Acc=1000,
        CheckLakeArea=1,
        fdr_path = '#',
        Acc_Thresthold_stream=500,
        max_memory=2048*3,
        Out_Sub_Reg_Folder="#",
        gis_platform='qgis',
    ):

        from subreg.defsubreg import (
            divide_domain_into_sub_regions,
        )

        divide_domain_into_sub_regions(
            input_geo_names=self.geofilenames,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            path_bkfwidthdepth = path_bkfwidthdepth,
            bkfwd_attributes = bkfwd_attributes,
            Min_Num_Domain=Min_Num_Domain,
            Max_Num_Domain=Max_Num_Domain,
            Initaial_Acc=Initaial_Acc,
            Delta_Acc=Delta_Acc,
            CheckLakeArea=CheckLakeArea,
            fdr_path = fdr_path,
            Acc_Thresthold_stream=Acc_Thresthold_stream,
            max_memory=max_memory,
            Out_Sub_Reg_Folder=Out_Sub_Reg_Folder,
            sub_reg_str_r = self.geofilenames["sub_reg_str_r"],
            sub_reg_str_v = self.geofilenames["sub_reg_str_v"],
            sub_reg_nfdr_grass = self.geofilenames["sub_reg_nfdr_grass"],
            sub_reg_nfdr_arcgis = self.geofilenames["sub_reg_nfdr_arcgis"],
            sub_reg_acc = self.geofilenames["sub_reg_acc"],
            sub_reg_dem = self.geofilenames["sub_reg_dem"],
            gis_platform=gis_platform,
        )

    def combine_sub_region_results_method(
        self,
        path_sub_region_info,
        sub_region_outputfolder,
        outputfolder,
        is_final_result,
        path_subregion_inlet,
        gis_platform="qgis",
    ):

        from subreg.defsubreg import (
            combine_sub_region,
        )

        combine_sub_region(
            path_sub_region_info = path_sub_region_info,
            sub_region_outputfolder = sub_region_outputfolder,
            outputfolder = outputfolder,
            is_final_result = is_final_result,
            qgis_prefix_path = self.qgispp,
            path_subregion_inlet = path_subregion_inlet,
        )



        
    def generate_raven_model_inputs_method(
        self,
        path_final_hru_info="#",
        startyear=-1,
        endYear=-1,
        CA_HYDAT="#",
        warmup=0,
        template_folder="#",
        lake_as_gauge=False,
        writeobsrvt=False,
        downloadobsdata=True,
        model_name="test",
        subbasingroup_nm_channel=["Allsubbasins"],
        subbasingroup_length_channel=[-1],
        subbasingroup_nm_lake=["AllLakesubbasins"],
        subbasingroup_area_lake=[-1],
        outputfolder="#",
        forcing_input_file="#",
    ):
        """Generate Raven input files.

        Function that used to generate Raven input files. All output will be stored in folder
        "<OutputFolder>/RavenInput".

        Parameters
        ----------

        path_final_hru_info     : string
            Path of the output shapefile from routing toolbox which includes
            all required parameters; Each row in the attribute table of this
            shapefile represent a HRU. Different HRU in the same subbasin has
            the same subbasin related attribute values.

            The shapefile should at least contains following columns
            ##############Subbasin related attributes###########################
            SubId           - integer, The subbasin Id
            DowSubId        - integer, The downstream subbasin ID of this
                                       subbasin
            IsLake          - integer, If the subbasin is a lake / reservior
                                       subbasin. 1 yes, <0, no
            IsObs           - integer, If the subbasin contains a observation
                                       gauge. 1 yes, < 0 no.
            RivLength       - float,   The length of the river in current
                                       subbasin in m
            RivSlope        - float,   The slope of the river path in
                                       current subbasin, in m/m
            FloodP_n        - float,   Flood plain manning's coefficient, in -
            Ch_n            - float,   main channel manning's coefficient, in -
            BkfWidth        - float,   the bankfull width of the main channel
                                       in m
            BkfDepth        - float,   the bankfull depth of the main channel
                                       in m
            HyLakeId        - integer, the lake id
            LakeVol         - float,   the Volume of the lake in km3
            LakeDepth       - float,   the average depth of the lake m
            LakeArea        - float,   the area of the lake in m2
            ############## HRU related attributes    ###########################
            HRU_S_mean      - float,   the slope of the HRU in degree
            HRU_A_mean      - float,   the aspect of the HRU in degree
            HRU_E_mean      - float,   the mean elevation of the HRU in m
            HRU_ID          - integer, the id of the HRU
            HRU_Area        - integer, the area of the HRU in m2
            HRU_IsLake      - integer, the 1 the HRU is a lake hru, -1 not
            LAND_USE_C      - string,  the landuse class name for this HRU, the
                                       name will be used in Raven rvh, and rvp
                                       file
            VEG_C           - string,  the Vegetation class name for this HRU, the
                                       name will be used in Raven rvh, and rvp
                                       file
            SOIL_PROF       - string,  the soil profile name for this HRU, the
                                       name will be used in Raven rvh, and rvp
                                       file
            HRU_CenX        - float,  the centroid coordinates for HRU in x
                                       dimension
            HRU_CenY        - float,  the centroid coordinates for HRU in y
                                       dimension
        lake_as_gauge   : Bool
            If "True", all lake subbasins will labeled as gauged
            subbasin such that Raven will export lake balance for
            this lake. If "False", lake subbasin will not be labeled
            as gauge subbasin.
        CA_HYDAT        : string,  optional
            path and filename of downloaded
            external database containing streamflow observations,
            e.g. HYDAT for Canada ("Hydat.sqlite3").
        startyear       : integer, optional
            Start year of simulation. Used to
            read streamflow observations from external databases.
        endYear         : integer, optional
            End year of simulation. Used to
            read streamflow observations from external databases.
        warmup          : integer, optional
            The warmup time (in years) used after
            startyear. Values in output file "obs/xxx.rvt" containing
            observations will be set to NoData value "-1.2345" from
            model start year to end of WarmUp year.
        template_folder : string, optional
            Input that is used to copy raven template files. It is a
            folder name containing raven template files. All
            files from that folder will be copied (unchanged)
            to the "<OutputFolder>/RavenInput".
        writeobsrvt     : Bool, optional
            Input that used to indicate if the observation data file needs
            to be generated.
        downloadobsdata : Bool, optional
            Input that used to indicate if the observation data will be Download
            from usgs website or read from hydat database for streamflow Gauge
            in US or Canada,respectively. If this parameter is False,
            while WriteObsrvt is True. The program will write the observation data
            file with "-1.2345" for each observation gauges.
        model_name      : string
           The Raven model base name. File name of the raven input will be
           Model_Name.xxx.
        subbasingroup_nm_channel       : List
            It is a list of names for subbasin groups, which are grouped based
            on channel length of each subbsin. Should at least has one name
        subbasingroup_length_channel   : List
            It is a list of float channel length thresthold in meter, to divide
            subbasin into different groups. for example, [1,10,20] will divide
            subbasins into four groups, group 1 with channel length (0,1];
            group 2 with channel length (1,10],
            group 3 with channel length (10,20],
            group 4 with channel length (20,Max channel length].
        subbasingroup_nm_lake          : List
            It is a list of names for subbasin groups, which are grouped based
            on Lake area of each subbsin. Should at least has one name
        subbasingroup_area_lake        : List
            It is a list of float lake area thresthold in m2, to divide
            subbasin into different groups. for example, [1,10,20] will divide
            subbasins into four groups, group 1 with lake area (0,1];
            group 2 with lake are (1,10],
            group 3 with lake are (10,20],
            group 4 with lake are (20,Max channel length].
        outputfolder                   : string
            Folder name that stores generated Raven input files. The raven
            input file will be generated in "<OutputFolder>/RavenInput"

        Notes
        -------
        Following ouput files will be generated in "<OutputFolder>/RavenInput"
        modelname.rvh              - contains subbasins and HRUs
        Lakes.rvh                  - contains definition and parameters of lakes
        channel_properties.rvp     - contains definition and parameters for channels
        xxx.rvt                    - (optional) streamflow observation for each gauge
                                     in shapefile database will be automatically
                                     generagted in folder "OutputFolder/RavenInput/obs/".
        obsinfo.csv                - information file generated reporting drainage area
                                     difference between observed in shapefile and
                                     standard database as well as number of missing
                                     values for each gauge

        Returns:
        -------
           None

        Examples
        -------

        """

        from hymodin.raveninput import (
            GenerateRavenInput,
        )

        GenerateRavenInput(
            Path_final_hru_info=path_final_hru_info,
            lenThres=1,
            iscalmanningn=-1,
            Startyear=startyear,
            EndYear=endYear,
            CA_HYDAT=CA_HYDAT,
            WarmUp=warmup,
            Template_Folder=template_folder,
            Lake_As_Gauge=lake_as_gauge,
            WriteObsrvt=writeobsrvt,
            DownLoadObsData=downloadobsdata,
            Model_Name=model_name,
            Old_Product=False,
            SubBasinGroup_NM_Channel=subbasingroup_nm_channel,
            SubBasinGroup_Length_Channel=subbasingroup_length_channel,
            SubBasinGroup_NM_Lake=subbasingroup_nm_lake,
            SubBasinGroup_Area_Lake=subbasingroup_area_lake,
            OutputFolder=outputfolder,
            Forcing_Input_File=forcing_input_file,
        )


        