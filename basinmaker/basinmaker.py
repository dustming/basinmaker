import os

class postproc:

    def __init__(
        self,
    ):
        self.qgispp = os.environ["QGIS_PREFIX_PATH"]
        print("asdf ")

    def generate_raven_model_inputs(
        self,
        path_final_hru_info="#",
        model_name="test",
        subbasingroup_nm_channel=["Allsubbasins"],
        subbasingroup_length_channel=[-1],
        subbasingroup_nm_lake=["AllLakesubbasins"],
        subbasingroup_area_lake=[-1],
        OutputFolder="#",
        aspect_from_gis = "grass",
    ):
        """Generate Raven input files.

        Function that used to generate Raven input files. All output will be stored in folder
        "<OutputFolder>/RavenInput".

        Parameters
        ----------
        OutputFolder                   : string
            is the folder that stores generated outputs
        path_final_hru_info     : string
            Path of the output HRUs shapefile from BasinMaker which includes all
            required parameters; Each row in the attribute table of this shapefile 
            represent a HRU.
        model_name      : string
           The Raven model base name. File name of the raven input will be
           Model_Name.xxx.
        subbasingroup_nm_channel       : List
            It is a list of names for subbasin groups, which are grouped based
            on channel length of each subbsin. Should at least has one name
        subbasingroup_length_channel   : List
            It is a list of float channel length thresthold in meter, to divide
            subbasin into different groups. for example, [1,10,20] will divide
            subbasins into four groups: 1)group 1 with channel length (0,1];
            2) group 2 with channel length (1,10];
            3) group 3 with channel length (10,20];
            and 4) group 4 with channel length (20,Max channel length].
        subbasingroup_nm_lake          : List
            It is a list of names for subbasin groups, which are grouped based
            on Lake area of each subbsin. Should at least has one name
        subbasingroup_area_lake        : List
            It is a list of float lake area thresthold in m2, to divide
            subbasin into different groups. for example, [1,10,20] will divide
            subbasins into four groups: 1) group 1 with lake area (0,1];
            2) group 2 with lake are (1,10],
            3) group 3 with lake are (10,20],
            4) group 4 with lake are (20,Max channel length].

        Returns
        -------

        Notes
        -------
        Following ouput files will be generated in "<OutputFolder>/RavenInput"
        
        | modelname.rvh              - contains subbasins and HRUs
        | Lakes.rvh                  - contains definition and parameters of lakes
        | channel_properties.rvp     - contains definition and parameters for channels

        Examples
        -------

        """

        from basinmaker.hymodin.raveninput import (
            GenerateRavenInput,
        )

        startyear=-1,
        endYear=-1,
        CA_HYDAT="#",
        warmup=0,
        template_folder="#",
        lake_as_gauge=False,
        writeobsrvt=False,
        downloadobsdata=False,
        forcing_input_file="#",


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
            aspect_from_gis = aspect_from_gis
        )


    def obtain_grids_polygon_from_netcdf_file_method(
        self,
        netcdf_path="#",
        output_folder="#",
        coor_x_nm="lon",
        coor_y_nm="lat",
        is_rotated_grid=1,
        r_coor_x_nm="rlon",
        r_coor_y_nm="rlat",
        spatial_ref="EPSG:4326",
        x_add=-360,
        y_add=0,
        gis_platform="qgis",
    ):


        from basinmaker.postprocessing.postprocessingfunctions import (
            obtain_grids_polygon_from_netcdf_file
        )

        obtain_grids_polygon_from_netcdf_file(
            netcdf_path=netcdf_path,
            output_folder=output_folder,
            coor_x_nm=coor_x_nm,
            coor_y_nm=coor_y_nm,
            is_rotated_grid=is_rotated_grid,
            r_coor_x_nm=r_coor_x_nm,
            r_coor_y_nm=r_coor_y_nm,
            spatial_ref=spatial_ref,
            x_add=x_add,
            y_add=y_add,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )


    def generate_area_weight_of_two_polygons_method(
        self,
        target_polygon_path="#",
        mapping_polygon_path="#",
        col_nm="HRU_ID",
        output_folder="#",
        gis_platform='qgis',
    ):
    
        from basinmaker.postprocessing.postprocessingfunctions import (
            generate_area_weight_of_two_polygons
        )

        generate_area_weight_of_two_polygons(
            target_polygon_path=target_polygon_path,
            mapping_polygon_path=mapping_polygon_path,
            col_nm=col_nm,
            output_folder=output_folder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )
        

    def simplify_routing_structure_by_filter_lakes(
        self,
        OutputFolder="#",
        Routing_Product_Folder = '#',
        Thres_Area_Conn_Lakes=-1,
        Thres_Area_Non_Conn_Lakes=-1,
        Selected_Lake_List_in=[],
        gis_platform="qgis",
    ):
        """Function that used to simplify the routing product by user
        provided lake area threstholds.

        Parameters
        ----------
        OutputFolder                   : string
            is the folder that stores generated outputs
        Routing_Product_Folder         : string
            is the folder where the input routing product is stored
        Thres_Area_Conn_Lakes          : float (optional)
            It is a lake area threshold for connected lakes in km2, 
            connected lake with lake area smaller than this value will be removed
        Thres_Area_Non_Conn_Lakes      : float (optional)
            It is a lake area threshold for non connected lakes in km2, 
            non connected lake with lake area smaller than this value will be removed
        Selected_Lake_List_in          : list
            A list of lake IDs from in the routing product. 
            Lakes with their lake ID in this list will be kept by the BasinMaker even 
            if their area smaller than the lake area threstholds
        gis_platform                   : string
            It is the parameter indicate which gis platform is used. It can be 
            either "qgis" or "arcgis".
            
        Returns
        -------
        
        Notes
        -----
        This function has no return values, the simplified routing product 
        will be saved in OutputFolder
        
        Examples
        -------


        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            simplify_routing_structure_by_filter_lakes_method,
        )
        Path_final_riv_ply = '#'
        Path_final_riv = '#'
        Path_Con_Lake_ply = '#'
        Path_NonCon_Lake_ply = '#'
        simplify_routing_structure_by_filter_lakes_method(
            Path_final_riv_ply=Path_final_riv_ply,
            Path_final_riv=Path_final_riv,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            Routing_Product_Folder = Routing_Product_Folder,
            Thres_Area_Conn_Lakes=Thres_Area_Conn_Lakes,
            Thres_Area_Non_Conn_Lakes=Thres_Area_Non_Conn_Lakes,
            Selected_Lake_List_in=Selected_Lake_List_in,
            OutputFolder=OutputFolder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )
    
        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = OutputFolder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )
        

    def simplify_routing_structure_by_drainage_area(
        self,
        OutputFolder="#",
        Routing_Product_Folder = '#',
        Drain_Area_Min=-1,
        gis_platform="qgis",
        qgis_prefix_path="#",
    ):
        """Function that used to simplify the routing product by
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

        Parameters
        ----------
        OutputFolder                   : string
            is the folder that stores generated outputs
        Routing_Product_Folder         : string
            is the folder where the input routing product is stored
        Drain_Area_Min                 : float (optional)
            It is a catchment drainage area thresthold, catchment with their 
            drainage area smaller than this thresthold will be removed.     
        gis_platform                   : string
            It is the parameter indicate which gis platform is used. It can be 
            either "qgis" or "arcgis".
            
        Returns
        -------
        
        Notes
        -----
        This function has no return values, the simplified routing product 
        will be saved in OutputFolder
        
        Examples
        -------


        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            simplify_routing_structure_by_drainage_area_method,
        )
        
        simplify_routing_structure_by_drainage_area_method(
            Routing_Product_Folder = Routing_Product_Folder,
            Area_Min=Drain_Area_Min,
            OutputFolder=OutputFolder,
            gis_platform=gis_platform,
            qgis_prefix_path=self.qgispp,
        )

        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = OutputFolder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )
        
    def select_part_of_routing_product(
        self,
        OutputFolder="#",
        Routing_Product_Folder = '#',
        mostdownids=[-1],
        mostupids=[-1],
        gis_platform="qgis",
    ):
        """Extract region of interest based on provided subbasin IDs
                
        Parameters
        ----------
        OutputFolder                   : string
            is the folder that stores generated outputs
        Routing_Product_Folder         : string
            is the folder where the input routing product is stored
        mostdownids                    : list
            A list of subbasin ID, the subbasin IDs in this list should 
            be the most downstream subbasin ID of each interested watershed. 
        mostupids                      : list 
            A list of subbasin ID, the subbasin IDs in this list should be 
            the most upstream subbasin ID of each interested watershed. 
            It is should be -1 when the entire watershed drainage to subbasin IDs 
            in mostdownid is needed. It should be some subbbasin ID, 
            when we needs an incomplete watershed, then the drainage area between 
            mostdownid and mostupid will be extracted.
        gis_platform                   : string
            It is the parameter indicate which gis platform is used. It can be 
            either "qgis" or "arcgis".
        Returns
        -------
        
        Notes
        -----
        This function has no return values, instead extected routing product 
        for region of interested will be saved in OutputFolder
        
        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            select_part_of_routing_product_method,
        )
        Path_Points="#",
        Gauge_NMS="#",
        Path_Catchment_Polygon="#",
        Path_River_Polyline="#",
        Path_Con_Lake_ply="#",
        Path_NonCon_Lake_ply="#",
        qgis_prefix_path="#",

        select_part_of_routing_product_method(
            Path_Points=Path_Points,
            Gauge_NMS=Gauge_NMS,
            OutputFolder=OutputFolder,
            mostdownid=mostdownid,
            mostupid=mostupid,
            Path_Catchment_Polygon=Path_Catchment_Polygon,
            Path_River_Polyline=Path_River_Polyline,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            qgis_prefix_path=qgis_prefix_path,
            Routing_Product_Folder = Routing_Product_Folder,
            gis_platform=gis_platform,
        )

        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = OutputFolder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )
        
    def generate_hrus(
        self,
        OutputFolder,
        Path_Subbasin_Ply,
        Landuse_info,
        Soil_info,
        Veg_info,
        Path_Connect_Lake_ply="#",
        Path_Non_Connect_Lake_ply="#",
        Path_Landuse_Ply="#",
        Path_Soil_Ply="#",
        Path_Veg_Ply="#",
        Path_Other_Ply_1="#",
        Path_Other_Ply_2="#",
        Inmportance_order = [],
        min_hru_area_pct_sub = 0.0,
        DEM="#",
        Project_crs="EPSG:3573",
        qgis_prefix_path="#",
        gis_platform="qgis",
    ):
        """Function that be used to overlay: subbasin polygon, lake polygon (optional)
        , Land use polygon (optional), soil type polygon(optional),
        vegetation polygon (optional), and two other user defined polygons
        (optional).Non-lake HRU polygons in a subbasin is defined by an unique
        combination of all user provided datasets.

        Parameters
        ----------
        OutputFolder                   : string
            is the folder that stores generated outputs
        Path_Subbasin_Ply                 : string
            It is the path of the subbasin polygon, which is generated by
            BasinMaker. 
        Landuse_info                      : string
            Path to a csv file that contains landuse information, including
            following attributes:

            | Landuse_ID (integer) -- the landuse ID in the landuse polygon,-1 for lake
            | LAND_USE_C (string) -- the landuse class name for each landuse type
        Soil_info                        : string
            Path to a csv file that contains soil information, including
            following attributes:
            
            | Soil_ID (integer) -- the soil ID  in the soil polygon,-1 for lake
            | SOIL_PROF (string) -- the soil profile name for each soil profile type            
        Veg_info                         : string
            Path to a csv file that contains vegetation information, including
            following attributes:
            
            | Veg_ID (integer) -- the vegetation ID  in the vegetation polygon,-1 for lake
            | VEG_C (string) -- the vegetation class name for each vegetation Type
        Path_Connect_Lake_ply            : string (Optional)
            Path to the connected lake's polygon
        Path_Non_Connect_Lake_ply        : string (Optional)
            Path to the non connected lake's polygon
        Path_Landuse_Ply                 : string (Optional)
            Path to the landuse polygon. when Path_Landuse_Ply is not
            provided. The Landuse ID in Landuse_info should be
            1: land, -1: lake
        Path_Soil_Ply                    : string (Optional)
            Path to the soil polygon. when soil polygon is not
            provided. The Soil ID in Soil_info should be the same
            as Landuse ID.
        Path_Veg_Ply                     : string (Optional)
            Path to the vegetation polygon. when Veg polygon is not
            provided. The Veg ID in Veg_info should be the same
            as Landuse ID.
        Path_Other_Ply_1                 : string (Optional)
            Path to the other polygon that will be used to define HRU,
            such as elevation band, or aspect.
        Path_Other_Ply_2                 : string (Optional)
            Path to the other polygon that will be used to define HRU,
            such as elevation band, or aspect.
        DEM                              : string (optional)
            the path to a raster elevation dataset, that will be used to
            calcuate average apspect, elevation and slope within each HRU.
            if no data is provided, basin average value will be used for
            each HRU.
        Project_crs                      : string
            the EPSG code of a projected coodinate system that will be used to
            calcuate HRU area and slope.
        Returns
        -------
        
        Notes
        -----
        This function has no return values, instead extected routing product 
        for region of interested will be saved in OutputFolder
        
        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            generate_hrus_method,
        )
        Sub_Lake_ID="HyLakeId",
        Sub_ID="SubId",
        Lake_Id="Hylak_id",
        Landuse_ID="Landuse_ID",
        Soil_ID="Soil_ID",
        Other_Ply_ID_1="O_ID_1",
        Veg_ID="Veg_ID",
        Other_Ply_ID_2="O_ID_2",

        generate_hrus_method(
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
            Inmportance_order = Inmportance_order,
            min_hru_area_pct_sub = min_hru_area_pct_sub,
            Project_crs=Project_crs,
            OutputFolder=OutputFolder,
            qgis_prefix_path=qgis_prefix_path,
            gis_platform = gis_platform,
        )


class dlidem:

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
        path_output_folder="#",
        path_working_folder="#",
    ):

        # define drived values
        # create folders
        self.path_output_folder = path_output_folder
        self.path_working_folder = path_working_folder

        #        os.makedirs(self.path_output_folder, exist_ok=True)
        os.makedirs(self.path_working_folder, exist_ok=True)

        # obtain qgis prefix path
        self.qgispp = os.environ["QGIS_PREFIX_PATH"]
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
            "sub_reg_str_r": "sub_reg_str_r",
            "sub_reg_str_v": "sub_reg_str_v",
            "sub_reg_nfdr_grass": "sub_reg_nfdr_grass",
            "sub_reg_nfdr_arcgis": "sub_reg_nfdr_arcgis",
            "sub_reg_acc": "sub_reg_acc",
            "sub_reg_dem": "sub_reg_dem",
            "problem_seg": "problem_seg",
            "lake_outflow_pourpoints": "lake_outflow_pourpoints",
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
        gis_platform="qgis",
    ):

        """Function that used to define project processing spatial extent (PSE).
        BasinMaker will not process grids or features outside the processing spatial extent.

        Parameters
        ----------
        mode                              : string (required)
            It is a string indicate which method to define project processing
            spatial extent
            
            | 'using_dem'            : the extent of input dem will be used
            | 'using_hybasin'        : the extent will be defined using hydrobasin product
            | 'using_outlet_pt'      : the extent will be defined with provided outlet point
            | 'using_provided_ply'   : the extent will be defined by provided polygon    
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

        Returns
        -------

        Notes
        -------
        Outputs are following files in grass database loacated in 
        os.path.join(grassdb, grass_location) 
        
        MASK.*        : raster  
            it is a mask raster stored in grass database, which indicate the PSE.
            The grass database is located at 
        dem.*         : raster 
            it is a dem raster stored in grass database, which is has the same extent 
            with MASK. The grass database is located at os.path.join(grassdb, grass_location)

        Examples
        -------
        """

        from basinmaker.extent.projectextent import define_project_extent
        up_hybasin_id = -1
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
        fdr_path="#",
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
        subreg_fdr_path="#",
        subreg_acc_path="#",
        subreg_str_r_path="#",
        subreg_str_v_path="#",
    ):
        """Function that used to Generate a routing structure with using user 
        provied flow accumulation thresthold without considering lake.

        Parameters
        ----------
        accthresold       : float
            It is the flow accumulation thresthold, used to determine
            subbsains and river network. Increasing of accthresold will
            increase the size of generated subbasins, reduce the number
            subbasins and reduce the number of generated stream segments
        fdr_path          : string (optional)
            It is a string indicate path of flow direction dataset
        mode              : string (required)
            It is a string indicate which dataset will be used to delineate
            watershed.
            | 'usingdem'             : dem is used for delineation
            | 'usingfdr'             : flow direction data is used for delineation
            | 'usingsubreg'          : predefined subregion inputs is used for delineation
        max_memroy        : integer
            It is the maximum memeory that allow to be used.

        Returns
        -------

        Notes
        -------
        Outputs are following files in the grass database

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


        Examples
        -------
        """
        from basinmaker.delineationnolake.watdelineationwithoutlake import (
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
        only_included_lake_at_river_interction = False,
        path_obsfile_in="#",
        obs_attributes=[],
        path_sub_reg_outlets_v="#",
        path_sub_reg_lake_r="#",
        path_sub_reg_lake_bd_r="#",
        search_radius=100,
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
    ):
        """ Update the subbasin delineation result generated by lake inflow and 
        outflow points and observation gauges will be added into as a new subbasin
        outlet. Result from this tool is not the final delineation result. Because 
        some lake may cover several subbasins. these subbasin are not megered into 
        one subbasin yet


        Parameters
        ----------
        path_lakefile_in                   : string (optional)
            It is a string to indicate the full path of the polygon 
            shapefile that include lake data
        lake_attributes                   : list (optional)
            the columns names in the path_lakefile_in that indicate following 
            items has to be included
            
            | 1) column name for the unique Id of each lake within the lake polygon shpfile
            | 2) column name for type of the lake should be integer
            | 3) column name for the volume of the lake in km3
            | 4) column name for the average depth of the lake in m
            | 5) column name for the area of the lake in km2
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
        obs_attributes                    : list (optional)
            the columns names in path_obsfile_in that indicate following items has 
            to be included
            
            | 1) column name for the unique Id of each observation point
            | 2) column name for the unique name of each observation point
            | 3) column name for the drainage area of each observation point in km3
            | 4) column name for the source of the observation point: 'CA' for observation in canada;
                 'US' for observation in US, or any other names                      
                  
        Returns
        -------

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


        Examples
        -------

        """
        from basinmaker.addlakeandobs.addlakeandobsintowatershed import (
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
            only_included_lake_at_river_interction = only_included_lake_at_river_interction,
            threshold_non_con_lake=threshold_non_con_lake,
            search_radius=search_radius,
            path_sub_reg_lake_r=path_sub_reg_lake_r,
            path_sub_reg_lake_bd_r=path_sub_reg_lake_bd_r,
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
        projection="EPSG:3573",
        k_in=-1,
        c_in=-1,
        gis_platform="qgis",
        lake_attributes=[],
        obs_attributes=[],
        outlet_obs_id=-1,
        path_sub_reg_outlets_v="#",
        output_folder="#",
    ):
        """Calculate hydrological paramters for each subbasin. 

        Parameters
        ----------
        path_bkfwidthdepth             : string
            It is a string to indicate the full path of the
            polyline shapefile that having bankfull width and
            depth data
        bkfwd_attributes               : list
            the columns names that indicate following items has to be included
            
            | 1) column name for the Bankfull width in m;
            | 2) column name for the Bankfull depth in m;
            | 3) column name for the annual mean discharge in m3/s;
        path_landuse                   : string
            It is a string to indicate the full path of the landuse data.
            It will be used to estimate the floodplain roughness
            coefficient. Should have the same projection with the DEM data
            in ".tif" format.
        path_landuse_info              : string
            It is a string to indicate the full path of the table in '.csv'
            format.The table describe the floodplain roughness coefficient
            correspond to a given landuse type. The table should have two
            columns:
             
            | RasterV: is the landuse value in the landuse raster for each land use type
            | MannV: is the roughness coefficient value for each landuse type.
        projection                     : string
            It is a string EPSG code to indicate a projected coordinate system,
            which wiil be used to calcuate area, slope and aspect.            
        lake_attributes                   : list (optional)
            the columns names in the path_lakefile_in that indicate following 
            items has to be included
            
            | 1) column name for the unique Id of each lake within the lake polygon shpfile
            | 2) column name for type of the lake should be integer
            | 3) column name for the volume of the lake in km3
            | 4) column name for the average depth of the lake in m
            | 5) column name for the area of the lake in km2
        obs_attributes                    : list (optional)
            the columns names in path_obsfile_in that indicate following items has 
            to be included
            
            | 1) column name for the unique Id of each observation point
            | 2) column name for the unique name of each observation point
            | 3) column name for the drainage area of each observation point in km3
            | 4) column name for the source of the observation point: 'CA' for observation in canada;
                 'US' for observation in US, or any other names       
        OutputFolder                  : string
            The path to a folder to save outputs

        Returns
        -------

        Notes
        -------
        Five vector files will be generated in the output folder. these files
        can be used to define final routing structure by "combine_catchments_covered_by_the_same_lake"
        or be used as input for other postprocessing tools. All files
        are stored at

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


        Examples
        -------

        """

        from basinmaker.addattributes.addattributestocatchments import add_attributes_to_catchments

        add_attributes_to_catchments(
            input_geo_names=self.geofilenames,
            path_bkfwidthdepth=path_bkfwidthdepth,
            bkfwd_attributes=bkfwd_attributes,
            path_landuse=path_landuse,
            path_landuse_info=path_landuse_info,
            projection=projection,
            k_in=k_in,
            c_in=k_in,
            out_cat_name=self.geofilenames["catchment_without_merging_lakes"],
            out_riv_name=self.geofilenames["river_without_merging_lakes"],
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
            obs_attributes=obs_attributes,
            lake_attributes=lake_attributes,
            outlet_obs_id=outlet_obs_id,
            path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            output_folder=output_folder,
        )

    def divide_domain_into_sub_regions_method(
        self,
        path_lakefile_in,
        lake_attributes,
        path_bkfwidthdepth="#",
        bkfwd_attributes="#",
        Min_Num_Domain=9,
        Max_Num_Domain=13,
        Initaial_Acc=5000,
        Delta_Acc=1000,
        CheckLakeArea=1,
        fdr_path="#",
        Acc_Thresthold_stream=500,
        max_memory=2048 * 3,
        Out_Sub_Reg_Folder="#",
        gis_platform="qgis",
    ):

        from basinmaker.subreg.defsubreg import (
            divide_domain_into_sub_regions,
        )

        divide_domain_into_sub_regions(
            input_geo_names=self.geofilenames,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            path_lakefile_in=path_lakefile_in,
            lake_attributes=lake_attributes,
            path_bkfwidthdepth=path_bkfwidthdepth,
            bkfwd_attributes=bkfwd_attributes,
            Min_Num_Domain=Min_Num_Domain,
            Max_Num_Domain=Max_Num_Domain,
            Initaial_Acc=Initaial_Acc,
            Delta_Acc=Delta_Acc,
            CheckLakeArea=CheckLakeArea,
            fdr_path=fdr_path,
            Acc_Thresthold_stream=Acc_Thresthold_stream,
            max_memory=max_memory,
            Out_Sub_Reg_Folder=Out_Sub_Reg_Folder,
            sub_reg_str_r=self.geofilenames["sub_reg_str_r"],
            sub_reg_str_v=self.geofilenames["sub_reg_str_v"],
            sub_reg_nfdr_grass=self.geofilenames["sub_reg_nfdr_grass"],
            sub_reg_nfdr_arcgis=self.geofilenames["sub_reg_nfdr_arcgis"],
            sub_reg_acc=self.geofilenames["sub_reg_acc"],
            sub_reg_dem=self.geofilenames["sub_reg_dem"],
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
        start_sub_id = 0,
        k =1,
        c = 1,
    ):

        from basinmaker.subreg.defsubreg import (
            combine_sub_region,
        )

        combine_sub_region(
            path_sub_region_info=path_sub_region_info,
            sub_region_outputfolder=sub_region_outputfolder,
            outputfolder=outputfolder,
            is_final_result=is_final_result,
            qgis_prefix_path=self.qgispp,
            path_subregion_inlet=path_subregion_inlet,
            start_sub_id = start_sub_id,
            k = k,
            c = c
        )

    def combine_catchments_covered_by_the_same_lake(
        self,
        Routing_Product_Folder='#',
        gis_platform="qgis",
    ):
        """Generate the final lake river routing structure by merging subbasin
        polygons that are covered by the same lake.

        Parameters
        ----------
        Routing_Product_Folder         : string
            is the folder where the input routing product is stored
        gis_platform                   : string
            It is the parameter indicate which gis platform is used. It can be 
            either "qgis" or "arcgis".
        Returns
        -------
        
        Notes
        -----
        This function has no return values, instead the generated routing product 
        will be saved in OutputFolder
        
        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = Routing_Product_Folder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )
        
