import os

class postprocess:

    def __init__(
        self,
    ):
        if os.getenv("QGIS_PREFIX_PATH"):
            self.qgispp = os.environ["QGIS_PREFIX_PATH"]
        else:
            self.qgispp = '#'

    def Generate_Raven_Model_Inputs(
        self,
        path_output_folder,
        path_hru_polygon="#",
        aspect_from_gis = "grass",
        model_name="test",
        subbasingroup_names_channel=["Allsubbasins"],
        subbasingroup_length_channel=[-1],
        subbasingroup_name_lake=["AllLakesubbasins"],
        subbasingroup_area_lake=[-1],
        lake_out_flow_method = 'broad_crest',
        model_time_step = 1,
    ):
        """This function is to generate Raven input files. A subfolder 'RavenInput'
        in path_output_folder is created to sotre resultant outputs.

        Parameters
        ----------
        path_output_folder                   : string
            is the folder path that stores generated Raven input fiels
        path_hru_polygon     : string
            is path of the output HRU shapefile from BasinMaker which includes all
            required hydrologic attributes and parameters; Each row in the attribute
            table of this shapefile represent a HRU.
        aspect_from_gis                  : string
            is a string indicating how aspect of each HRU in path_hru_polygon is
            caculated. "grass" represent the aspect is calcuated by GRASS GIS;
            "arcgis" represent the aspect is calcuated by ArcGIS.
        model_name      : string (optional)
            is the The Raven model base name.
        subbasingroup_names_channel       : list (optional)
            is a list of names for subbasin groups, which are grouped based
            on channel length of each subbsin. Should at least has one name
        subbasingroup_length_channel   : list (optional)
            is a list of channel length thresthold in meter, that divides
            subbasin into different groups. For example, [1,10,20] divides
            subbasins into four groups:

            | 1) group 1 with channel length range from (0,1];
            | 2) group 2 with channel length range from (1,10];
            | 3) group 3 with channel length range from (10,20];
            | 4) group 4 with channel length range from (20,Max channel length].
        subbasingroup_name_lake          : list (optional)
            is a list of names for subbasin groups, which are grouped based
            on area of lakes in lake subbasins. Should at least has one name
        subbasingroup_area_lake        : list (optional)
            is a list of lake area thresthold in m2, that divides
            subbasin into different groups. For example, [1,10,20] divides
            subbasins into four groups:

            | 1) group 1 with lake area range from (0,1];
            | 2) group 2 with lake are range from (1,10],
            | 3) group 3 with lake are range from (10,20],
            | 4) group 4 with lake are range from (20,max_lake_area].
        Returns
        -------

        Notes
        -------
        The following ouput files will be generated in "<path_output_folder>/RavenInput"

        | modelname.rvh              - contains subbasins and HRUs
        | Lakes.rvh                  - contains definition and parameters of lakes
        | channel_properties.rvp     - contains definition and parameters for channels

        Examples
        -------

        """

        from basinmaker.hymodin.raveninput import (
            GenerateRavenInput,
        )

        startyear=-1
        endYear=-1
        CA_HYDAT="#"
        warmup=0
        template_folder="#"
        lake_as_gauge=False
        writeobsrvt=False
        downloadobsdata=False
        forcing_input_file="#"


        GenerateRavenInput(
            Path_final_hru_info=path_hru_polygon,
            lenThres=0.2,
            iscalmanningn=1,
            Startyear=startyear,
            EndYear=endYear,
            CA_HYDAT=CA_HYDAT,
            WarmUp=warmup,
            time_step = model_time_step,
            Template_Folder=template_folder,
            Lake_As_Gauge=lake_as_gauge,
            WriteObsrvt=writeobsrvt,
            DownLoadObsData=downloadobsdata,
            Model_Name=model_name,
            Old_Product=False,
            SubBasinGroup_NM_Channel=subbasingroup_names_channel ,
            SubBasinGroup_Length_Channel=subbasingroup_length_channel ,
            SubBasinGroup_NM_Lake=subbasingroup_name_lake ,
            SubBasinGroup_Area_Lake=subbasingroup_area_lake ,
            OutputFolder=path_output_folder,
            Forcing_Input_File=forcing_input_file,
            aspect_from_gis = aspect_from_gis,
            lake_out_flow_method = lake_out_flow_method,
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


    def Remove_Small_Lakes(
        self,
        path_output_folder="#",
        routing_product_folder  = '#',
        connected_lake_area_thresthold=-1,
        non_connected_lake_area_thresthold =-1,
        selected_lake_ids=[],
        area_thresthold = 10*30*30/1000/1000,
        gis_platform="qgis",
    ):
        """This function is to simplify the hydrologic routing network by removing
        lakes.

        Parameters
        ----------
        path_output_folder                   : string
            is the folder path that stores generated outputs
        routing_product_folder         : string
            is the folder path where the input hydrologic routing network is stored
        gis_platform                   : string
            is the parameter indicating which gis platform is used. It can be
            either "qgis" or "arcgis".
        connected_lake_area_thresthold                 : float (optional)
            is a lake area thresthold for connected lakes in km2.
            Connected lake with lake area below which lake will be removed
        non_connected_lake_area_thresthold             : float (optional)
            is a lake area thresthold for non-connected lakes in km2
            Non connected lake with lake area below below which lake will be removed
        selected_lake_ids          : list (optional)
            A list of lake IDs from in the hydrologic routing network (Column 'HyLakeId').
            Lakes with their lake ID in this list will be kept by the BasinMaker even
            if their area smaller than the lake area threstholds
        area_thresthold             : float (optional)
            This parameter sets a subbasin area threshold in square kilometers (km²). 
            Subbasins with an area below this threshold will be merged with downstream 
            or upstream subbasins, depending on their location within the river network. 
            This can be useful for simplifying the river network and reducing computational
            requirements. The parameter is optional and has a default value of 0.009 km2. 
            If the parameter is set to 0, no subbasin merging will be performed.

        Returns
        -------

        Notes
        -----
        This function has no return values, The simplified hydrological routing
        network will be generated in the path_output_folder including following
        files:

        finalcat_info.shp                     : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            subbasin polygons which respect the lake inflow and outflow routing
            structures. This layer contains all the necessary information for
            hydrologic routing through the lake-river network.
        finalcat_info_riv.shp                 : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            river network polylines in the routing network.
        catchment_without_merging_lakes.shp             : shapefile
            The GIS layer containing subbasin polygons of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network subbasin
            polygons covered by the same lake are not merged into one lake subbasin
            yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        river_without_merging_lakes.shp                 : shapefile
            The GIS layer containing river polylines of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network, the
            river polylines covered by the same lake are not merged into one river
            segment yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        sl_connected_lake.shp                           : shapefile
            the GIS layer containing the lake polygons of lakes that are connected
            by the river_without_merging_lakes.shp
        sl_non_connected_lake.shp                       : shapefile
            the GIS layer containing the lake polygons of lakes that are not connected
            by the river_without_merging_lakes.shp
        poi                                            : shapefile
            It is the point shapefile that represent the point of interest
            after snap to river network.

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
            Routing_Product_Folder = routing_product_folder ,
            Thres_Area_Conn_Lakes=connected_lake_area_thresthold ,
            Thres_Area_Non_Conn_Lakes=non_connected_lake_area_thresthold ,
            Selected_Lake_List_in=selected_lake_ids ,
            OutputFolder=path_output_folder ,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = path_output_folder,
            qgis_prefix_path=self.qgispp,
            area_thresthold = area_thresthold,
            gis_platform=gis_platform,
        )


    def Decrease_River_Network_Resolution(
        self,
        path_output_folder="#",
        routing_product_folder = '#',
        gis_platform="qgis",
        minimum_subbasin_drainage_area=-1,
        area_thresthold = 10*30*30/1000/1000,
    ):
        """ This function is to simplify the hydrologic routing network by
        removing subbasins/river reaches with their drainage area below user provided
        drainage area thresthold.

        Parameters
        ----------
        path_output_folder                   : string
            is the folder path that stores generated outputs
        routing_product_folder         : string
            is the folder path where the input hydrologic routing network is stored
        gis_platform                   : string
            is the parameter indicating which gis platform is used. It can be
            either "qgis" or "arcgis".
        minimum_subbasin_drainage_area                 : float
            is a subbasin drainage area thresthold, subbasin with their
            drainage area smaller than this thresthold will be removed.
        area_thresthold             : float (optional)
            This parameter sets a subbasin area threshold in square kilometers (km²). 
            Subbasins with an area below this threshold will be merged with downstream 
            or upstream subbasins, depending on their location within the river network. 
            This can be useful for simplifying the river network and reducing computational
            requirements. The parameter is optional and has a default value of 0.009 km2. 
            If the parameter is set to 0, no subbasin merging will be performed.

        Returns
        -------

        Notes
        -----
        This function has no return values, The simplified hydrological routing
        network will be generated in the path_output_folder including following
        files:

        finalcat_info.shp                     : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            subbasin polygons which respect the lake inflow and outflow routing
            structures. This layer contains all the necessary information for
            hydrologic routing through the lake-river network.
        finalcat_info_riv.shp                 : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            river network polylines in the routing network.
        catchment_without_merging_lakes.shp             : shapefile
            The GIS layer containing subbasin polygons of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network subbasin
            polygons covered by the same lake are not merged into one lake subbasin
            yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        river_without_merging_lakes.shp                 : shapefile
            The GIS layer containing river polylines of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network, the
            river polylines covered by the same lake are not merged into one river
            segment yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        sl_connected_lake.shp                           : shapefile
            the GIS layer containing the lake polygons of lakes that are connected
            by the river_without_merging_lakes.shp
        sl_non_connected_lake.shp                       : shapefile
            the GIS layer containing the lake polygons of lakes that are not connected
            by the river_without_merging_lakes.shp
        poi                                      : shapefile
            It is the point shapefile that represent the point of interest
            after snap to river network.

        Examples
        -------


        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            simplify_routing_structure_by_drainage_area_method,
        )

        simplify_routing_structure_by_drainage_area_method(
            Routing_Product_Folder = routing_product_folder,
            Area_Min=minimum_subbasin_drainage_area,
            OutputFolder=path_output_folder,
            gis_platform=gis_platform,
            qgis_prefix_path=self.qgispp,
        )

        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = path_output_folder,
            qgis_prefix_path=self.qgispp,
            area_thresthold = area_thresthold,
            gis_platform=gis_platform,
        )

    def Add_Point_Of_Interest_Sites_In_Routing_Product(
        self,
        routing_product_folder= '#',
        path_to_points_of_interest_points = '#',
        path_output_folder="#",
        clean_exist_pois = True,
        gis_platform="purepy",
        area_thresthold = 10*30*30/1000/1000,
    ):
        """ The function allows the user to modify point of interest (POI) 
        related attributes in the routing product. The function is particularly
        useful when the user wants to add new POI sites, remove existing POI 
        sites,or modify the location of existing POI sites in the developed 
        routing products. If the user observes a small subbasin that cannot 
        be removed by the Decrease_River_Network_Resolution function due 
        to it representing an unexpected POI, they can use this function to 
        remove the POI from the routing product and then 
        apply the Decrease_River_Network_Resolution function.

        Parameters
        ----------
        path_output_folder                   : string
            is the folder path that stores generated outputs
        routing_product_folder         : string
            is the folder path where the input hydrologic routing network is stored
        path_to_points_of_interest_points                : string
            is the path to the point shapefile that contains the point of 
            interest (POI) sites. The shapefile must have an attribute table that includes 
            the following columns:
              - Obs_NM (string): This column should contain the name or ID of the POI site.
              - DA_Obs (float): This column should contain the drainage area of the POI site.
              - SRC_obs (string): This column should contain the source of the POI site.
        gis_platform                   : string
            is the parameter indicating which gis platform is used. Currenly, only "purepy" 
            is allowed for this parameter.
        clean_exist_pois                : boolean
            Indicate if user want to remove all existing POI in the input routing product.
        area_thresthold             : float (optional)
            This parameter sets a subbasin area threshold in square kilometers (km²). 
            Subbasins with an area below this threshold will be merged with downstream 
            or upstream subbasins, depending on their location within the river network. 
            This can be useful for simplifying the river network and reducing computational
            requirements. The parameter is optional and has a default value of 0.009 km2. 
            If the parameter is set to 0, no subbasin merging will be performed.

        Returns
        -------

        Notes
        -----
        This function has no return values, The modified hydrological routing
        network will be generated in the path_output_folder including following
        files:

        finalcat_info.shp                     : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            subbasin polygons which respect the lake inflow and outflow routing
            structures. This layer contains all the necessary information for
            hydrologic routing through the lake-river network.
        finalcat_info_riv.shp                 : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            river network polylines in the routing network.
        catchment_without_merging_lakes.shp             : shapefile
            The GIS layer containing subbasin polygons of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network subbasin
            polygons covered by the same lake are not merged into one lake subbasin
            yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        river_without_merging_lakes.shp                 : shapefile
            The GIS layer containing river polylines of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network, the
            river polylines covered by the same lake are not merged into one river
            segment yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        sl_connected_lake.shp                           : shapefile
            the GIS layer containing the lake polygons of lakes that are connected
            by the river_without_merging_lakes.shp
        sl_non_connected_lake.shp                       : shapefile
            the GIS layer containing the lake polygons of lakes that are not connected
            by the river_without_merging_lakes.shp
        poi                                      : shapefile
            It is the point shapefile that represent the point of interest
            after snap to river network.

        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            add_point_of_interest_sites_in_routing_product_method,
        )

        add_point_of_interest_sites_in_routing_product_method(
            routing_product_folder= routing_product_folder,
            path_to_points_of_interest_points = path_to_points_of_interest_points,
            path_output_folder = path_output_folder,
            clean_exist_pois = clean_exist_pois,
            gis_platform=gis_platform,
        )

        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = path_output_folder,
            qgis_prefix_path=self.qgispp,
            area_thresthold = area_thresthold,
            gis_platform=gis_platform,
        )

    def Select_Subregion_Of_Routing_Structure(
        self,
        path_output_folder="#",
        routing_product_folder = '#',
        gis_platform="qgis",
        most_down_stream_subbasin_ids=[-1],
        most_up_stream_subbasin_ids=[-1],
    ):
        """Select subregion of hydrologic routing network based on provided subbasin IDs

        Parameters
        ----------
        path_output_folder                   : string
            is the folder path that stores generated outputs
        routing_product_folder         : string
            is the folder path where the input hydrologic routing network is stored
        gis_platform                   : string
            is the parameter indicating which gis platform is used. It can be
            either "qgis" or "arcgis".
        most_down_stream_subbasin_ids                    : list
            A list of subbasin ID, the subbasin IDs in this list should
            be the most downstream subbasin ID of each interested watershed.
        most_up_stream_subbasin_ids                      : list
            A list of subbasin ID, the subbasins that drainage to the SubID in this list will be excluded.
            Value [-1] is required to indicate no upstream subbasin needs to be removed.

        Returns
        -------

        Notes
        -----
        This function has no return values, The extracted hydrological routing
        network will be generated in the path_output_folder including following
        files:

        finalcat_info.shp                     : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            subbasin polygons which respect the lake inflow and outflow routing
            structures. This layer contains all the necessary information for
            hydrologic routing through the lake-river network.
        finalcat_info_riv.shp                 : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            river network polylines in the routing network.
        catchment_without_merging_lakes.shp             : shapefile
            The GIS layer containing subbasin polygons of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network subbasin
            polygons covered by the same lake are not merged into one lake subbasin
            yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        river_without_merging_lakes.shp                 : shapefile
            The GIS layer containing river polylines of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network, the
            river polylines covered by the same lake are not merged into one river
            segment yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        sl_connected_lake.shp                           : shapefile
            the GIS layer containing the lake polygons of lakes that are connected
            by the river_without_merging_lakes.shp
        sl_non_connected_lake.shp                       : shapefile
            the GIS layer containing the lake polygons of lakes that are not connected
            by the river_without_merging_lakes.shp
        poi                                      : shapefile
            It is the point shapefile that represent the point of interest
            after snap to river network.

        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            select_part_of_routing_product_method,
        )
        Path_Points="#"
        Gauge_NMS=["#"]
        Path_Catchment_Polygon="#"
        Path_River_Polyline="#"
        Path_Con_Lake_ply="#"
        Path_NonCon_Lake_ply="#"

        select_part_of_routing_product_method(
            Path_Points=Path_Points,
            Gauge_NMS=Gauge_NMS,
            OutputFolder=path_output_folder,
            mostdownid=most_down_stream_subbasin_ids,
            mostupid=most_up_stream_subbasin_ids,
            Path_Catchment_Polygon=Path_Catchment_Polygon,
            Path_River_Polyline=Path_River_Polyline,
            Path_Con_Lake_ply=Path_Con_Lake_ply,
            Path_NonCon_Lake_ply=Path_NonCon_Lake_ply,
            qgis_prefix_path=self.qgispp,
            Routing_Product_Folder = routing_product_folder,
            gis_platform=gis_platform,
        )

        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = path_output_folder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

    def Generate_HRUs(
        self,
        path_output_folder,
        gis_platform,
        path_subbasin_polygon,
        path_landuse_info,
        path_soil_info,
        path_veg_info,
        projected_epsg_code='EPSG:3573',
        path_connect_lake_polygon="#",
        path_non_connect_lake_polygon="#",
        path_landuse_polygon="#",
        path_soil_polygon="#",
        path_vegetation_polygon="#",
        path_other_polygon_1="#",
        path_other_polygon_2="#",
        importance_order = [],
        min_hru_subbasin_area_ratio = 0.0,
        path_to_dem="#",
        pixel_size = 30,
        area_ratio_thresholds = [0,0,0]
    ):
        """This function is to generate HRU map based on subbasin polygon, lake
        polygon (optional), Land use polygon (optional), soil type polygon(optional),
        vegetation polygon (optional), and two other user defined polygons
        (optional).

        Parameters
        ----------
        path_output_folder                   : string
            is the folder path that stores generated outputs
        gis_platform                   : string
            is the parameter indicating which gis platform is used. It can be
            either "qgis" or "arcgis".
        path_subbasin_polygon                 : string
            is the path of the subbasin polygon, which is generated by
            BasinMaker.
        path_landuse_info                      : string
            Path to a csv file that contains landuse information, including
            the following attributes:

            | Landuse_ID (integer) -- the landuse ID in the landuse polygon,-1 for lake
            | LAND_USE_C (string) -- the landuse class name for each landuse type
        path_soil_info                        : string
            is the path to a csv file that contains soil information, including
            following attributes:

            | Soil_ID (integer) -- the soil ID  in the soil polygon,-1 for lake
            | SOIL_PROF (string) -- the soil profile name for each soil profile type
        path_veg_info                         : string
            is the path to a csv file that contains vegetation information, including
            following attributes:

            | Veg_ID (integer) -- the vegetation ID  in the vegetation polygon,-1 for lake
            | VEG_C (string) -- the vegetation class name for each vegetation Type
        projected_epsg_code                     : string (optional)
            is a EPSG code to indicate a projected coordinate system.
            If the routing network generated by basinmaker is under a
            projected system, please assgin the EPSG code of that projected
            system to projected_epsg_code. If the routing network generated by
            basinmaker is under a geographic system, please use EPSG code from
            any projected system. The coordinate system of routing network from
            basinmaker is determined by the coordinate system  of input DEM.
        path_connect_lake_polygon            : string (Optional)
            is the path to the connected lake's polygon
        path_non_connect_lake_polygon        : string (Optional)
            is the path to the non connected lake's polygon
        path_landuse_polygon                 : string (Optional)
            is the path to the landuse polygon. when path_landuse_polygon is not
            provided. The Landuse ID in path_landuse_info should be
            1: land, -1: lake
        path_soil_polygon                    : string (Optional)
            is the path to the soil polygon. when soil polygon is not
            provided. The Soil ID in path_soil_info should be the same
            as Landuse ID.
        path_vegetation_polygon                     : string (Optional)
            is the path to the vegetation polygon. when vegetation polygon is not
            provided. The Veg ID in path_veg_info should be the same
            as Landuse ID.
        path_other_polygon_1                 : string (Optional)
            is the path to the other polygon that will be used to define HRU,
            such as elevation band, or aspect.
        DEM                              : string (optional)
            is the path to a raster elevation dataset, that is used to
            calcuate average apspect, elevation and slope within each HRU.
            if no data is provided, subbasin averaged value will be used for
            each HRU.
        area_ratio_thresholds                  : list (optional)
            It is a list of HRU area thresthold for landuse, soil and other_1 layer respectively. In BasinMaker HRU calculation, each layer will firstly be overlaid to the subbasin map. Attributes falling into each subbasin with their area ratios (i.e., the intersected area to the subbasin area) smaller than the defined threshold values will then be dissolved into the largest part. For example, if forest area ratio in a subbasin, say subbasin #10, is 0.05, while we set the area threshold for land cover is 0.1. The forest polygons will then be dissolved to the largest land cover type in subbasin #10. 
        Returns
        -------

        Notes
        -----
        This function has no return values, but a HRU map saved in the
        path_output_folder

        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            generate_hrus_method,
        )
        Sub_Lake_ID="HyLakeId"
        Sub_ID="SubId"
        Lake_Id="Hylak_id"
        Landuse_ID="Landuse_ID"
        Soil_ID="Soil_ID"
        Other_Ply_ID_1="O_ID_1"
        Veg_ID="Veg_ID"
        Other_Ply_ID_2="O_ID_2"

        os.environ["GISDBASE"] = path_output_folder

        generate_hrus_method(
            Path_Subbasin_Ply=path_subbasin_polygon,
            Landuse_info=path_landuse_info,
            Soil_info=path_soil_info,
            Veg_info=path_veg_info,
            Sub_Lake_ID=Sub_Lake_ID,
            Sub_ID=Sub_ID,
            Path_Connect_Lake_ply=path_connect_lake_polygon,
            Path_Non_Connect_Lake_ply=path_non_connect_lake_polygon,
            Lake_Id=Lake_Id,
            Path_Landuse_Ply=path_landuse_polygon,
            Landuse_ID=Landuse_ID,
            Path_Soil_Ply=path_soil_polygon,
            Soil_ID=Soil_ID,
            Path_Veg_Ply=path_vegetation_polygon,
            Veg_ID=Veg_ID,
            Path_Other_Ply_1=path_other_polygon_1,
            Other_Ply_ID_1=Other_Ply_ID_1,
            Path_Other_Ply_2=path_other_polygon_2,
            Other_Ply_ID_2=Other_Ply_ID_2,
            DEM=path_to_dem,
            Inmportance_order = importance_order,
            min_hru_area_pct_sub = min_hru_subbasin_area_ratio,
            Project_crs=projected_epsg_code,
            OutputFolder=path_output_folder,
            qgis_prefix_path=self.qgispp,
            gis_platform = gis_platform,
            pixel_size = pixel_size,
            area_ratio_thresholds = area_ratio_thresholds,
        )


class delineate:

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
        fac_thresold = 5000,
    ):

        # define drived values
        # create folders
        self.path_output_folder = path_output_folder
        self.path_working_folder = path_working_folder
        self.fac_thresold = fac_thresold
        #        os.makedirs(self.path_output_folder, exist_ok=True)
        os.makedirs(self.path_working_folder, exist_ok=True)

        file = open(os.path.join(path_working_folder, "log.txt"), 'w')
        file.close()

        # obtain qgis prefix path
        if os.getenv("QGIS_PREFIX_PATH"):
            self.qgispp = os.environ["QGIS_PREFIX_PATH"]
        else:
            self.qgispp = '#'
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
    def Define_Project_Spatial_Extent(
        self,
        mode,
        path_to_dem_input,
        watershed_outlet_coordinates = [-1, -1],
        path_to_spatial_extent_polygon="#",
        buffer_distance=0.0,
        path_to_hydrobasin_polygon="#",
        hydrobasin_id_of_watershed_outlet=-1,
        path_to_snap_raster = "#",
        gis_platform="qgis",
    ):

        """This function is to define project spatial extent (PSE). Domain ouside
        of the PSE is not processed by BasinMaker functions.

        Parameters
        ----------
        mode                              : string (required)
            is a string indicating which to define PSE

            | 'using_dem'            : the extent of input dem is used
            | 'using_hybasin'        : the extent is defined by subbasins that are
                                       drainage to the provided watershed outlet
                                       subbasins ID in HydroBASINS product
            | 'using_outlet_pt'      : the extent is defined by the watershed
                                       generated from input dem and
                                       the watershed outlet coordinates
            | 'using_provided_ply'   : the extent of provided polygon is used
        path_to_dem_input            : string (required)
            is the path to input dem
        watershed_outlet_coordinates : list (optional)
            is a list that indicate the outlet coordinates of the
            region of interest in [lat, lon]. It is needed when mode = 'using_outlet_pt'.
        path_to_spatial_extent_polygon                  : string (optional)
            is the path of a polygon shapefile, the extent of which will be used
            as PSE.
        buffer_distance                  : float (optional)
            is a float number to enlarge the PSE. It is needed when mode = 'using_hybasin'
            or mode = 'using_provided_ply'. It is the distance around the PSE from
            HydroBASINS or provided polygons that will be buffered. The unit is the
            same with the spatial unit of input DEM.
        path_to_hydrobasin_polygon                      : string (optional)
            is a path to the HydroBASINs product. It is needed when mode = 'using_hybasin'
        hydrobasin_id_of_watershed_outlet               : int (optional)
            is a HydroBASINS subbasin ID of the watershed outlet. It is needed
            when mode = 'using_hybasin'

        Returns
        -------

        Notes
        -------
        Outputs are following files in GRASS GIS database loacated in
        os.path.join(path_working_folder,'grassdb')

        MASK.*        : raster/shp
            it is a mask raster stored in grass database, which indicate the PSE.
        dem.*         : raster/shp
            it is a dem raster stored in grass database, which is has the same extent
            with MASK.

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
            path_dem_in=path_to_dem_input,
            outlet_pt=watershed_outlet_coordinates,
            path_extent_ply=path_to_spatial_extent_polygon,
            buffer_distance=buffer_distance,
            hybasin_ply=path_to_hydrobasin_polygon,
            down_hybasin_id=hydrobasin_id_of_watershed_outlet,
            up_hybasin_id=up_hybasin_id,
            path_to_snap_raster = path_to_snap_raster,
            mask=self.geofilenames["mask"],
            dem=self.geofilenames["dem"],
            gis_platform=gis_platform,
        )

    def Delineation_Initial_Subbasins_Without_Lakes(
        self,
        fac_thresold,
        mode = 'using_dem',
        path_flow_dirction="#",
        path_flow_accumulation='#',
        max_memroy=1024 * 4,
        gis_platform="qgis",
        subreg_fdr_path="#",
        subreg_acc_path="#",
        subreg_str_r_path="#",
        subreg_str_v_path="#",
    ):
        """Function that used to generate a initial routing structure with
        user provied flow accumulation thresthold without considering lake.

        Parameters
        ----------
        fac_thresold       : float
            is the flow accumulation thresthold, used to determine
            subbsains and river network. Increasing of this paramter will
            increase the size of generated subbasins, reduce the number
            subbasins and reduce the number of generated stream reaches
        mode              : string (required)
            is a string indicate which dataset will be used to delineate
            watershed.

            | 'using_dem' : dem is used for initial subbasin delineation
            | 'using_fdr' : flow direction data is used for subbasin delineation
        path_flow_dirction          : string (optional)
            is a path indicating the path of flow direction input dataset
            only needed when mode = 'using_fdr'
        max_memroy        : integer
            is the maximum memeory that allow to be used in MB.

        Returns
        -------

        Notes
        -------
        Outputs are following files in GRASS GIS database loacated in
        os.path.join(path_working_folder,'grassdb')

        fdr_grass              : raster
            is a raster represent flow direction dataset, which is
            using 1 - 8 to represent different directions
        fdr_arcgis             : raster
            is a raster represent flow direction dataset, which is
            using 1,2,4,...64,128 to represent different directions
        str_v                  : vector
            is a river network in vector format
        str_r                  : raster
            is a river network in raster format
        cat_no_lake            : raster
            is the raster represent the delineated subbasins without
            considering lakes
        acc                    : raster
            is the raster represent the flow accumulation


        Examples
        -------
        """
        from basinmaker.delineationnolake.watdelineationwithoutlake import (
            watershed_delineation_without_lake,
        )
        if mode == 'using_dem':
            mode = 'usingdem'
        elif mode == 'using_sub_reg':
            mode = 'usingsubreg'
        elif mode == 'using_fac':
            mode = 'usingfac'
        else:
            mode = 'usingfdr'

        watershed_delineation_without_lake(
            mode=mode,
            input_geo_names=self.geofilenames,
            acc_thresold=fac_thresold,
            fdr_path=path_flow_dirction,
            fac_path = path_flow_accumulation,
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

    def Add_New_Subbasin_Outlet_Points(
        self,
        path_lake_polygon="#",
        lake_attributes=[],
        connected_lake_area_thresthold=0,
        non_connected_lake_area_thresthold=0,
        only_included_lake_at_river_interction = False,
        path_point_of_interest="#",
        point_of_interest_attributes=[],
        path_sub_reg_outlets_v="#",
        path_sub_reg_lake_r="#",
        path_sub_reg_lake_bd_r="#",
        search_radius=100,
        mode="#",
        max_memroy=1024 * 4,
        gis_platform="qgis",
    ):
        """ Update the subbasin delineation result by adding lake inflow and
        outflow points and observation gauges as a new subbasin outlets. The output
        is not the final delineation result. because:

        | 1) Hydrologcial related attributes for each subbasin are not calcuated yet.
        | 2) Some lakes may cover several subbasins. The output needs to be finalized
             by combing those subbasins with the same lake as one subbasin only.


        Parameters
        ----------
        path_lake_polygon                 : string (optional)
            is a path of the lake polygon shapefile
        lake_attributes                   : list (optional)
            the columns names in the input lake polygon that indicate following
            items (mandatory). It is needed only when path_lake_polygon_in is
            provided. Columns (2-4 in following list) in lake poylon can be fill
            with any number, when these infomation is not avaiable.

            | 1) column name for the unique Id of each lake, datatype is integer
            | 2) column name for type of the lake, datatype is integer
            | 3) column name for the volume of lakes in km3, datatype is float
            | 4) column name for the average depth of lakes in m, datatype is float
            | 5) column name for the area of lakes in km2, datatype is float
        connected_lake_area_thresthold                 : float (optional)
            is a lake area thresthold for connected lakes in km2.
            Connected lake with lake area below this value will not be considerd
        non_connected_lake_area_thresthold             : float (optional)
            is a lake area thresthold for non-connected lakes in km2
            Non connected lake with lake area below this value will not be considered
        path_point_of_interest
            is a path of the point shapefile that indicate points of interest,which
            can include different observation gauges.
        point_of_interest_attributes                    : list (optional)
            the columns names in the point of interest shapefile that indicate
            following items (mandatory). It is needed only when path_point_of_interest
            is provided. Columns (2-4 in following list) in point shapefile can be
            fill with any value, when these infomation is not avaiable.

            | 1) column name for the unique Id of each observation point, datatype is integer
            | 2) column name for the unique name of each observation point, datatype is string
            | 3) column name for the drainage area of each observation point in km3, datatype is float
            | 4) column name for the source of the observation point: 'CA' for observation in canada;
                 'US' for observation in US, or any user-provided names, datatype is string
        max_memroy        : integer (optional)
            is the maximum memeory that allow to be used in MB.

        Returns
        -------

        Notes
        -------
        Output raster and vector files that will be used by next step are list as
        following. All files are saved in GRASS GIS database loacated in
        os.path.join(path_working_folder,'grassdb')

        selected_lakes                    : raster
            it is a raster represent all lakes that are selected by two lake
            area threstholds
        sl_nonconnect_lake       : raster
            it is a raster represent all non connected lakes that are selected
            by lake area threstholds
        sl_connected_lake           : raster
            it is a raster represent all connected lakes that are selected
            by lake area threstholds
        river_without_merging_lakes                         : raster/vector
            it is the updated river segment for each subbasin
        catchment_without_merging_lakes                     : raster/vector
            it is a raster represent updated subbasins after adding lake inflow
            and outflow points as new subbasin outlet.
        snapped_obs_points                                  : raster/vector
            it is a name of the point gis file represent successfully sanpped
            point of interest points


        Examples
        -------

        """
        from basinmaker.addlakeandobs.addlakeandobsintowatershed import (
            add_lakes_and_obs_into_existing_watershed_delineation,
        )

        add_lakes_and_obs_into_existing_watershed_delineation(
            input_geo_names=self.geofilenames,
            path_lakefile_in=path_lake_polygon,
            lake_attributes=lake_attributes,
            path_obsfile_in=path_point_of_interest,
            obs_attributes=point_of_interest_attributes,
            path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            threshold_con_lake=connected_lake_area_thresthold,
            only_included_lake_at_river_interction = only_included_lake_at_river_interction,
            threshold_non_con_lake=non_connected_lake_area_thresthold,
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
            fac_thresold = self.fac_thresold,
            max_memroy=max_memroy,
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
        )

    def Generate_Hydrologic_Routing_Attributes(
        self,
        path_output_folder ="#",
        projected_epsg_code ="EPSG:3573",
        path_bkfwidthdepth_polyline ="#",
        bkfwd_attributes=[],
        path_k_c_zone_polygon = '#',
        k =-1,
        c =-1,
        path_landuse ="#",
        path_landuse_and_manning_n_table="#",
        lake_attributes =[],
        point_of_interest_attributes =[],
        outlet_obs_id=-1,
        path_sub_reg_outlets_v="#",
        gis_platform="qgis",
    ):
        """Calculate hydrological paramters for each subbasin.

        Parameters
        ----------
        path_output_folder                  : string
            The path to a folder to save outputs
        projected_epsg_code                     : string (optional)
            is a EPSG code to indicate a projected coordinate system.
            If the routing network generated by basinmaker is under a
            projected system, please assgin the EPSG code of that projected
            system to projected_epsg_code. If the routing network generated by
            basinmaker is under a geographic system, please use EPSG code from
            any projected system. The coordinate system of routing network from
            basinmaker is determined by the coordinate system  of input DEM.
        path_bkfwidthdepth_polyline             : string (optional)
            is a path of the polyline shapefile that contains bankfull width (w) and
            depth (d) data.Following the methodology in Andreadis et al. (2013),
            the w, and d of each subbasin can be calculated by these two equation
            w= 7.2 Q**0.5 and d=0.27Q**0.39 , respectively. And the the bankfull
            discharge Q of each subbasin can be estimated using the relationship
            between Q and the drainage area (DA in [km2]), which is the Q=k×DA**c.
            If the bankfull width and depth polyline data is provided, basinmaker
            will use it to estimate the k and c for the entire watershed.
        bkfwd_attributes               : list (optional)
            the columns names that indicate following items (mandatory).It is only
            needed with path_bkfwidthdepth is provided

            | 1) column name for the Bankfull width in m;
            | 2) column name for the Bankfull depth in m;
            | 3) column name for the annual mean discharge in m3/s;
            | 4) column name for the drainage area in km2;
        k                              : float (optional)
            the coefficient in Q=k×DA**c.if k and c is provided, the
            path_bkfwidthdepth_polyline will not be used
        c                              : float (optional)
            the coefficient in Q=k×DA**c.if k and c is provided, the
            path_bkfwidthdepth_polyline will not be used
        path_landuse                   : string (optional)
            is a path of the landuse raster.It is used to estimate the floodplain
            Manning's coefficient. Require the same projection with the DEM data
            and formatted in ".tif".
        path_landuse_and_manning_n_table              : string (optional)
            is a path of the table in '.csv' format.The table describe the
            floodplain Manning's coefficient correspond to a given landuse type.
            The table should have two columns:

            | RasterV: is the landuse value in the landuse raster for each land use type
            | MannV: is the roughness coefficient value for each landuse type.
        lake_attributes                   : list (optional)
            the columns names in the input lake polygon that indicate following
            items (mandatory). It is needed only when path_lake_polygon_in is
            provided. Columns (2-4 in following list) in lake poylon can be fill
            with any number, when these infomation is not avaiable.

            | 1) column name for the unique Id of each lake, datatype is integer
            | 2) column name for type of the lake, datatype is integer
            | 3) column name for the volume of lakes in km3, datatype is float
            | 4) column name for the average depth of lakes in m, datatype is float
            | 5) column name for the area of lakes in km2, datatype is float
        point_of_interest_attributes                    : list (optional)
            the columns names in the point of interest shapefile that indicate
            following items (mandatory). It is needed only when path_point_of_interest
            is provided. Columns (2-4 in following list) in point shapefile can be
            fill with any value, when these infomation is not avaiable.

            | 1) column name for the unique Id of each observation point, datatype is integer
            | 2) column name for the unique name of each observation point, datatype is string
            | 3) column name for the drainage area of each observation point in km3, datatype is float
            | 4) column name for the source of the observation point: 'CA' for observation in canada;
                 'US' for observation in US, or any user-provided names, datatype is string

        Returns
        -------

        Notes
        -------
        Five vector files will be generated in the output folder. these files
        can be further finalized as hydrologcial routing network by
        function "combine_catchments_covered_by_the_same_lake"
        or be used as input for BasinMaker post processint tools.
        The attributes included in each GIS file can be found in
        http://hydrology.uwaterloo.ca/basinmaker/index.html

        catchment_without_merging_lakes.shp             : shapefile
            The GIS layer containing subbasin polygons of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network subbasin
            polygons covered by the same lake are not merged into one lake subbasin
            yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        river_without_merging_lakes.shp                 : shapefile
            The GIS layer containing river polylines of an incomplete hydrologic
            routing network. In this incomplete hydrologic routing network, the
            river polylines covered by the same lake are not merged into one river
            segment yet. This incomplete hydrologic routing network is only intended as
            input to customize the routing network with our BasinMaker GIS toolbox
            (for example by defining new lake area thresholds and/or a new catchment
            minimum drainage area threshold)
        sl_connected_lake.shp                           : shapefile
            the GIS layer containing the lake polygons of lakes that are connected
            by the river_without_merging_lakes.shp
        sl_non_connected_lake.shp                       : shapefile
            the GIS layer containing the lake polygons of lakes that are not connected
            by the river_without_merging_lakes.shp
        poi                                       : shapefile
            It is the point shapefile that represent the point of interest
            after snap to river network.


        Examples
        -------

        """

        from basinmaker.addattributes.addattributestocatchments import add_attributes_to_catchments

        add_attributes_to_catchments(
            input_geo_names=self.geofilenames,
            path_bkfwidthdepth=path_bkfwidthdepth_polyline,
            bkfwd_attributes=bkfwd_attributes,
            path_landuse=path_landuse,
            path_landuse_info=path_landuse_and_manning_n_table,
            projection=projected_epsg_code,
            path_k_c_zone_polygon = path_k_c_zone_polygon,
            k_in=k,
            c_in=c,
            out_cat_name=self.geofilenames["catchment_without_merging_lakes"],
            out_riv_name=self.geofilenames["river_without_merging_lakes"],
            grassdb=self.grassdb,
            grass_location=self.grass_location_geo,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
            obs_attributes=point_of_interest_attributes,
            lake_attributes=lake_attributes,
            outlet_obs_id=outlet_obs_id,
            path_sub_reg_outlets_v=path_sub_reg_outlets_v,
            output_folder=path_output_folder,
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

    def Combine_Subbasins_Covered_by_The_Same_Lake(
        self,
        routing_product_folder='#',
        gis_platform="qgis",
        area_thresthold = 10*30*30/1000/1000,
    ):
        """Finalize a incomplete hydrologic routing network by merging subbasin
        polygons that are covered by the same lake.

        Parameters
        ----------
        routing_product_folder         : string
            is the folder where the input routing product is stored
        gis_platform                   : string
            It is the parameter indicate which gis platform is used. It can be
            either "qgis" or "arcgis".
        Returns
        -------

        Notes
        -----
        This function has no return values, two vector files will be generated in
        the routing_product_folder

        finalcat_info.shp                     : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            subbasin polygons which respect the lake inflow and outflow routing
            structures. This layer contains all the necessary information for
            hydrologic routing through the lake-river network.
        finalcat_info_riv.shp                 : shapefile
            The finalized hydrologic routing network. the GIS layer containing
            river network polylines in the routing network.

        Examples
        -------

        """
        from basinmaker.postprocessing.postprocessingfunctions import (
            combine_catchments_covered_by_the_same_lake_method,
        )

        combine_catchments_covered_by_the_same_lake_method(
            Routing_Product_Folder = routing_product_folder,
            qgis_prefix_path=self.qgispp,
            gis_platform=gis_platform,
            area_thresthold = area_thresthold,
        )
