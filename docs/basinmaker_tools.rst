==================================
BasinMaker developer documentation
==================================


BasinMaker - postprocessing tools
=================================

Extract the region of interest
------------------------------

.. autofunction:: basinmaker.basinmaker.postprocess.Select_Subregion_Of_Routing_Structure( path_output_folder,routing_product_folder,gis_platform,most_down_stream_subbasin_ids = [],most_up_stream_subbasin_ids = [])

Filter lakes
------------

.. autofunction:: basinmaker.basinmaker.postprocess.Remove_Small_Lakes(path_output_folder,routing_product_folder,gis_platform,connected_lake_area_thresthold = -1,non_connected_lake_area_thresthold = -1,selected_lake_ids=[],area_thresthold = 0.009)

Increase catchment area
-----------------------

.. autofunction:: basinmaker.basinmaker.postprocess.Decrease_River_Network_Resolution(path_output_folder,routing_product_folder,gis_platform,minimum_subbasin_drainage_area,area_thresthold)

Modify point of interest in the routing product
-----------------------

.. autofunction:: basinmaker.basinmaker.postprocess.Add_Point_Of_Interest_Sites_In_Routing_Product(path_output_folder,routing_product_folder,gis_platform,clean_exist_pois,area_thresthold)

Generate HRUs
-------------

.. autofunction:: basinmaker.basinmaker.postprocess.Generate_HRUs(path_output_folder,gis_platform,path_subbasin_polygon,path_landuse_info,path_soil_info,path_veg_info,prjected_epsg_code='EPSG:3573',path_connect_lake_polygon="#",path_non_connect_lake_polygon="#",path_landuse_polygon='#',path_soil_polygon="#",path_vegetation_polygon="#",path_other_polygon_1="#",area_ratio_thresholds = [0,0,0],path_to_dem = "#")

Generate Raven input files
--------------------------

.. autofunction:: basinmaker.basinmaker.postprocess.Generate_Raven_Model_Inputs(path_output_folder,path_hru_polygon,aspect_from_gis,model_name="test",subbasingroup_nm_channel=["Allsubbasins"],subbasingroup_length_channel=[-1],subbasingroup_nm_lake=["AllLakesubbasins"],subbasingroup_area_lake=[-1])



BasinMaker - delineate lake-river routing product tools
=======================================================

Define project spatial extent
-----------------------------

.. autofunction:: basinmaker.basinmaker.delineate.Define_Project_Spatial_Extent(mode,path_to_dem_input,watershed_outlet_coordinates=[-1,-1],path_to_spatial_extent_polygon = '#',buffer_distance=0.0,path_to_hydrobasin_polygon='#',hydrobasin_id_of_watershed_outlet=-1)


Delineate routing structure without lakes
-----------------------------------------

.. autofunction:: basinmaker.basinmaker.delineate.Delineation_Initial_Subbasins_Without_Lakes(fac_thresold,mode = 'using_dem',path_flow_dirction_in = '#',max_memroy = 4096)


Add lake control points and points of interest
-------------------------------

.. autofunction:: basinmaker.basinmaker.delineate.Add_New_Subbasin_Outlet_Points(path_lake_polygon = '#',lake_attributes=[],connected_lake_area_thresthold = 0,non_connected_lake_area_thresthold = 0,path_point_of_interest = '#',point_of_interest_attributes = [],max_memroy = 4096)


Guidance on Points of Interest input layer preparation
-------------------------------

For those running BasinMaker with QGIS & GRASS, points of interest will be snapped by this function automatically to 
the closest delineated river channel and as such, no input layer preparation is strictly required. However, in our 
experience this auto-snapping approach is only modestly successful. Please inspect snapped results carefully. If users 
ensure that points of interest associated with lake levels are located within a lake polygon, then this function will 
include the point of interest at the lake (and lake subbasin) outlet.


For those running BasinMaker in ArcGIS pro, the user needs to carefully prepare the points of interest shapefile 
as input. Specifically, the location of each non-lake point of interest should be carefully snapped to the the appropriate 
eventual delineated river channel. To do this, we recommend users build a temporary channel raster and snap to that.
We also recommend users ensure that points of interest associated with lake levels are located within a lake polygon. 
These POI will be ignored by this function and must be added separately from this function using the function called 
"Add_Point_Of_Interest_Sites_In_Routing_Product".
  
A robust new function that can help users automatically snap points of interest to delineated river 
channels/lakes is coming in the next version of BasinMaker.


Add hydrology related attributes
--------------------------------

.. autofunction:: basinmaker.basinmaker.delineate.Generate_Hydrologic_Routing_Attributes(path_output_folder,prjected_epsg_code = "EPSG:3573",path_bkfwidthdepth_polyline="#",bkfwd_attributes=[],k = -1,c=-1,path_landuse="#",path_landuse_and_manning_n_table="#",lake_attributes=[],point_of_interest_attributes=[])


Combine catchment covered by the same lake
------------------------------------------

.. autofunction:: basinmaker.basinmaker.delineate.Combine_Subbasins_Covered_by_The_Same_Lake(routing_product_folder,gis_platform="qgis")




