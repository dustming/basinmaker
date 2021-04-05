================
BasinMaker tools
================


BasinMaker - postprocessing tools
=================================

Extract the region of interest
------------------------------

.. autofunction:: basinmaker.basinmaker.postproc.select_part_of_routing_product

Filter lakes
------------

.. autofunction:: basinmaker.basinmaker.postproc.simplify_routing_structure_by_filter_lakes

Increase catchment area
-----------------------

.. autofunction:: basinmaker.basinmaker.postproc.simplify_routing_structure_by_drainage_area

Generate HRUs
-------------

.. autofunction:: basinmaker.basinmaker.postproc.generate_hrus

Generate Raven input files
--------------------------

.. autofunction:: basinmaker.basinmaker.postproc.generate_raven_model_inputs



BasinMaker - delineate lake-river routing product tools
=======================================================

Define project spatial extent
-----------------------------

.. autofunction:: basinmaker.basinmaker.dlidem.Define_Project_Spatial_Extent(mode,path_to_dem_input,watershed_outlet_coordinates=[-1,-1],path_to_spatial_extent_polygon = '#',buffer_distance=0.0,path_to_hydrobasin_polygon='#',hydrobasin_id_of_watershed_outlet=-1)


Delineate routing structure without lakes
-----------------------------------------

.. autofunction:: basinmaker.basinmaker.dlidem.Delineation_Initial_Subbasins_Without_Lakes(flow_accumulation_thresold,mode = 'using_dem',path_flow_dirction_in = '#',max_memroy = 4096)


Add lake and obs control points
-------------------------------

.. autofunction:: basinmaker.basinmaker.dlidem.Add_New_Subbasin_Outlet_Points(path_lake_polygon = '#',lake_attributes=[],connected_lake_area_thresthold = 0,non_connected_lake_area_thresthold = 0,path_point_of_interest = '#',point_of_interest_attributes = [],max_memroy = 4096)


Add hydrology related attributes
--------------------------------

.. autofunction:: basinmaker.basinmaker.dlidem.Generate_Hydrologic_Routing_Attributes(path_output_folder,prjected_epsg_code = "EPSG:3573",path_bkfwidthdepth_polyline="#",bkfwd_attributes=[],k = -1,c=-1,path_landuse="#",path_landuse_and_manning_n_table="#",lake_attributes=[],point_of_interest_attributes=[])


Combine catchment covered by the same lake
------------------------------------------

.. autofunction:: basinmaker.basinmaker.dlidem.Combine_Subbasins_Covered_by_The_Same_Lake(routing_product_folder,gis_platform="qgis")




