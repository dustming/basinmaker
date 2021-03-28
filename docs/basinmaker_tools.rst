================
BasinMaker tools
================


BasinMaker - postprocessing tools
=================================

Extract the region of interest
------------------------------

.. autofunction:: basinmaker.postproc.select_part_of_routing_product

Filter lakes
------------

.. autofunction:: basinmaker.postproc.simplify_routing_structure_by_filter_lakes

Increase catchment area
-----------------------

.. autofunction:: basinmaker.postproc.simplify_routing_structure_by_drainage_area

Generate HRUs
-------------

.. autofunction:: basinmaker.postproc.generate_hrus

Generate Raven input files
--------------------------

.. autofunction:: basinmaker.postproc.generate_raven_model_inputs



BasinMaker - delineate lake-river routing product tools
=======================================================

Define project spatial extent
-----------------------------

.. autofunction:: basinmaker.basinmaker.define_project_extent_method(mode,path_dem_in,outlet_pt=[-1,-1],path_extent_ply = '#',buffer_distance=0.0,hybasin_ply='#',down_hybasin_id=-1)


Delineate routing structure without lakes
-----------------------------------------

.. autofunction:: basinmaker.basinmaker.watershed_delineation_without_lake_method(acc_thresold,fdr_path,mode,max_memroy)


Add lake and obs control points
-------------------------------

.. autofunction:: basinmaker.basinmaker.watershed_delineation_add_lake_control_points(path_lakefile_in,lake_attributes,threshold_con_lake,threshold_non_con_lake,path_obsfile_in,obs_attributes,max_memroy)


Add hydrology related attributes
--------------------------------

.. autofunction:: basinmaker.basinmaker.add_attributes_to_catchments_method(path_bkfwidthdepth="#",bkfwd_attributes=[],path_landuse="#",path_landuse_info="#",projection="EPSG:3573",lake_attributes=[],obs_attributes=[],OutputFolder="#")


Combine catchment covered by the same lake
------------------------------------------

.. autofunction:: basinmaker.basinmaker.combine_catchments_covered_by_the_same_lake




