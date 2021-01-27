# Example - extract part of the existing lake-river routing structure 
 
# Overview
In this example, the part existing lake-river routing structure will be extracted based on either:1) a point shpfile, all catchments drainage to the catchment contain provided point will be extracted; 2) most down stream subbasin id and most upstream subbasin id. all catchment located between these two subbasin ids will be extracted; 3)a gauge name, all catchment drainage to the catchment contain the gauge name in the attribute table will be extracted.  
 
```
python example_select_part_of_routing_product.py
```

# Explanation of script 

## Define the working and output folder 
The working folder is the folder where basinmaker will save all geo spatial files generated during the watershed delineation; while the output folder is the folder the basinmaker will save final outputs.

```
#############################################
# define working folder, output folder amd data folder  
#############################################
num  = str(np.random.randint(1, 10000 + 1))
path_output_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_" +num,"output")
path_working_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_" +num,"work")
datafolder = os.path.join("../../tests/testdata", "existing_lake_river_routing_structure")
HRU_Folder = os.path.join("../../tests/testdata/", "HRU")
```


## Initialize basinmaker 
The basinmaker can be initialized with path to the working folder 

```
#############################################
# initialize basinmaker with working folder    
#############################################
basinmaker = basinmaker(
    path_working_folder=path_working_folder
)
```


## extract part of the existing routing structure  

### parameters 

* \<OutputFolder\> 
      is the folder that stores generated outputs

* \<Path_Catchment_Polygon\> 
      Path to the catchment polygon which is the routing product before merging lakes catchment and need to be processed before used. It is the input for simplify the routing product based on lake area or drianage area.

* \<Path_River_Polyline\> 
      Path to the river polyline which is the routing product before merging lakes catchments and need to be processed before used. It is the input for simplify the routing product based on lake area or drianage area.

* \<Path_Con_Lake_ply\> 
      Path to a connected lake polygon. Connected lakes are lakes that are connected by Path_final_riv.

* \<Path_NonCon_Lake_ply\> 
      Path to a non connected lake polygon. Connected lakes are lakes that are not connected by Path_final_riv.

* \<Gauge_NMS\> 
      a list contain gauge name      

* \<mostdownid\> 
      the most downstream subbasin id       

* \<mostupid\> 
      the most upstream subbasin id       

* \<Path_Points\> 
      the path to the point shp file       

* \<gis_platform\> 
      It is the parameter indicate which gis platform is used. For now only "qgis" is supported   

```
#############################################
# obtain part of the routing sturcture by gauge name or
# subbasin ids   
#############################################
basinmaker.select_part_of_routing_product_methods(
    OutputFolder=path_output_folder,
    Path_Catchment_Polygon=os.path.join(
        datafolder, "catchment_without_merging_lakes.shp"
    ),
    Path_River_Polyline=os.path.join(datafolder, "river_without_merging_lakes.shp"),
    Path_Con_Lake_ply=os.path.join(datafolder, "sl_connected_lake.shp"),
    Path_NonCon_Lake_ply=os.path.join(datafolder, "sl_non_connected_lake.shp"),
    Gauge_NMS= ['#'],  #['02KB001'],
    mostdownid = 1388,
    mostupid = 1198,
    Path_Points = '#',
    gis_platform="qgis",
)
```

### outputs  
* \<sl_connected_lake\> 
      is selected non connected lakes, located in output folder  

* \<sl_non_connected_lake\> 
      is selected connected lakes, located in output folder 

* \<river_without_merging_lakes\> 
      is the river network after add lake and gauge control points, located in output folder  

* \<catchment_without_merging_lakes\> 
      is catchment polygons after add lake and gauge control points, located in output folder 

* \<obs_gauges\> 
      is the gauge points after snapped to closest river network, located in output folder             


## combine catchment covered by the same lake
In this step, catchments covered by the same lake will be combined together. and the final lake-river routing network will be saved in the output folder 

### parameters 

* \<output_folder\> 
      is the path to the output folder  

* \<Path_final_info\> 
      is path to the catchment polygon represent watershed delineation after add lake and gauges control points and update hydrological attributes  

* \<Path_final_inforiv\> 
      is path to the catchment river polyline represent watershed delineation after add lake and gauges control points and update hydrological attributes  

* \<gis_platform\> 
      It is the parameter indicate which gis platform is used. For now only "qgis" is supported  
 
```
#############################################
# combine catchments covered by the same lakes 
#############################################
basinmaker.combine_catchments_covered_by_the_same_lake_method(
    OutputFolder=path_output_folder,
    Path_final_rivply=os.path.join(
        path_output_folder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(path_output_folder, "river_without_merging_lakes.shp"),
    gis_platform="qgis",
)
```
### outputs  
* \<finalcat_info\> 
      is catchment of final lake-river routing structure in output folder 
* \<finalcat_info_riv\> 
      is the river of the final lake-river routing structure in output folder  
