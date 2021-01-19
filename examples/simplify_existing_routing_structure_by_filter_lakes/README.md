# Example - simplify the existing lake-river routing structure by filtering lakes
 
# Overview
In this example, the existing lake-river routing structure will be simplified by removing lakes. Lake area threshold will be needed. Lakes smaller than lake area threshold will be removed from the existing lake-river routing structure 
```
python example_simplify_routing_structure_by_filter_lakes.py
```

# Explanation of script 

## Define the working and output folder 
The working folder is the folder where basinmaker will save all geo spatial files generated during the watershed delineation; while the output folder is the folder the basinmaker will save final outputs.

```
#############################################
# define working folder, output folder amd data folder  
#############################################
num  = str(np.random.randint(1, 10000 + 1))
path_output_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_merit" +num,"output")
path_working_folder = os.path.join(tempfile.gettempdir(), "basinmaker_exp_merit" +num,"work")
datafolder = os.path.join("../../tests/testdata", "Required_data_to_start_from_dem")
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

## simplify existing routing structure by filtering lakes 
The basinmaker can be initialized with path to the working folder 

### parameters 

* \<OutputFolder\> 
      is the folder that stores generated outputs

* \<Path_final_riv_ply\> 
      Path to the catchment polygon which is the routing product before merging lakes catchment and need to be processed before used. It is the input for simplify the routing product based on lake area or drianage area.

* \<Path_final_riv\> 
      Path to the river polyline which is the routing product before merging lakes catchments and need to be processed before used. It is the input for simplify the routing product based on lake area or drianage area.

* \<Path_Con_Lake_ply\> 
      Path to a connected lake polygon. Connected lakes are lakes that are connected by Path_final_riv.

* \<Path_NonCon_Lake_ply\> 
      Path to a non connected lake polygon. Connected lakes are lakes that are not connected by Path_final_riv.

* \<Thres_Area_Conn_Lakes\> 
      It is the lake area threshold for connated lakes, in km2      

* \<Thres_Area_Non_Conn_Lakes\> 
      It is the lake area threshold for non connated lakes, in km2

* \<Selection_Method\> 
      It is a string indicate lake selection methods "ByArea" means lake in the routing product will be selected based on two lake area thresthold Thres_Area_Conn_Lakes and Thres_Area_Non_Conn_Lakes; "ByLakelist" means lake in the routing product will be selected based on user provided hydrolake id, in Selected_Lake_List_in

* \<Selected_Lake_List_in\> 
      A list of lake ids that will be keeped in the routing product. Lakes not in the list will be removed from routing product.

* \<gis_platform\> 
    It is the parameter indicate which gis platform is used. For now only "qgis" is supported   


```
#############################################
# simplify existing routing structure by filtering lakes   
#############################################
basinmaker.simplify_routing_structure_by_filter_lakes_method(
    OutputFolder=path_output_folder,
    Path_final_riv_ply=os.path.join(
        datafolder, "catchment_without_merging_lakes.shp"
    ),
    Path_final_riv=os.path.join(datafolder, "river_without_merging_lakes.shp"),
    Path_Con_Lake_ply=os.path.join(datafolder, "sl_connected_lake.shp"),
    Path_NonCon_Lake_ply=os.path.join(
        datafolder, "sl_non_connected_lake.shp"
    ),
    Thres_Area_Conn_Lakes=5,
    Thres_Area_Non_Conn_Lakes=1,
    Selection_Method="ByArea",
    Selected_Lake_List_in=[],
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

* \<Path_final_rivply\> 
      is path to the catchment polygon represent watershed delineation after add lake and gauges control points and update hydrological attributes  

* \<Path_final_riv\> 
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
