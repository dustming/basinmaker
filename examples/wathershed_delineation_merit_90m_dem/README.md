# Example - delineate lake river routing structure with MERIT 90 m DEM
 
# Overview
In this example, the MERIT 90 m DEM will be used to delineate the lake river routing structure. The extent of the watershed boundary will be estimated by an provided outlet coordinates. The basinmaker is called in the in the example python script "example_watershed_delineation_with_lakes_merit_90m.py". After setup required basinmaker working environments (instruction can be found at [here](https://github.com/dustming/basinmaker/wiki/Installation-of-the-BasinMaker)). This script can be run by following commmand:
```
python example_watershed_delineation_with_lakes_merit_90m.py
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


## Define project spatial extent
The first step of watershed delineation is to define the project spatial extent. Once the project spatial extent is defined. The basinmaker will delineate lake-river routing structure within this extent. 

In this example, the project spatial extent will be defined using a wathershed outlet coordinates and a MERIT 90m DEM.

### parameters 

* \<define_project_extent_method\> 
      is the function that can be used to define project extent 
* \<mode\> 
      is a parameter indicate which method will be used to define project extent. In this case, "using_outlet_pt" is used. 
* \<path_dem_in\> 
      is a full path to the DEM dataset. 
* \<outlet_pt\> 
      is the coordinates of the watershed outlet point in [lon,lat]
* \<buffer_distance\> 
      is the parameter that used to buffer the extent defined by \<outlet_pt\>      
 
```
#############################################
# define extent of the processing domain  
#############################################
basinmaker.define_project_extent_method(
    mode="using_outlet_pt",
    path_dem_in=os.path.join(datafolder, "DEM_big_merit.tif"),
    outlet_pt=[-92.387, 49.09],
    buffer_distance=0.00,
)
```
### outputs  
* \<dem\> 
      is the dem within the project extent, located in working folder 
* \<mask\> 
      is the a mask to represent project extent, located in working folder  


## Delineate watershed without considering lakes 
After obtain dem within the project extent. basinmaker will first generate a watershed delineation without considering lakes. 

### parameters 

* \<acc_thresold\> 
      is the flow accumulation threshold, the density of the generated river network and the size of the catchment will be decreasing with increasing of this parameter   
* \<mode\> 
      is a parameter indicate which method will be used to delineate watershed without considering lakes. In this example dem is used, with "usingdem".
* \<max_memroy\> 
      It is the maximum memory allow to be used by basinmaker  
* \<gis_platform\> 
      It is the parameter indicate which gis platform is used. For now only "qgis" is supported 
 
```
#############################################
# generate a watershed delineation without considering lakes 
#############################################
basinmaker.watershed_delineation_without_lake_method(
    acc_thresold=500,
    mode="usingdem",
    max_memroy=1024 * 4,
    gis_platform="qgis",
)
```
### outputs  
* \<fdr_grass\> 
      is the flow direction dataset, which is using 1 - 8 to represent different directions, located in working folder 
* \<fdr_arcgis\> 
      is the flow direction dataset, which is using 1,2,4,...64,128 to represent different directions, located in working folder  
* \<str_v\> 
      is the generated river network in vector format, located in working folder 
* \<str_r\> 
      is the generated river network in raster format, located in working folder  
* \<cat_no_lake\> 
      is the delineated watersheds without considering lakes, located in working folder 
* \<acc\> 
      is the a flow accumulation dataset, located in working folder  


## add lakes and gauges control points into existing watershed delineation  
After obtain the watershed delineation without considering lakes. Here, the basinmaker will add lakes and gauge control points into existing watershed delineation.

Gauges will be snapped to the closest river network and then added as a catchment outlets.

each lake's inflow and outflow points will be identified and added as a catchment outlets. 

### parameters 

* \<path_lakefile_in\> 
      is a full path to the lake polygon file. user do not needs to do any preprocessing such as clipped and reproject. basinmaker will handle them. 
* \<lake_attributes\> 
      is a list of column names in the lake polygon file.lake_attributes[0] should be the column name represnt the unique lake id (integer). lake_attributes[1] should be the column name represent the lake type (integer). lake_attributes[2] should be the column name represent the lake area in km2 (float). lake_attributes[3] should be the column name represent the lake volumn in km3 (float). lake_attributes[4] should be the column name represent the lake depth in m (float).
* \<threshold_con_lake\> 
      is a lake area threshold in km2, all connected lakes with lake area smaller than this threshold will be removed 
* \<threshold_non_con_lake\> 
      is a lake area threshold in km2, all non connected lakes with lake area smaller than this threshold will be removed         
* \<path_obsfile_in\> 
      is a full path to the gauge point shp file. 
* \<obs_attributes\> 
      is a list of column names in the obs point file.obs_attributes[0] should be the column name represnt the unique obs id (integer). obs_attributes[1] should be the column name represent the obs name (string). obs_attributes[2] should be the column name represent the gauge drainage area in km2 (float). obs_attributes[3] should be the column name represent the source of the observation gauges,"CA" or "US" (character).      
* \<max_memroy\> 
      It is the maximum memory allow to be used by basinmaker  
* \<gis_platform\> 
      It is the parameter indicate which gis platform is used. For now only "qgis" is supported  
```
basinmaker.watershed_delineation_add_lake_control_points(
    path_lakefile_in=os.path.join(datafolder, "hylake.shp"),
    lake_attributes=["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    threshold_con_lake = 0,
    threshold_non_con_lake = 0,
    path_obsfile_in=os.path.join(datafolder, "obs.shp"),
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    max_memroy=1024 * 4,
    gis_platform="qgis",
)
```
### outputs  
* \<selected_lakes\> 
      is all selected lakes, located in working folder 
* \<sl_nonconnect_lake\> 
      is selected non connected lakes, located in working folder  
* \<sl_connected_lake\> 
      is selected connected lakes, located in working folder 
* \<river_without_merging_lakes\> 
      is the river network after add lake and gauge control points, located in working folder  
* \<catchment_without_merging_lakes\> 
      is catchment polygons after add lake and gauge control points, located in working folder 
* \<snapped_obs_points\> 
      is the gauge points after snapped to closest river network, located in working folder             

