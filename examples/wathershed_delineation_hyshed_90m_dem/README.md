# Example - delineate lake river routing structure with HydroSHEDS 90 m DEM
 
# Overview
In this example, the HydroSHEDS 90 m DEM will be used to delineate the lake river routing structure. The extent of the watershed boundary will be estimated by an provided outlet coordinates. The basinmaker is called in the in the example python script "example_watershed_delineation_with_lakes_merit_90m.py". After setup required basinmaker working environments (instruction can be found at [here](https://github.com/dustming/basinmaker/wiki/Installation-of-the-BasinMaker)). This script can be run by following commmand:
```
python example_watershed_delineation_with_lakes_hyshed_90m.py
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
      is a parameter indicate which method will be used to define project extent. In this case, "using_hybasin" is used. 
* \<path_dem_in\> 
      is a full path to the DEM dataset. 
* \<hybasin_ply\> 
      is the to the HydroBASINS products
* \<down_hybasin_id\> 
      is the most downstream hydroBASIN ID in the region of interest.
* \<buffer_distance\> 
      is the parameter that used to buffer the extent defined by \<down_hybasin_id\>      
* \<gis_platform\> 
    It is the parameter indicate which gis platform is used. For now only "qgis" is supported   
```
#############################################
# define extent of the processing domain  
# please linkt to your local hydrbasin level 07 
# product
#############################################
basinmaker.define_project_extent_method(
    mode="using_hybasin",
    path_dem_in=os.path.join(datafolder, "hyshed_90_dem.tif"),
    buffer_distance=0.01,
    hybasin_ply="C:/Users/dustm/OneDrive - University of Waterloo/Documents/ProjectData/Geo_Data_Base/Shapefiles/HydroBASINS/hybas_na_lev07_v1c.shp",
    down_hybasin_id=7070356140,
    gis_platform="qgis",
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


## add hydrological attributes to catchments without merging lakes   
After adding lakes and gauge control points into existing watershed delineation. Here, the basinmaker will add hydrological attributes for each generated catchment.

### parameters 

* \<path_bkfwidthdepth\> 
      is a full path to the bankfull width and depth dataset. user do not needs to do any preprocessing such as clipped and reproject. basinmaker will handle them. the dataset in this example can be found [here](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/wrcr.20440) 
* \<bkfwd_attributes\> 
      is a list of column names in the bankfull width and depth file. bkfwd_attributes[0] should be the column name represnt the width of the channel in m (float). bkfwd_attributes[1] should be the column name represent depth of channel in m (float). bkfwd_attributes[2] should be the column name represent the annual mean flow of the channel in m3/s (float). bkfwd_attributes[3] should be the column name represent the drainage area of the channel km2 (float).
* \<path_landuse\> 
      is the full path of the table in '.csv' format.The table describe the floodplain roughness coefficient correspond to a given landuse type. The table should have two columns: RasterV and MannV. RasterV is the landuse value in the landuse raster for each land use type and the MannV is the roughness coefficient value for each landuse type. 
* \<path_landuse_info\> 
      is a lake area threshold in km2, all non connected lakes with lake area smaller than this threshold will be removed  
* \<lake_attributes\> 
      is a list of column names in the lake polygon file.lake_attributes[0] should be the column name represnt the unique lake id (integer). lake_attributes[1] should be the column name represent the lake type (integer). lake_attributes[2] should be the column name represent the lake area in km2 (float). lake_attributes[3] should be the column name represent the lake volumn in km3 (float). lake_attributes[4] should be the column name represent the lake depth in m (float).       
* \<obs_attributes\> 
      is a list of column names in the obs point file. obs_attributes[0] should be the column name represnt the unique obs id (integer). obs_attributes[1] should be the column name represent the obs name (string). obs_attributes[2] should be the column name represent the gauge drainage area in km2 (float). obs_attributes[3] should be the column name represent the source of the observation gauges,"CA" or "US" (character).      
* \<outlet_obs_id\> 
      is the one gauge id from obs_attributes[0] in the provided gauge file. if it is larger than zero. basinmaker will remove the catchment do not drainage to this gauge 
* \<path_sub_reg_outlets_v\> 
      is the path to the sub region outlet file. it is only needed when basinmaker has divide the project extent into several sub regions. and then process each sub region in a parallel approach.     
* \<gis_platform\> 
      It is the parameter indicate which gis platform is used. For now only "qgis" is supported  
* \<output_folder\> 
      is the path to the output folder  
```
#############################################
# add hydrological attributes to existing watershed delineation  
#############################################
basinmaker.add_attributes_to_catchments_method(
    path_bkfwidthdepth=os.path.join(datafolder, "bkf_wd.shp"),
    bkfwd_attributes=["WIDTH", "DEPTH", "Q_Mean", "UP_AREA"],
    path_landuse=os.path.join(datafolder, "landuse_modis_250.tif"),
    path_landuse_info=os.path.join(datafolder, "Landuse_info3.csv"),
    gis_platform="qgis",
    obs_attributes=["Obs_ID", "STATION_NU", "DA_obs", "SRC_obs"],
    lake_attributes =["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"] ,
    outlet_obs_id=-1,
    path_sub_reg_outlets_v="#",
    output_folder=path_output_folder,
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