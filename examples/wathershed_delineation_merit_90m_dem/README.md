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
           
## delineate watershed without considering lakes 
After obtain the project extent. 

In this example, the project spatial extent will be defined using a wathershed outlet coordinates and a MERIT 90m DEM.

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
# generate a watershed delineation without considering lakes 
#############################################
basinmaker.watershed_delineation_without_lake_method(
    acc_thresold=500,
    mode="usingdem",
    max_memroy=1024 * 4,
    gis_platform="qgis",
)
```


