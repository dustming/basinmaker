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

## Define project extent
The first step of watershed delineation is to define the project spatial extent. Once the project spatial extent is defined. The basinmaker will delineate lake-river routing structure within this extent. 

In this example, the project spatial extent will be defined using a wathershed outlet coordinates and a MERIT 90m DEM.

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


Before introducing the methodology and application procedure of the BasinMaker, the differences between routing structure without considering lake and routing structure with lakes defined by BasinMaker will be introduce in this section. 

A routing structure without considering lakes is showed in Figure 1A. Catchments defined by river reaches in this routing structure are showed in Figure 1A. In hydrological models, the streamflow is only explicitly simulated at the outlet of each catchment. It is not explicitly simulated at the mid of the river segment or any point inside each catchment. Thus, when using this routing structure to build hydrological model, we could obtain the simulated streamflow at end of each reach or catchment outlet in Figure 1A, but the inflow and the outflow of each lake showed in Figure 1 cannot be simulated. Thus, the impact of lake on the routing process modeling such as flow attenuation can not be correctly modeled. 

An example of routing structure with lakes generated by the BasinMaker is shown in Figure 1B. Lakes are divided into two categories by the BasinMaker: (1) connected lakes (CL), which indicates lakes directly connected with the river network; and (2) non-connected lakes (NCL), which denotes lakes not connected with the river network (Figure 1B). The resolution of the river network is defined by user provided flow accumulation threshold. Both CL and NCL defined here are within the drainage area of the watershed and water released from both CL and NCL will move to the outlet of the watershed. Thus, the definition of NCL here is not the same as the definition of none contributing area, in which runoff generated will not move to the watershed outlet.  

For both CL and NCL, the BasinMaker will use lake polygon to identify each lake’s inlets and outlet and then represent each lake as a lake catchment (LC) shown Figure 1B. A lake catchment will fully cover the lake polygon and its outlet is the same as the outlet of the lake (Figure 1B). At the same time, each inlet of the CL will also be identified as a catchment outlet (Figure 1B). In this way, both inflow and outflow of each lake can be explicitly simulated by the semi-distributed hydrological models. 

<figure>
    <p align="center">
    <img src="https://github.com/dustming/RoutingTool/wiki/Figures/Figure1.png" width="100%" height="100%" />
    </p>
    <font size="1">
    <figcaption width="50%"> <b>Figure 1</b>: Lakes in the generated routing network by the lake river routing toolbox. A is the predefined river network and catchment boundary (Catchment boundary (River)) defined by the predefined river network without considering lakes. Figure B, the generated lake river routing structure by this toolbox using predefined river network and lake's polygons.<br>
    </figcaption>
    </font>
</figure>

# Documentation 

The BasinMaker is developed within python3 environment and using several basic raster and vector functions in QGIS and GRASS GIS. The list of functions from GRASS and QGIS and have been used by BasinMaker can be found in [here](https://github.com/dustming/RoutingTool/wiki/Installation-of-the-BasinMaker#functions-from-qgis-and-grass-gis) 

* The installation of the toolbox can be found in [here](https://github.com/dustming/RoutingTool/wiki/Installation-of-the-BasinMaker)

* List of all tools provided by BasinMaker can be found in [here](https://github.com/dustming/RoutingTool/wiki/List-of-tools-provided-by-BasinMaker)

# Example usage

* The application procedure to delineate a lake river routing structure from DEM can be found in [here](https://github.com/dustming/RoutingTool/wiki/Application-procedure-for-watershed-delineation-with-lakes)

* The application procedure to simplify an existing routing network can be found in [here](https://github.com/dustming/RoutingTool/wiki/Application-procedure-for-post-processing-tools)

* The application procedure to working with large domain with high resolution DEM [here]()


# Citation
Han, M., Mai, J., Tolson, B. A., Craig, J. R., Gaborit, É., Liu, H., and Lee, K. (2020a): Subwatershed-based lake and river routing products for hydrologic and land surface models applied over Canada, Canadian Water Resources Journal, 0, 1-15. ([publication](https://www.tandfonline.com/doi/ref/10.1080/07011784.2020.1772116?scroll=top))

Han, M. et al. (2020b): An automated GIS toolbox for watershed delineation with lakes
In preparation.

Han, M., Mai, J., Tolson, B. A., Craig, J. R., Gaborit, É., Liu, H., and Lee, K. (2020c): A catchment-based lake and river routing product for hydrologic and land surface models in Canada (Dataset) Zenodo. ([dataset](https://zenodo.org/record/3667677#.X7xD0c1KiUk))

# Acknowledgments
 
Funded under Lake Futures project of Global Water Futures program.
<p align="left">
   <img src="https://github.com/dustming/Crop_Growth_Library/wiki/figures/logos/GWF_logo.png" width="10%"  />
   <img src="https://github.com/dustming/Crop_Growth_Library/wiki/figures/logos/LAKEFUTURES_logo.png" width="10%"  />
</p>
