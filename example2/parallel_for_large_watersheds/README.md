# Table of content 
1. [Overview]()

2. [Divide Domain into sub-regions]()

3. [Parallel define lake river routing structure for each sub-region]()

4. [Combine lake river routing structure of each sub-region]()

# Overview

When we working with an extremely large region and the resolution of the DEM is very high. For example, we are trying to delineate the lake river routing structure including all lakes for the entire great lake watershed with a DEM resolution of 90 m. The procedure showed in here can help reduce the processing time and set up a workflow when part of the routing structure needs to be adjusted. It is very common that we may find the routing network needs to be adjusted. For example, the location of the gauge station is incorrect or part of the river network is incorrect. Instead of re-delineate the routing structure for the whole domain, the procedure in this section allows us to only adjust the region where it goes wrong. 

The procedure can be divided into three steps: The first step is to divide the whole domain into several smaller subregions; the second step is to define the lake river routing structure for each sub-region in a parallel way; Finally, the lake river routing network for the whole domain can be defined by combine lake river routing structure of each sub-region. 

In the following example, when we apply this approach to delineate lake river routing for a region with 22243473 grids are following:

- Define sub-region takes 8.5 mins  
- Generate lake river routing structure for each sub-region, take 56 mins 
- Combine lake river routing structure for each sub-region, take 1.2 mins 

Please note that the processing time in step two include: 1) It includes the time used to clip, reproject and rasterize global input products (e.g. HydroLakes polygons) into the sub-region extent; 2) It includes the time define lake-river routing structure; and It includes the time used to estimate routing parameters for each catchments in the lake-river routing structure.  
  
# Divide Domain into sub-regions
