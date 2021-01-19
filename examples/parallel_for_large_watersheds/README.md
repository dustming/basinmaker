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

## example usage 
```
python Parallel_Step_1_Define_Sub_Region.py
```
## Explanation of script

### parameters
* \<path_lakefile_in\> 
      is a full path to the lake polygon file. user do not needs to do any preprocessing such as clipped and reproject. basinmaker will handle them. 

* \<lake_attributes\> 
      is a list of column names in the lake polygon file.lake_attributes[0] should be the column name represnt the unique lake id (integer). lake_attributes[1] should be the column name represent the lake type (integer). lake_attributes[2] should be the column name represent the lake area in km2 (float). lake_attributes[3] should be the column name represent the lake volumn in km3 (float). lake_attributes[4] should be the column name represent the lake depth in m (float).
      
* \<path_bkfwidthdepth\> 
      is a full path to the bankfull width and depth dataset. user do not needs to do any preprocessing such as clipped and reproject. basinmaker will handle them. the dataset in this example can be found [here](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/wrcr.20440) 

* \<bkfwd_attributes\> 
      is a list of column names in the bankfull width and depth file. bkfwd_attributes[0] should be the column name represnt the width of the channel in m (float). bkfwd_attributes[1] should be the column name represent depth of channel in m (float). bkfwd_attributes[2] should be the column name represent the annual mean flow of the channel in m3/s (float). bkfwd_attributes[3] should be the column name represent the drainage area of the channel km2 (float).
      
* \<Min_Num_Domain\> 
      is the minimum number of sub region 

* \<Max_Num_Domain\> 
      is the maximum number of sub region  

* \<Initaial_Acc\> 
      is the initial flow accumulation threshold to divide domain into sub regions 

* \<Delta_Acc\> 
      is the change of flow accumulation threshold if the number of sub region is largher or smaller than Max_Num_Domain and Min_Num_Domain

* \<CheckLakeArea\> 
      is the lake area threshold to determine subregion. lake with lake area smaller than this value will not be considered during define sub region process.  

* \<Acc_Thresthold_stream\> 
      is the flow accumulation threshold will be used to delineate watershed within each sub region 

* \<max_memory\> 
       is the max memeory allow to be used     
      
* \<Out_Sub_Reg_Folder\> 
      is the output folder to save generated sub regions              
      
```
# divide domain in to subregions 
basinmaker.divide_domain_into_sub_regions_method(
    path_lakefile_in = in_lake,
    lake_attributes = ["Hylak_id", "Lake_type", "Lake_area", "Vol_total", "Depth_avg"],
    path_bkfwidthdepth=os.path.join(datafolder, "bkf_wd.shp"),
    bkfwd_attributes=["WIDTH", "DEPTH", "Q_Mean", "UP_AREA"],
    Min_Num_Domain=9,
    Max_Num_Domain=30,
    Initaial_Acc=250000,
    Delta_Acc=20000,
    CheckLakeArea=10,
    fdr_path = '#',
    Acc_Thresthold_stream=2000,
    max_memory=2048*3,
    Out_Sub_Reg_Folder=Out_Sub_Reg_Folder,
)
```

### outputs  
* \<HyMask_region_xxxxx.shp\> 
      it is the buffered sub region mask   

* \<HyMask_region_xxxxx_nobuffer.shp\> 
      is the sub region make without buffer 

* \<outlet_pt_info.shp\> 
      is the sub region outlet point   

* \<sub_reg_inlet.shp\> 
      is the sub region inlet  

* \<sub_reg_nfdr_arcgis.pack\> 
      is the updated flow direction in 1,2,4,8...              

* \<sub_reg_nfdr_grass.pack\> 
      is updated flow direction in 1,2,3,4...   

* \<Sub_Reg_Outlet_v.pack\> 
      is the packed sub region outlet   

* \<sub_reg_str_r.pack\> 
      is the stream raster generated by Acc_Thresthold_stream   

* \<sub_reg_str_v.pack\> 
      is the sub region inlet  

* \<sub_reg_nfdr_arcgis.pack\> 
      is the stream raster generated by Acc_Thresthold_stream              



# parallel generate lake-river routing structure in each sub region 

In this example, for cores will be used to parallel process sub regions. number of cores can be modified the variable 'ncores' in "Parallel_Step_2_Main_Function.py" . 
The explanation of the script in "" can be found in [here](https://github.com/dustming/basinmaker/tree/master/examples/wathershed_delineation_merit_90m_dem).  

## example usage 
```
python Parallel_Step_2_Main_Function.py
```



