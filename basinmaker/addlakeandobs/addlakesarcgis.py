from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *


def add_lakes_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_lakefile_in,
    lake_attributes,
    threshold_con_lake,
    threshold_non_con_lake,
    only_included_lake_at_river_interction = False,
    remove_lake_inlets = False,
    path_sub_reg_lake_r="#",
    path_sub_reg_lake_bd_r="#",
    sl_connected_lake="sl_connected_lake",
    sl_non_connected_lake="sl_nonconnect_lake",
    sl_lakes="selected_lakes",
    sl_str_connected_lake="str_sl_connected_lake",
    nfdr_arcgis="narcgis_fdr",
    nfdr_grass="ngrass_fdr",
    cat_add_lake="cat_add_lake",
    pourpoints_with_lakes="pourpoints_with_lakes",
    cat_use_default_acc="cat_use_default_acc",
    lake_outflow_pourpoints="lake_outflow_pourpoints",
    problem_seg="problem_seg",
    max_memroy=1024 * 4,
):
    work_folder = grassdb
    # define required input files names
    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]

    # define internal file names
    lake_inflow_pourpoints = Internal_Constant_Names["lake_inflow_pourpoints"]
    catchment_pourpoints_outside_lake = Internal_Constant_Names[
        "catchment_pourpoints_outside_lake"
    ]
    cat_add_lake_old_fdr = Internal_Constant_Names["cat_add_lake_old_fdr"]
    str_connected_lake = Internal_Constant_Names["str_connected_lake"]
    alllake = Internal_Constant_Names["all_lakes"]
    lake_boundary = Internal_Constant_Names["lake_boundary"]
    connected_lake = Internal_Constant_Names["connect_lake"]
    non_connected_lake = Internal_Constant_Names["nonconnect_lake"]
    lakes_lg_cl_thres = "lakes_lg_cl_thres"
    lakes_lg_ncl_thres = "lakes_lg_ncl_thres"


    if not os.path.exists(work_folder):
        os.makedirs(work_folder)
    arcpy.env.workspace = work_folder
     
    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    cellSize = float(arcpy.GetRasterProperties_management(dem+'.tif', "CELLSIZEX").getOutput(0))
    SptailRef = arcpy.Describe(dem+'.tif').spatialReference
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize    
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(dem + '.tif').extent
    
    ## processing lakes move to a function later 
    
    ## obtain lakes fully contained in the mask region 
    ## any lake overlay with the mask boundary will be removed 
    arcpy.Project_management(
        path_lakefile_in,
        "lake_proj.shp", 
        arcpy.SpatialReference(int(SptailRef.factoryCode)),
        )  
        
    arcpy.FeatureToLine_management(mask + '.shp', mask + '_line.shp')
    arcpy.Clip_analysis("lake_proj.shp", mask + '.shp', "lake_clip.shp", "")
    arcpy.Intersect_analysis(["lake_proj.shp",mask + '_line.shp'], 'lake_inter.shp')
    inter_lake = pd.DataFrame.spatial.from_featureclass(os.path.join(work_folder,'lake_inter.shp'))
    all_cliped_lakes = pd.DataFrame.spatial.from_featureclass(os.path.join(work_folder,'lake_clip.shp'))
    lakeids_inter = inter_lake[lake_attributes[0]].values
    select_lake = all_cliped_lakes.loc[~all_cliped_lakes[lake_attributes[0]].isin(lakeids_inter)]
    select_lake.spatial.to_featureclass(location=os.path.join(work_folder,alllake + ".shp"),overwrite=True,sanitize_columns=False)    
    ###
    arcpy.FeatureToLine_management(alllake + '.shp', alllake + '_line.shp')
    
    arcpy.Dissolve_management(alllake + '_line.shp', lake_boundary + ".shp", [lake_attributes[0]])
    
    arcpy.PolygonToRaster_conversion(alllake + ".shp", lake_attributes[0], alllake + ".tif",
                                 "MAXIMUM_COMBINED_AREA",lake_attributes[0], cellSize)
                                 
    arcpy.PolylineToRaster_conversion(lake_boundary + ".shp", lake_attributes[0], lake_boundary + ".tif",
                                 "MAXIMUM_COMBINED_LENGTH",lake_attributes[0], cellSize)
                                 
    return
