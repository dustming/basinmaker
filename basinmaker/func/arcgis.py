from arcgis import GIS
from arcgis.features import GeoAccessor, GeoSeriesAccessor
import arcpy
from arcpy import env
from arcpy.sa import *
import numpy as np
import os
import pandas as pd
import tempfile
import arcpy.cartography as CA
from json import load, JSONEncoder
import json
import shutil
import copy
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
pd.options.mode.chained_assignment = None
#####
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")

def select_feature_by_attributes_arcgis(input,Attri_NM,Attri_v,output):
    where_clause = '"%s" IN' % (Attri_NM)
    where_clause = where_clause + " ("
    for i in range(0,len(Attri_v)):
        if i == 0:
            where_clause = where_clause + str(Attri_v[i])
        else:
            where_clause = where_clause + "," + str(Attri_v[i])
    where_clause = where_clause + ")"

    arcpy.Select_analysis(input, output, where_clause)
    return
##################


def Remove_Unselected_Lake_Attribute_In_Finalcatinfo_Arcgis(finalcat_ply, Conn_Lake_Ids):
    """Functions will set lake id not in Conn_Lake_Ids to -1.2345 in attribute
        table of Path_Finalcatinfo
    ----------

    Notes
    -------

    Returns:
    -------
        None, the attribute table of Path_shpfile will be updated
    """

    mask1 = np.logical_not(finalcat_ply['HyLakeId'].isin(Conn_Lake_Ids))
    mask2 = finalcat_ply['Lake_Cat'] != 2
    mask = np.logical_and(mask1,mask2)

    finalcat_ply.loc[mask,'HyLakeId'] = 0
    finalcat_ply.loc[mask,'LakeVol'] = 0
    finalcat_ply.loc[mask,'LakeArea'] = 0
    finalcat_ply.loc[mask,'LakeDepth'] = 0
    finalcat_ply.loc[mask,'Laketype'] =0
    finalcat_ply.loc[mask,'Lake_Cat'] = 0

    return finalcat_ply


def save_modified_attributes_to_outputs(mapoldnew_info,tempfolder,OutputFolder,cat_name,riv_name,Path_final_riv,dis_col_name='SubId'):

    mapoldnew_info.spatial.to_featureclass(location=os.path.join(tempfolder,'updateattri.shp'),overwrite=True,sanitize_columns=False)
    arcpy.Dissolve_management(os.path.join(tempfolder,'updateattri.shp'), os.path.join(OutputFolder,cat_name), [dis_col_name])
    arcpy.JoinField_management(os.path.join(OutputFolder,cat_name), dis_col_name, os.path.join(tempfolder,'updateattri.shp'), dis_col_name)
    arcpy.DeleteField_management(os.path.join(OutputFolder,cat_name),
        ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters"]
    )


    if riv_name != '#':
        arcpy.CalculateGeometryAttributes_management(os.path.join(OutputFolder, cat_name), [["centroid_x", "CENTROID_X"], ["centroid_y", "CENTROID_Y"]],coordinate_system = arcpy.SpatialReference(4326))
        cat_colnms = mapoldnew_info.columns
        drop_cat_colnms = cat_colnms[cat_colnms.isin(["SHAPE","SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters","Old_DowSubId"])]
        cat_pd = mapoldnew_info.drop(columns=drop_cat_colnms)

        if Path_final_riv != '#':
            riv_pd = pd.DataFrame.spatial.from_featureclass(Path_final_riv)
            riv_pd['Old_SubId'] = riv_pd['SubId']
            # remove all columns
            riv_pd = riv_pd[['SHAPE','Old_SubId']]
            riv_pd = pd.merge(riv_pd, cat_pd, on='Old_SubId', how='left')

            riv_pd = riv_pd.drop(columns=['Old_SubId'])
            mask = np.logical_or(riv_pd['RivLength'] > 0,riv_pd['Lake_Cat'] > 0)
            riv_pd = riv_pd[mask]
            if len(riv_pd) > 0:
                riv_pd.drop(columns=['centroid_x','centroid_y'])
                riv_pd.spatial.to_featureclass(location=os.path.join(tempfolder,'riv_attri.shp'),overwrite=True,sanitize_columns=False)
                arcpy.Dissolve_management(os.path.join(tempfolder,'riv_attri.shp'), os.path.join(OutputFolder,riv_name), ["SubId"])
                arcpy.JoinField_management(os.path.join(OutputFolder,riv_name), "SubId", os.path.join(tempfolder,'riv_attri.shp'), "SubId")
                arcpy.DeleteField_management(os.path.join(OutputFolder,riv_name),
                    ["SubId_1", "Id","nsubid2", "nsubid","ndownsubid","Old_SubId","Old_DowSub","Join_Count","TARGET_FID","Id","SubID_Oldr","HRU_ID_N_1","HRU_ID_N_2","facters",'centroid_x','centroid_y']
                )

    # if "finalcat_info" in cat_name:
    #     create_geo_jason_file(os.path.join(OutputFolder,cat_name))

def clean_attribute_name_arcgis(table,names):
    remove_column_names = table.columns[np.logical_not(np.isin(table.columns,names))]
    table = table.drop(columns=remove_column_names)
    return table


####
def create_geo_jason_file(Input_Polygon_path):

    arcpy.env.overwriteOutput = True
    product_dir = os.path.dirname(Input_Polygon_path)
    Names_in = os.path.basename(Input_Polygon_path).split('_')
    n_charc = len(Names_in)
    version  = Names_in[n_charc - 1][0:4]
    TOLERANCEs = [0.0001,0.0005,0.001,0.005,0.01,0.05]


    head_name_cat = "finalcat_info"
    head_name_riv = "finalcat_info_riv"
    head_name_slake = "sl_connected_lake"
    head_name_nlake = "sl_non_connected_lake"

    Input_file_name = []
    Output_file_name = []
    if 'v' in version:
        Input_file_name = [
                           head_name_cat + "_"+version+'.shp',
                           head_name_riv + "_"+version+'.shp',
                           head_name_slake + "_"+version+'.shp',
                           head_name_nlake + "_"+version+'.shp',
                           ]
        Output_file_name = [
                           head_name_cat + "_"+version+'.geojson',
                           head_name_riv + "_"+version+'.geojson',
                           head_name_slake + "_"+version+'.geojson',
                           head_name_nlake + "_"+version+'.geojson',
                           ]
    else:
        Input_file_name = [
                           head_name_cat +'.shp',
                           head_name_riv +'.shp',
                           head_name_slake +'.shp',
                           head_name_nlake +'.shp',
                           ]
        Output_file_name = [
                           head_name_cat +'.geojson',
                           head_name_riv +'.geojson',
                           head_name_slake +'.geojson',
                           head_name_nlake +'.geojson',
                           ]
    created_jason_files = []
    created_jason_files_lake_riv = []

    for i  in  range(0,len(Input_file_name)):
        input_path = os.path.join(product_dir,Input_file_name[i])
        output_jason_path = os.path.join(product_dir,Output_file_name[i])
        if not os.path.exists(input_path):
            continue
        created_jason_files.append(output_jason_path)

        if 'finalcat_info_riv' in Input_file_name[i] or 'connected_lake' in Input_file_name[i]:
            created_jason_files_lake_riv.append(output_jason_path)

        # reproject to WGS84
        input_wgs_84 = os.path.join(tempfile.gettempdir(),"input_wgs_84.shp")
        arcpy.Project_management(input_path, input_wgs_84, arcpy.SpatialReference(int(4326)))


        if 'finalcat_info' in Input_file_name[i] or "finalcat_info_riv" in Input_file_name[i]:
            arcpy.AddField_management(input_wgs_84, 'rvhName', "TEXT")
            arcpy.CalculateField_management(input_wgs_84, 'rvhName', "'sub' + str(int(\"!SubId!\"))", "PYTHON3")

        arcpy.RepairGeometry_management(input_wgs_84)

        for TOLERANCE in TOLERANCEs:
            input_wgs_84_simplify = os.path.join(tempfile.gettempdir(),"input_wgs_84_simplify.shp")
            if arcpy.Exists(input_wgs_84_simplify):
                arcpy.Delete_management(input_wgs_84_simplify)


            if "finalcat_info_riv" not in Input_file_name[i]:
                CA.SimplifyPolygon(input_wgs_84, input_wgs_84_simplify, "POINT_REMOVE", TOLERANCE)
            else:
                CA.SimplifyLine(input_wgs_84, input_wgs_84_simplify, "POINT_REMOVE", TOLERANCE)

            arcpy.FeaturesToJSON_conversion(input_wgs_84_simplify,
                                            output_jason_path,"FORMATTED","NO_Z_VALUES","NO_M_VALUES","GEOJSON")

            json_file_size = os.stat(output_jason_path).st_size/1024/1024 #to MB
            if json_file_size <= 100:
                break

    if len(created_jason_files_lake_riv) > 1 and os.stat(os.path.join(product_dir,Output_file_name[1])).st_size/1024/1024 < 500:
        for i in range(0,len(created_jason_files_lake_riv)):
            injson2 = load(open(created_jason_files_lake_riv[i]))
            if 'finalcat_info_riv' in created_jason_files_lake_riv[i]:
                new_features = []
                for element in injson2["features"]:
                    if element["properties"]["Lake_Cat"] == 0:
                        new_features.append(element)
                injson2["features"] = new_features

            if i == 0:
                output_jason_lake_riv = injson2
            else:
                output_jason_lake_riv['features'] += injson2['features']

        with open(os.path.join(product_dir,'routing_product_lake_river.geojson'), 'w', encoding='utf-8') as f:
            json.dump(output_jason_lake_riv, f, ensure_ascii=False, indent=4)
    else:
        shutil.copy(created_jason_files_lake_riv[0], os.path.join(product_dir,'routing_product_lake_river.geojson'))



    return

def read_table_as_pandas(table,col_names,work_folder):
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arr = arcpy.da.TableToNumPyArray(table, col_names)
    df = pd.DataFrame(arr)
    return df


def find_zonal_maximum_or_minimum_point_on_the_raster(zonal_raster,val_raster,outname,statistics_type,cellSize,SptailRef,work_folder,snapRaster):
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arcpy.env.overwriteOutput = True
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(snapRaster).extent
    arcpy.env.snapRaster =  snapRaster

    maxacc_value_in_zone = ZonalStatistics(zonal_raster,"VALUE",val_raster,
                                     statistics_type, "NODATA")
    outlet_points = Con(Abs(val_raster - maxacc_value_in_zone) < 0.000001,zonal_raster)
    outlet_points.save(outname+"_r")
    arcpy.RasterToPoint_conversion(outname+"_r", outname + "_v", "VALUE")
    outlets =  pd.DataFrame.spatial.from_featureclass(outname + "_v")

    # check if each zone only has one outlet points
    dupli_mask = outlets["grid_code"].duplicated()
    outlets_dup = outlets[dupli_mask]
    if len(outlets_dup) == 0:
        ExtractMultiValuesToPoints(outname + "_v", [[val_raster,val_raster]], "NONE")
        return outname+"_r",outname + "v"
    else:
        with open(os.path.join(work_folder,'log.txt'), 'a') as logfile:
            logfile.write("Check SubId and DownSubId of following subbasins \n")
            for item in outlets_dup['grid_code'].values:
                logfile.write("%i\n" % item)
        # remove duplicated outlet points based on

        ExtractMultiValuesToPoints(outname + "_v", [[val_raster,val_raster]], "NONE")

        outlets_n =  pd.DataFrame.spatial.from_featureclass(outname + "_v")
        outlets_n = outlets_n.sort_values(by=val_raster, ascending=True)
        if statistics_type == "MAXIMUM":
            outlets_n = outlets_n.drop_duplicates(subset=['grid_code'], keep='last')
        else:
            outlets_n = outlets_n.drop_duplicates(subset=['grid_code'], keep='first')
        outlets_n.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb",outname + "_v"),overwrite=True,sanitize_columns=False)
        arcpy.PointToRaster_conversion(outname + "_v", 'grid_code', outname + "_r")


        return outname + "_r",outname + "_v"


def pre_process_lake_polygon(path_lakefile_in,alllake,lake_attributes,lake_boundary,mask,cellSize,SptailRef,work_folder,snapRaster):

    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arcpy.env.overwriteOutput = True
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(snapRaster).extent
    arcpy.env.snapRaster =  snapRaster

    arcpy.Project_management(
        path_lakefile_in,
        "lake_proj",
        arcpy.SpatialReference(int(SptailRef.factoryCode)),
        )
    arcpy.FeatureToLine_management(mask + '_ply', mask + '_line')
    arcpy.Clip_analysis("lake_proj", mask + '_ply', "lake_clip", "")

    arcpy.analysis.Buffer(in_features = "lake_clip", out_feature_class = "lake_clip_bf",buffer_distance_or_field = str(cellSize*1/4) + "  Meters",method="GEODESIC")

    arcpy.Intersect_analysis(["lake_clip_bf",mask + '_line'], 'lake_inter')
    inter_lake = pd.DataFrame.spatial.from_featureclass('lake_inter')
    all_cliped_lakes = pd.DataFrame.spatial.from_featureclass('lake_clip_bf')
    lakeids_inter = np.unique(inter_lake[lake_attributes[0]].values)
    select_lake = all_cliped_lakes[~all_cliped_lakes[lake_attributes[0]].isin(lakeids_inter)]
    select_lake.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb",alllake+"_v"),overwrite=True,sanitize_columns=False)
    arcpy.FeatureToLine_management(alllake+"_v", lake_boundary+"_v","0.001 Meters")
#    arcpy.Dissolve_management(alllake + '_line', lake_boundary+"_v", [lake_attributes[0]])
    arcpy.PolygonToRaster_conversion(alllake+"_v", lake_attributes[0], alllake+"_rt",
                                 "MAXIMUM_COMBINED_AREA",None, snapRaster)
    arcpy.PolylineToRaster_conversion(lake_boundary+"_v", lake_attributes[0], lake_boundary+"_r",
                                 "MAXIMUM_COMBINED_LENGTH",None, snapRaster)

    OutRas = Con(IsNull(lake_boundary+"_r"), alllake+"_rt", lake_boundary+"_r")
    OutRas.save(alllake+"_r")

    lake_at_lake_boundary = all_cliped_lakes[all_cliped_lakes[lake_attributes[0]].isin(lakeids_inter)]
    lake_at_lake_boundary.spatial.to_featureclass(location=os.path.join(work_folder,"lakes_at_lake_watershed_boundary.shp"),overwrite=True,sanitize_columns=False)
#    arcpy.management.Delete(r"'lake_proj';'lake_clip_bf';'lake_clip';'lake_inter'")
    return


def define_cl_and_ncl_lakes(str_r,alllake,str_connected_lake,sl_connected_lake,sl_non_connected_lake,sl_lakes,lake_attributes,threshold_con_lake,threshold_non_con_lake,cellSize,SptailRef,work_folder,snapRaster):

    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arcpy.env.overwriteOutput = True
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(snapRaster).extent
    arcpy.env.snapRaster =  snapRaster
    lowerLeft = arcpy.Point(arcpy.Describe(snapRaster).extent.XMin,arcpy.Describe(snapRaster).extent.YMin)

    str_connected_lake_var = Con(IsNull(str_r),-9999,alllake+"_r")
    str_connected_lake_var.save(str_connected_lake)

    str_connected_lake_arr = arcpy.RasterToNumPyArray(str_connected_lake_var,nodata_to_value=-9999)
    all_lake_arr = arcpy.RasterToNumPyArray(alllake+"_r",nodata_to_value=-9999)
    Connect_Lake_Ids = np.unique(str_connected_lake_arr)
    Connect_Lake_Ids = Connect_Lake_Ids[Connect_Lake_Ids!=-9999]
    cl_lake_arr = copy.deepcopy(all_lake_arr)
    ncl_lake_arr = copy.deepcopy(all_lake_arr)
    sl_lake_arr = copy.deepcopy(all_lake_arr)
    mask = np.isin(all_lake_arr, Connect_Lake_Ids)
    ncl_lake_arr[mask] = -9999
    cl_lake_arr[~mask] = -9999

    all_lakes_table = pd.DataFrame.spatial.from_featureclass(alllake + "_v")
    un_sl_con_lakeids = all_lakes_table[
        all_lakes_table[lake_attributes[2]] < threshold_con_lake
    ][lake_attributes[0]].values
    un_sl_noncon_lakeids = all_lakes_table[
        all_lakes_table[lake_attributes[2]] < threshold_non_con_lake
    ][lake_attributes[0]].values

    if len(un_sl_con_lakeids) > 0:
        mask = np.isin(cl_lake_arr, un_sl_con_lakeids)
        cl_lake_arr[mask] = -9999
    if len(un_sl_noncon_lakeids) > 0:
        mask = np.isin(ncl_lake_arr, un_sl_noncon_lakeids)
        ncl_lake_arr[mask] = -9999

    ncl_lake_r = arcpy.NumPyArrayToRaster(ncl_lake_arr,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    ncl_lake_r.save(sl_non_connected_lake+"_r")
    cl_lake_r = arcpy.NumPyArrayToRaster(cl_lake_arr,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    cl_lake_r.save(sl_connected_lake+"_r")

    cl_lake_ids = np.unique(cl_lake_arr)
    ncl_lake_ids = np.unique(ncl_lake_arr)

    sl_lake_ids = np.concatenate((cl_lake_ids, ncl_lake_ids))
    sl_lake_ids = sl_lake_ids[sl_lake_ids!=-9999]
    mask = np.isin(sl_lake_arr, sl_lake_ids)
    sl_lake_arr[~mask] = -9999
    sl_lake_r = arcpy.NumPyArrayToRaster(sl_lake_arr,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    sl_lake_r.save(sl_lakes+"_r")
    cl_lake_v = all_lakes_table[all_lakes_table[lake_attributes[0]].isin(cl_lake_ids)].copy(deep=True)
    cl_lake_v.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb",sl_connected_lake+"_v"),overwrite=True,sanitize_columns=False)

    ncl_lake_v = all_lakes_table[all_lakes_table[lake_attributes[0]].isin(ncl_lake_ids)].copy(deep=True)
    ncl_lake_v.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb",sl_non_connected_lake+"_v"),overwrite=True,sanitize_columns=False)

    return


def return_routing_info_using_str_v(str_v,cellSize,SptailRef,work_folder,snapRaster):
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arcpy.env.overwriteOutput = True
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(snapRaster).extent
    arcpy.env.snapRaster =  snapRaster
    lowerLeft = arcpy.Point(arcpy.Describe(snapRaster).extent.XMin,arcpy.Describe(snapRaster).extent.YMin)

    str_v_pd = pd.DataFrame.spatial.from_featureclass(str_v)
    str_v_pd["SubId"] = str_v_pd["grid_code"]
    str_v_pd["DowSubId"] = np.nan
    str_v_pd["n_up_sub"] = np.nan
    for idx in str_v_pd.index:
        to_node = str_v_pd.loc[idx,"to_node"]
        from_node = str_v_pd.loc[idx,"from_node"]
        dowsubinfo = str_v_pd[str_v_pd["from_node"] == to_node]
        upsubinfo = str_v_pd[str_v_pd["to_node"] == from_node]
        if len(dowsubinfo) == 1:
            if dowsubinfo["grid_code"].values[0] == str_v_pd.loc[idx,"SubId"]:
               str_v_pd.loc[idx,"DowSubId"] = -1.2345
            else:
                str_v_pd.loc[idx,"DowSubId"] = dowsubinfo["grid_code"].values[0]
        else:
            str_v_pd.loc[idx,"DowSubId"] = -1

        if len(upsubinfo) > 0:
            str_v_pd.loc[idx,"n_up_sub"] = len(dowsubinfo)
        else:
            str_v_pd.loc[idx,"n_up_sub"] = 0

    str_v_pd.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb",str_v + "_rout"),overwrite=True,sanitize_columns=False)
    str_v_pd = str_v_pd[str_v_pd['DowSubId'] != -1.2345]

    return str_v_pd[['SubId','DowSubId','n_up_sub']]


def arcgis_raster_setnull_array(input,output,values,cellSize,SptailRef,work_folder,snapRaster):
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arcpy.env.overwriteOutput = True
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(snapRaster).extent
    arcpy.env.snapRaster =  snapRaster
    lowerLeft = arcpy.Point(arcpy.Describe(snapRaster).extent.XMin,arcpy.Describe(snapRaster).extent.YMin)

    input_arr = arcpy.RasterToNumPyArray(input,nodata_to_value=-9999)
    if (len(values) > 0):
        mask = np.isin(input_arr, values)
        input_arr[mask] = -9999

    input_raster = arcpy.NumPyArrayToRaster(input_arr,lowerLeft,cellSize,
                                     value_to_nodata=-9999)
    input_raster.save(output)
    return input_raster,input_arr


def create_pour_points_with_lakes(str_r,str_v,cat_no_lake,sl_lakes,sl_connected_lake,acc,pourpoints_with_lakes,
                                  lake_inflow_pourpoints,lake_outflow_pourpoints,
                                  catchment_pourpoints_outside_lake,cellSize,SptailRef,work_folder,snapRaster):

    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")
    arcpy.env.overwriteOutput = True
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(snapRaster).extent
    arcpy.env.snapRaster =  snapRaster
    lowerLeft = arcpy.Point(arcpy.Describe(snapRaster).extent.XMin,arcpy.Describe(snapRaster).extent.YMin)


    routing_info  = return_routing_info_using_str_v(str_v,cellSize,SptailRef,work_folder,snapRaster)

    find_zonal_maximum_or_minimum_point_on_the_raster(sl_lakes+"_r",acc,lake_outflow_pourpoints,"MAXIMUM",cellSize,SptailRef,work_folder,snapRaster)
    find_zonal_maximum_or_minimum_point_on_the_raster(str_r,acc,"sub_outlet","MAXIMUM",cellSize,SptailRef,work_folder,snapRaster)
    find_zonal_maximum_or_minimum_point_on_the_raster(str_r,acc,"sub_inlet","MINIMUM",cellSize,SptailRef,work_folder,snapRaster)



    ExtractMultiValuesToPoints("sub_inlet" + "_v", [[sl_lakes+"_r","sl_cl_id"],[str_r,"IL_SubId"],[acc,"MaxAcc_cat"]], "NONE")
    Overlay_str_lake = Combine([str_r,sl_lakes+"_r"])
    Overlay_str_lake.save("unique_lake_str")

    uniq_lake_str_pd = read_table_as_pandas("unique_lake_str",['value', 'str_r', 'selected_lakes_r'],work_folder)
    str_start_pt_lakeid  =  pd.DataFrame.spatial.from_featureclass("sub_inlet" + "_v")
    str_start_pt_lakeid = str_start_pt_lakeid[['IL_SubId', 'sl_cl_id','MaxAcc_cat']]
    str_start_pt_lakeid = str_start_pt_lakeid.fillna(-1)
    str_start_pt_lakeid['SubId'] = str_start_pt_lakeid['IL_SubId']

    routing_info = routing_info.merge(str_start_pt_lakeid, on="SubId", how="left")
    (
        str_id_within_lakes,
        non_lake_inflow_segs,
        str_id_lake_inlfow,
    ) = return_non_lake_inflow_segs_and_segs_within_lakes(
        uniq_lake_str_pd['value'].values, uniq_lake_str_pd['str_r'].values, uniq_lake_str_pd['value'].values, uniq_lake_str_pd['selected_lakes_r'].values, routing_info, str_start_pt_lakeid
    )

    # obtain sub outlet outside lakes
    arcgis_raster_setnull_array("sub_outlet_r",catchment_pourpoints_outside_lake,str_id_within_lakes,cellSize,SptailRef,work_folder,snapRaster)

    # for each unique lake and str value, only keep lake inflow lake_str value
    temp,lake_inflow_unique_lakestr_array = arcgis_raster_setnull_array("unique_lake_str","unique_lake_str_inflow",non_lake_inflow_segs,cellSize,SptailRef,work_folder,snapRaster)

    # find the inflow points within lake domain
    find_zonal_maximum_or_minimum_point_on_the_raster("unique_lake_str_inflow",acc,"lake_inlet_inlake","MINIMUM",cellSize,SptailRef,work_folder,snapRaster)

    # get the lake_Str id of lake inflow strs
    lake_inflow_unique_lakestr_ids = np.unique(lake_inflow_unique_lakestr_array)
    lake_inflow_unique_lakestr_ids = lake_inflow_unique_lakestr_ids[lake_inflow_unique_lakestr_ids!=-9999]
    # get lake inflow str value
    lake_inlet_subids = uniq_lake_str_pd[uniq_lake_str_pd["value"].isin(lake_inflow_unique_lakestr_ids)]["str_r"].values

    # buffer the inflow points within lake domain
    arcpy.analysis.Buffer(in_features = "lake_inlet_inlake_v", out_feature_class = "lake_inlet_inlake_vbf",buffer_distance_or_field = str(2.5*cellSize) + "  Meters",method="PLANAR")
    # obtain str within the buffer

    arcpy.PolygonToRaster_conversion("lake_inlet_inlake_vbf","grid_code", "str_around_lake_inlets",
                                 "MAXIMUM_COMBINED_AREA",None, snapRaster)

    # remove str that are not lake inflow strs
    excludeids = routing_info[~routing_info["SubId"].isin(lake_inlet_subids)]["SubId"].values
    arcgis_raster_setnull_array("str_r","lake_inflow_str",excludeids,cellSize,SptailRef,work_folder,snapRaster)

    remove_value_not_on_lake_inflow_str = Con(IsNull("lake_inflow_str"),-9999,"str_around_lake_inlets")
    remove_value_inside_lakes = Con(IsNull(sl_lakes+"_r"),remove_value_not_on_lake_inflow_str,-9999)
    remove_value_inside_lakes_no999 = SetNull(remove_value_inside_lakes,remove_value_inside_lakes,"VALUE = -9999")
    remove_value_inside_lakes_no999.save("lake_inflow_seg_outside_lake")
    find_zonal_maximum_or_minimum_point_on_the_raster("lake_inflow_seg_outside_lake",acc,lake_inflow_pourpoints,"MAXIMUM",cellSize,SptailRef,work_folder,snapRaster)

    arcpy.RasterToPoint_conversion(catchment_pourpoints_outside_lake, catchment_pourpoints_outside_lake + "_v", "VALUE")

    cat_outlet =  pd.DataFrame.spatial.from_featureclass(catchment_pourpoints_outside_lake + "_v")
    lake_inlets =  pd.DataFrame.spatial.from_featureclass(lake_inflow_pourpoints + "_v")
    lake_outlets =  pd.DataFrame.spatial.from_featureclass(lake_outflow_pourpoints + "_v")

    pourpoints_withlake_var = pd.concat([cat_outlet, lake_inlets,lake_outlets], ignore_index=True)
    pourpoints_withlake_var = pourpoints_withlake_var.reset_index()
    pourpoints_withlake_var["SubId"] = pourpoints_withlake_var.index + 1
    pourpoints_withlake_var = pourpoints_withlake_var[["SubId","SHAPE"]]
    pourpoints_withlake_var.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb",pourpoints_with_lakes + "_v"),overwrite=True,sanitize_columns=False)

    Str_Id = uniq_lake_str_pd['str_r'].values
    CL_Id = uniq_lake_str_pd['selected_lakes_r'].values

    Lakes_WIth_Multi_Outlet, Remove_Str = Check_If_Lake_Have_Multi_OutLet(
        CL_Id, Str_Id, routing_info
    )

    return Lakes_WIth_Multi_Outlet, Remove_Str


def calculate_bkf_width_depth(catinfo):
    idx = catinfo.index
    for i in range(0, len(idx)):
        idx_i = idx[i]
        da = catinfo.loc[idx_i, "DrainArea"] / 1000 / 1000  # m2 to km2
        catid = catinfo.loc[idx_i, "SubId"]
        k = catinfo.loc[idx_i, "k"]
        c =  catinfo.loc[idx_i, "c"]

        if k < 0:
            k_sub =  0.00450718
            c_sub =  0.98579699
        else:
            k_sub = k
            c_sub = c
            if k_sub < 0:
                k_sub =  0.00450718
                c_sub =  0.98579699

        q = func_Q_DA(da, k_sub, c_sub)
        catinfo.loc[idx_i, "BkfWidth"] =  max(7.2 * q ** 0.5,min_bkf_width)
        catinfo.loc[idx_i, "BkfDepth"] =  max(0.27 * q ** 0.3,min_bkf_depth)
        catinfo.loc[idx_i, "Q_Mean"] = q

        if catinfo.loc[idx_i, "Lake_Cat"] < 2:
            if catinfo.loc[idx_i, "Max_DEM"] < 0:
                catinfo.loc[idx_i, "Max_DEM"] = catinfo.loc[idx_i, "MeanElev"]
                catinfo.loc[idx_i, "Min_DEM"] = catinfo.loc[idx_i, "MeanElev"]

    catinfo_riv = catinfo.loc[(catinfo["Lake_Cat"] < 2) & (catinfo["RivLength"] != -1.2345)].copy(deep=True)

    Seg_IDS = catinfo_riv["Seg_ID"].values
    Seg_IDS = np.unique(Seg_IDS)

    for iseg in range(0, len(Seg_IDS)):
        i_seg_id = Seg_IDS[iseg]
        i_seg_info = catinfo_riv[catinfo_riv["Seg_ID"] == i_seg_id]
        if len(i_seg_info["RivLength"].values[i_seg_info["RivLength"].values > 0]) > 0:
            length_seg = np.sum(i_seg_info["RivLength"].values[i_seg_info["RivLength"].values > 0])
        else:
            length_seg = 1
        qmean_seg = np.average(i_seg_info["Q_Mean"].values[i_seg_info["Q_Mean"].values > 0])
        width_seg = np.average(i_seg_info["BkfWidth"].values[i_seg_info["BkfWidth"].values > 0])
        depth_Seg = np.average(i_seg_info["BkfDepth"].values[i_seg_info["BkfDepth"].values > 0])
        floodn_Seg = np.average(i_seg_info["FloodP_n"].values[i_seg_info["FloodP_n"].values > 0])
        basslp_Seg = np.average(i_seg_info["BasSlope"].values[i_seg_info["BasSlope"].values > 0])
        basasp_Seg = np.average(i_seg_info["BasAspect"].values[i_seg_info["BasAspect"].values > 0])
        baselv_Seg = np.average(i_seg_info["MeanElev"].values[i_seg_info["MeanElev"].values > 0])

        if len(i_seg_info["Max_DEM"].values[i_seg_info["Max_DEM"].values > -999999999999999999]) > 0:
            max_elve_seg = np.max(i_seg_info["Max_DEM"].values[i_seg_info["Max_DEM"].values > -9999999999999999])
            min_elve_seg = np.min(i_seg_info["Min_DEM"].values[i_seg_info["Min_DEM"].values > -9999999999999999])
        else:
            max_elve_seg = baselv_Seg
            min_elve_seg = baselv_Seg

        slope_seg = (max_elve_seg - min_elve_seg) / length_seg
        if slope_seg < 0.000000001:
            slope_seg = min_riv_slope  #### Needs to update later

        n_seg = calculateChannaln(width_seg, depth_Seg, qmean_seg, slope_seg)

        for i in range(0, len(i_seg_info)):
            subid = i_seg_info["SubId"].values[i]
            max_elve_rch = i_seg_info["Max_DEM"].values[i]
            min_elve_rch = i_seg_info["Min_DEM"].values[i]
            length_rch = max(i_seg_info["RivLength"].values[i],min_riv_lenth)
            qmean_rch = i_seg_info["Q_Mean"].values[i]
            width_rch = i_seg_info["BkfWidth"].values[i]
            depth_rch = i_seg_info["BkfDepth"].values[i]
            floodn_rch = i_seg_info["FloodP_n"].values[i]

            if min_elve_rch < -2000:
                if i_seg_info["MeanElev"].values[i] > -2000:
                    min_elve_rch = i_seg_info["MeanElev"].values[i]
                else:
                    min_elve_rch = baselv_Seg

            if max_elve_rch < -2000:
                if i_seg_info["MeanElev"].values[i] > -2000:
                    max_elve_rch = i_seg_info["MeanElev"].values[i]
                else:
                    max_elve_rch = baselv_Seg

            slope_rch = i_seg_info["RivSlope"].values[i]

            if slope_rch < min_riv_slope or slope_rch > max_riv_slope:
                if slope_seg >= min_riv_slope and slope_seg <= max_riv_slope:
                    slope_rch = slope_seg

            slope_rch = max(slope_rch,min_riv_slope)
            slope_rch = min(slope_rch,max_riv_slope)

            n_rch = calculateChannaln(width_rch, depth_rch, qmean_rch, slope_rch)

            if n_rch < min_manning_n or n_rch > max_manning_n:
                if n_seg >= min_manning_n and n_seg <= max_manning_n:
                    n_rch = n_seg

            n_rch = max(n_rch,min_manning_n)
            n_rch = min(n_rch,max_manning_n)

            catinfo.loc[catinfo["SubId"] == subid, "RivSlope"] = slope_rch
            catinfo.loc[catinfo["SubId"] == subid, "Ch_n"] = n_rch
            catinfo.loc[catinfo["SubId"] == subid, "RivLength"] = length_rch


            # if subid == 504:
            #     print(floodn_rch)


            if floodn_rch < 0:
                if floodn_Seg > 0:
                    floodn_rch = floodn_Seg
                else:
                    floodn_rch = DEFALUT_FLOOD_N

            floodn_rch = max(floodn_rch,n_rch)
            catinfo.loc[catinfo["SubId"] == subid, "FloodP_n"] = floodn_rch
            catinfo.loc[catinfo["SubId"] == subid, "Max_DEM"] = max_elve_rch
            catinfo.loc[catinfo["SubId"] == subid, "Min_DEM"] = min_elve_rch

            # if subid == 504:
            #     print(floodn_rch)
            #     asdf

            if i_seg_info["BasSlope"].values[i] <= 0:
                if basslp_Seg > 0:
                    catinfo.loc[catinfo["SubId"] == subid, "BasSlope"] = basslp_Seg
                else:
                    catinfo.loc[catinfo["SubId"] == subid, "BasSlope"] = 0

            if i_seg_info["BasAspect"].values[i] <= 0:
                if basasp_Seg > 0:
                    catinfo.loc[catinfo["SubId"] == subid, "BasAspect"] = basasp_Seg
                else:
                    catinfo.loc[catinfo["SubId"] == subid, "BasAspect"] = 0

            if i_seg_info["MeanElev"].values[i] < 0:
                if baselv_Seg > 0:
                    catinfo.loc[catinfo["SubId"] == subid, "MeanElev"] = baselv_Seg
                else:
                    catelv = catinfo['MeanElev'].values
                    catelv = catelv[catelv > 0]
                    catinfo.loc[catinfo["SubId"] == subid, "MeanElev"] =np.average(catelv)

            # print(subid,n_rch,n_seg,min_manning_n,max_manning_n)

    mask1 = catinfo["DA_Chn_Slp"] < min_riv_slope
    mask2 = catinfo["DA_Chn_Slp"] > max_riv_slope
    mask = np.logical_or(mask1,mask2)
    catinfo.loc[mask,"DA_Chn_Slp"] = catinfo.loc[mask,"RivSlope"]

    return catinfo
