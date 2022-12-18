from basinmaker.func.arcgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
from basinmaker.addlakeandobs.modifyfdr import modify_lakes_flow_direction


def add_obs_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_obsfile_in,
    obs_attributes=[],
    search_radius=100,
    path_sub_reg_outlets_v="#",
    max_memroy=1024 * 4,
    pourpoints_add_obs="pourpoints_add_obs",
    snapped_obs_points="snapped_obs_points",
):

    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]
    pourpoints_with_lakes = input_geo_names["pourpoints_with_lakes"]
    lake_outflow_pourpoints = input_geo_names["lake_outflow_pourpoints"]
    cat_add_lake = input_geo_names["cat_add_lake"]

    # define internal file names
    obsname = Internal_Constant_Names["obs"]

    work_folder = grassdb
    # read and define arcgis work enviroments
    if not os.path.exists(work_folder):
        os.makedirs(work_folder)
    arcpy.env.workspace = os.path.join(work_folder,"arcgis.gdb")

    arcpy.env.overwriteOutput = True
    arcpy.CheckOutExtension("Spatial")
    cellSize = float(arcpy.GetRasterProperties_management(dem, "CELLSIZEX").getOutput(0))
    SptailRef = arcpy.Describe(dem).spatialReference
    arcpy.env.XYTolerance = cellSize
    arcpy.arcpy.env.cellSize = cellSize
    arcpy.env.outputCoordinateSystem = arcpy.SpatialReference(int(SptailRef.factoryCode)) ### WGS84
    arcpy.env.extent = arcpy.Describe(dem).extent
    arcpy.env.snapRaster =  dem


    arcpy.Project_management(
        path_obsfile_in,
        "obs_proj",
        arcpy.SpatialReference(int(SptailRef.factoryCode)),
        )
    arcpy.arcpy.analysis.PairwiseClip("obs_proj", mask + '_ply', "obs_clip", "")

    obs_snap = pd.DataFrame.spatial.from_featureclass("obs_clip")
    obs_snap["tSubId"] =  1 + obs_snap.index
    obs_snap.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","obs_vt"),overwrite=True,sanitize_columns=False)

    outSnapPour = SnapPourPoint("obs_vt", acc, 0 * cellSize,obs_attributes[0])
    obs_nolake = Con(IsNull("selected_lakes_r"),outSnapPour)
    obs_nolake.save("obs_r")
    arcpy.RasterToPoint_conversion("obs_r", "obs_v", "VALUE")

    arcpy.management.CalculateGeometryAttributes("str_v_rout", [["Length", "LENGTH"]], "METERS")

    str_v_pd = pd.DataFrame.spatial.from_featureclass("str_v_rout")
    mask1 = str_v_pd["n_up_sub"] == 0
    mask2 = str_v_pd["Length"] > 30*cellSize
    mask = np.logical_and(mask1,mask2)
    headwater_str = str_v_pd.loc[mask]
    headwater_str.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","head_str"),overwrite=True,sanitize_columns=False)

    arcpy.FeatureVerticesToPoints_management("head_str", "str_nodes", "BOTH_ENDS")
    ExtractValuesToPoints("str_nodes",acc,"str_nodes_Acc","INTERPOLATE","VALUE_ONLY")
    head_sub_nodes = pd.DataFrame.spatial.from_featureclass("str_nodes_Acc")
    head_sub_nodes = head_sub_nodes.sort_values(by='RASTERVALU', ascending=False)
    head_sub_nodes = head_sub_nodes.drop_duplicates(subset=['grid_code'], keep='last')
    head_sub_nodes["tSubId"] = head_sub_nodes.index + 1
    head_sub_nodes.spatial.to_featureclass(location=os.path.join(work_folder,"arcgis.gdb","head_node"),overwrite=True,sanitize_columns=False)

    headsnap = SnapPourPoint("head_node", acc, 5 * cellSize,"tSubId")
    head_nolake = Con(IsNull("selected_lakes_r"),headsnap)
    head_nolake.save("head_node_r")
    arcpy.RasterToPoint_conversion(head_nolake, "head_node_v", "VALUE")

    return
