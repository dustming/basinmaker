from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from processing_functions_attribute_table import *
from utilities import *
import sqlite3
from addlakeandobs.defineroutinginfoqgis import generate_routing_info_of_catchments
from addlakeandobs.definelaketypeqgis import generate_stats_list_from_grass_raster

def define_pour_points_with_lakes(
    grass,
    con,
    garray,
    str_r="str_grass_r",
    cat_no_lake = "cat_no_lake",
    sl_lakes = "sl_lakes",
    sl_connected_lake = "sl_connected_lake", 
    sl_str_connected_lake = "sl_str_connected_lake",
    acc ="acc",
    lake_pourpoints = "lake_pourpoints",
):
    
    # define catchment pourpoints and routing info 
    routing_info = generate_routing_info_of_catchments(
        grass, con, cat=cat_no_lake, acc=acc, str=str_r,Name='cat1'
    )  
    routing_info = routing_info.fillna(-1)
    
    # define lake outlets 
    grass.run_command(
        "r.stats.zonal",
        base=sl_lakes,
        cover=acc,
        method="max",
        output="sl_lakes_maxacc",
        overwrite=True,
    )
    exp = "'%s' =if(%s == int('%s'),%s,null())" % (
        "pourpoints_sl_lakes",
        acc,
        "sl_lakes_maxacc",
        sl_lakes,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)      

    # add lake inlets 

    #### overlay str and lakes
    exp = "%s = if(isnull('%s'),null(),%s)" % ("connect_lake_str",sl_str_connected_lake, str_r)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    
    #### create a unique id for overlaied lake and river
    grass.run_command(
        "r.cross",
        input=["connect_lake_str", sl_str_connected_lake],
        output="unique_lake_str",
        flags="z",
        overwrite=True,
    )      

    ##### obtain lake id and coresponding overlaied str id
    riv_lake_id, str_id = generate_stats_list_from_grass_raster(
        grass, mode=2, input_a="unique_lake_str", input_b="connect_lake_str"
    )
    
    riv_lake_id2, cl_id = generate_stats_list_from_grass_raster(
        grass, mode=2, input_a="unique_lake_str", input_b=sl_str_connected_lake
    )
    
    ### get the lake raster value at the begining of each river segments,if the river 
    ### the beginning of the river segment located within the lake, than it can not 
    ### be used as the lake inflow segment.
    grass.run_command("v.what.rast", map="cat1_IL_v", raster=sl_connected_lake, column="sl_cl_id")
    sqlstat = "SELECT IL_SubId, sl_cl_id FROM %s" % ("cat1_IL_v")
    str_start_pt_lakeid = pd.read_sql_query(sqlstat, con)
    str_start_pt_lakeid = str_start_pt_lakeid.fillna(-1)

    str_id_within_lakes,non_lake_inflow_segs = return_non_lake_inflow_segs_and_segs_within_lakes(riv_lake_id,str_id,riv_lake_id2,cl_id,routing_info,str_start_pt_lakeid)

    # remove cat outlet that within the lake 
    grass.run_command("g.copy", rast=("cat1_OL", "cat1_OL_outlake"), overwrite=True)
    grass.run_command(
        "r.null", map="cat1_OL_outlake", setnull=str_id_within_lakes, overwrite=True
    )
    
    # remove non lake inflow river segment 
    grass.run_command("g.copy", rast=("unique_lake_str", "unique_lake_str_inflow"), overwrite=True)
    grass.run_command(
        "r.null", map="unique_lake_str_inflow", setnull=non_lake_inflow_segs, overwrite=True
    )    
    
    ## find lake inflow points from inflow segments 
    
    # extent the inflow seg to lake outside 
    # grass.run_command(
    #     "r.grow",
    #     input="unique_lake_str_inflow",
    #     output="unique_lake_str_inflow_grow",
    #     radius=1.5,
    #     overwrite=True,
    # )
    # 
    # # set non river cells to non 
    # exp = "%s = if(isnull('%s'),null(),%s)" % ("extented_lake_inflow_seg","cat1_acc_riv", "unique_lake_str_inflow_grow")
    # grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    # 
    # 
    
    
        
    
    