from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3


def define_pour_points_with_lakes(
    grass,
    con,
    garray,
    remove_lake_inlets = False,
    str_r="str_grass_r",
    cat_no_lake="cat_no_lake",
    sl_lakes="sl_lakes",
    sl_connected_lake="sl_connected_lake",
    sl_str_connected_lake="sl_str_connected_lake",
    acc="acc",
    pourpoints_with_lakes="pourpoints_with_lakes",
    lake_inflow_pourpoints="lake_inflow_pourpoints",
    lake_outflow_pourpoints="lake_outflow_pourpoints",
    catchment_pourpoints_outside_lake="catchment_pourpoints_outside_lake",
):

    # define catchment pourpoints and routing info
    routing_info = generate_routing_info_of_catchments(
        grass, con, cat=cat_no_lake, acc=acc, str=str_r, Name="cat1"
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
    exp = "'%s' =if(float(%s) == float(%s),%s,null())" % (
        lake_outflow_pourpoints,
        acc,
        "sl_lakes_maxacc",
        sl_lakes,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    # add lake inlets

    #### overlay str and lakes
    exp = "%s = if(isnull(int(%s)),null(),%s)" % (
        "connect_lake_str",
        sl_str_connected_lake,
        str_r,
    )
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
    grass.run_command(
        "v.what.rast", map="cat1_IL_v", raster=sl_connected_lake, column="sl_cl_id"
    )
    sqlstat = "SELECT IL_SubId, sl_cl_id FROM %s" % ("cat1_IL_v")
    str_start_pt_lakeid = pd.read_sql_query(sqlstat, con)
    str_start_pt_lakeid = str_start_pt_lakeid.fillna(-1)

    (
        str_id_within_lakes,
        non_lake_inflow_segs,
        str_id_lake_inlfow,
    ) = return_non_lake_inflow_segs_and_segs_within_lakes(
        riv_lake_id, str_id, riv_lake_id2, cl_id, routing_info, str_start_pt_lakeid
    )

    # remove cat outlet that within the lake
    grass.run_command(
        "g.copy", rast=("cat1_OL", catchment_pourpoints_outside_lake), overwrite=True
    )
    grass.run_command(
        "r.null",
        map=catchment_pourpoints_outside_lake,
        setnull=str_id_within_lakes,
        overwrite=True,
    )

    # remove non lake inflow river segment
    grass.run_command(
        "g.copy", rast=("unique_lake_str", "unique_lake_str_inflow"), overwrite=True
    )
    grass.run_command(
        "r.null",
        map="unique_lake_str_inflow",
        setnull=non_lake_inflow_segs,
        overwrite=True,
    )

    ## find lake inflow points from inflow segments

    # get all str that are lake inflow str
    all_river_ids, temp = generate_stats_list_from_grass_raster(
        grass, mode=1, input_a=str_r
    )
    str_id_non_lake_inlfow = [x for x in all_river_ids if x not in str_id_lake_inlfow]
    grass.run_command("g.copy", rast=(str_r, "lake_inflow_str"), overwrite=True)
    grass.run_command(
        "r.null", map="lake_inflow_str", setnull=str_id_non_lake_inlfow, overwrite=True
    )

    # find the lake inlet point within the lake

    grass.run_command(
        "r.stats.zonal",
        base="unique_lake_str_inflow",
        cover="cat1_acc_riv",
        method="min",
        output="lake_inflow_seg_minacc",
        overwrite=True,
    )
    ### Find the grid that equal to the max acc, thus this is the outlet grids
    exp = "%s =if(float(%s) == float(%s),%s,null())" % (
        "lake_inflow_IL",
        "cat1_acc_riv",
        "lake_inflow_seg_minacc",
        "unique_lake_str_inflow",
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    # extent the inflow seg to lake outside
    grass.run_command(
        "r.grow",
        input="lake_inflow_IL",
        output="lake_inflow_IL_grow",
        radius=1.5,
        overwrite=True,
    )

    # set non river cells to non
    exp = "%s = if(isnull('%s'),null(),int(%s))" % (
        "extented_lake_inflow_seg",
        "lake_inflow_str",
        "lake_inflow_IL_grow",
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    grass.run_command(
        "r.stats.zonal",
        base="extented_lake_inflow_seg",
        cover="cat1_acc_riv",
        method="min",
        output="extented_lake_inflow_seg_minacc",
        overwrite=True,
    )
    ### Find the grid that equal to the max acc, thus this is the outlet grids
    exp = "%s =if(float(%s) == float(%s),%s,null())" % (
        lake_inflow_pourpoints,
        "cat1_acc_riv",
        "extented_lake_inflow_seg_minacc",
        "extented_lake_inflow_seg",
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    if remove_lake_inlets:
        exp = "%s =if(float(%s) > -1000000,null(),null())" % (
            lake_inflow_pourpoints+'t',
            lake_inflow_pourpoints,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

        exp = "%s = int(%s)" % (
            lake_inflow_pourpoints,
            lake_inflow_pourpoints+'t',
        )

        grass.run_command("r.mapcalc", expression=exp, overwrite=True)


    #### create a unique id for overlaied lake and river
    grass.run_command(
        "r.cross",
        input=[
            lake_inflow_pourpoints,
            catchment_pourpoints_outside_lake,
            lake_outflow_pourpoints,
        ],
        output=pourpoints_with_lakes,
        overwrite=True,
    )

    exp = "%s = %s + 1" % (pourpoints_with_lakes, pourpoints_with_lakes)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    grass.run_command(
        "r.to.vect",
        input=pourpoints_with_lakes,
        output=pourpoints_with_lakes,
        type="point",
        overwrite=True,
        flags="v",
    )

    grass.run_command(
        "r.to.vect",
        input=lake_inflow_pourpoints,
        output=lake_inflow_pourpoints,
        type="point",
        overwrite=True,
        flags="v",
    )
    grass.run_command(
        "r.to.vect",
        input=catchment_pourpoints_outside_lake,
        output=catchment_pourpoints_outside_lake,
        type="point",
        overwrite=True,
        flags="v",
    )
    grass.run_command(
        "r.to.vect",
        input=lake_outflow_pourpoints,
        output=lake_outflow_pourpoints,
        type="point",
        overwrite=True,
        flags="v",
    )

    #
    ##### obtain lake id and coresponding overlaied str id
    CL_Id, Str_Id = generate_stats_list_from_grass_raster(
        grass, mode=2, input_a=sl_connected_lake, input_b=str_r
    )

    Lakes_WIth_Multi_Outlet, Remove_Str = Check_If_Lake_Have_Multi_OutLet(
        CL_Id, Str_Id, routing_info
    )

    return Lakes_WIth_Multi_Outlet, Remove_Str
