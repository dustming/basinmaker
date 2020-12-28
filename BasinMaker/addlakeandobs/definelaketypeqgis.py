from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from utilities import *
import sqlite3


def define_connected_and_non_connected_lake_type(
    grass,
    con,
    garray,
    str_r="str_grass_r",
    lake="alllake",
    connected_lake="connect_lake",
    non_connected_lake="nonconnect_lake",
    str_connected_lake="str_connected_lake",
):

    #### overlay str and lakes
    exp = "%s = if(isnull('%s'),null(),%s)" % (str_connected_lake, str_r, lake)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    ### obtain connected lake ids
    Connect_Lake_Ids, temp = generate_stats_list_from_grass_raster(
        grass, mode=1, input_a=str_connected_lake
    )
    #### create non connected lake raster
    grass.run_command("g.copy", rast=(lake, non_connected_lake), overwrite=True)
    grass.run_command(
        "r.null", map=non_connected_lake, setnull=Connect_Lake_Ids, overwrite=True
    )
    #### create potential connected lake raster
    exp = "%s = if(isnull('%s'),%s,null())" % (connected_lake, non_connected_lake, lake)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    return


def generate_stats_list_from_grass_raster(
    grass, mode=1, input_a="P_Connect_Lake", input_b="str_grass_r", method="max"
):

    list_a = []
    list_b = []
    if mode == 1:
        p = grass.pipe_command("r.category", map=input_a, separator=" ", quiet=True)
    elif mode == 2:
        p = grass.pipe_command(
            "r.stats", input=[input_a, input_b], flags="n", quiet=True
        )
    else:
        print("error mode opton did not exist")

    for line in p.stdout:
        line_str = line.decode("utf8").rstrip("\r\n").split(" ")
        list_a.append(int(line_str[0]))
        if mode != 1:
            list_b.append(int(line_str[1]))
    p.wait()

    return list_a, list_b
