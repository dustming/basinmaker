from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
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
    exp = "%s = if(isnull(int(%s)),null(),%s)" % (str_connected_lake, str_r, lake)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    ### obtain connected lake ids
    Connect_Lake_Ids, temp = generate_stats_list_from_grass_raster(
        grass, mode=1, input_a=str_connected_lake
    )
    #### create non connected lake raster
#    grass.run_command("g.copy", rast=(lake, non_connected_lake), overwrite=True)


    non_connected_lake_array = garray.array(mapname=lake)
    if (len(Connect_Lake_Ids) > 0):
        mask = np.isin(non_connected_lake_array, Connect_Lake_Ids)
        non_connected_lake_array[mask] = -9999

    temparray = garray.array()
    temparray[:, :] = non_connected_lake_array[:, :]
    temparray.write(mapname=non_connected_lake, overwrite=True)
    grass.run_command("r.null", map=non_connected_lake, setnull=[-9999, 0])
    exp = "%s = int(%s)" % (
        non_connected_lake,
        non_connected_lake,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    # grass.run_command(
    #     "r.null", map=non_connected_lake, setnull=Connect_Lake_Ids, overwrite=True
    # )
    #### create potential connected lake raster
    exp = "%s = if(isnull(int(%s)),%s,null())" % (connected_lake, non_connected_lake, lake)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    return
