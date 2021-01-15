from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3


def select_lakes_by_area_r(
    grass,
    con,
    str_r,
    lakes="alllake",
    lakes_lg_cl_thres = 'lakes_lg_cl_thres',
    lakes_lg_ncl_thres = 'lakes_lg_ncl_thres',
    connected_lake="connect_lake",
    non_connected_lake="nonconnect_lake",
    str_connected_lake="str_connected_lake",
    sl_connected_lake="sl_connect_lake",
    sl_non_connected_lake="sl_nonconnect_lake",
    sl_lakes="selected_lakes",
    sl_str_connected_lake="str_sl_connected_lake",
):

    exp = "%s = if(isnull('%s'),null(),int(%s))" % (
        sl_connected_lake,
        lakes_lg_cl_thres,
        connected_lake,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    exp = "%s = if(isnull('%s'),null(),int(%s))" % (
        sl_non_connected_lake,
        lakes_lg_ncl_thres,
        non_connected_lake,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
        
    exp = "%s = if(isnull('%s'),int(%s),int(%s))" % (
        sl_lakes,
        sl_connected_lake,
        sl_non_connected_lake,
        sl_connected_lake,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    exp = "%s = if(isnull('%s'),null(),%s)" % (sl_str_connected_lake, str_r, sl_connected_lake)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)
    
    return
