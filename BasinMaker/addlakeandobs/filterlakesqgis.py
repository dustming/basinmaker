from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from utilities import *
import sqlite3


def select_lakes_by_area_r(
    grass,
    con,
    lake_v_path,
    threshold_con_lake,
    threshold_non_con_lake,
    lakes="alllake",
    connected_lake="connect_lake",
    non_connected_lake="nonconnect_lake",
    str_connected_lake="str_connected_lake",
    sl_connected_lake="sl_connect_lake",
    sl_non_connected_lake="sl_nonconnect_lake",
    sl_lakes="selected_lakes",
    sl_str_connected_lake="str_sl_connected_lake",
):

    alllakinfo = Dbf_To_Dataframe(lake_v_path)

    Lakeid_lt_CL_Remove = alllakinfo.loc[alllakinfo["Lake_area"] < threshold_con_lake][
        "Hylak_id"
    ].values
    Lakeid_lt_NCL_Remove = alllakinfo.loc[
        alllakinfo["Lake_area"] < threshold_non_con_lake
    ]["Hylak_id"].values

    grass_raster_setnull(
        grass,
        connected_lake,
        Lakeid_lt_CL_Remove.tolist(),
        True,
        sl_connected_lake,
    )
    grass_raster_setnull(
        grass,
        non_connected_lake,
        Lakeid_lt_NCL_Remove.tolist(),
        True,
        sl_non_connected_lake,
    )

    exp = "%s = if(isnull('%s'),int(%s),int(%s))" % (
        sl_lakes,
        sl_connected_lake,
        sl_non_connected_lake,
        sl_connected_lake,
    )
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    grass_raster_setnull(
        grass,
        str_connected_lake,
        Lakeid_lt_CL_Remove.tolist(),
        True,
        sl_str_connected_lake,
    )

    return
