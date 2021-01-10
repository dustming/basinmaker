from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3


def add_gauge_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    catinfo,
    input_geo_names,
    obs_attributes=[],
):

    outlet_pt_info = input_geo_names["outlet_pt_info"]
    snapped_obs_points = input_geo_names["snapped_obs_points"]

    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    con = sqlite3.connect(
        os.path.join(grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db")
    )

    grass.run_command(
        "v.what.rast",
        map=outlet_pt_info,
        raster=snapped_obs_points,
        column="obsid_pour",
    )

    grass_raster_v_db_join(
        grass,
        map=outlet_pt_info,
        column="obsid_pour",
        other_table=snapped_obs_points,
        other_column=obs_attributes[0] + "n",
    )

    ### read catchment
    sqlstat = "SELECT SubId,%s,%s,%s,%s FROM %s" % (
        obs_attributes[0],
        obs_attributes[1],
        obs_attributes[2],
        obs_attributes[3],
        outlet_pt_info,
    )
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-9999)

    for i in range(0, len(outletinfo)):
        catid = outletinfo["SubId"].values[i]
        catrow = catinfo["SubId"] == catid
        obsid = outletinfo[obs_attributes[0]].values[i]

        if obsid > 0:
            catinfo.loc[catrow, "IsObs"] = obsid
            catinfo.loc[catrow, "DA_Obs"] = outletinfo[obs_attributes[2]].values[i]
            catinfo.loc[catrow, "Obs_NM"] = outletinfo[obs_attributes[1]].values[i]
            catinfo.loc[catrow, "SRC_obs"] = outletinfo[obs_attributes[3]].values[i]

    PERMANENT.close()
    return catinfo
