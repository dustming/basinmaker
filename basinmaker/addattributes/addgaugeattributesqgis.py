from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
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
    obsname = "obs"  # input_geo_names["obsname"]

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
        raster=obsname,
        column="obsid_pour",
    )

    grass.run_command(
        "v.out.ogr",
        input=outlet_pt_info,
        output=os.path.join(grassdb, outlet_pt_info + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )

    ### read catchment
    sqlstat = "SELECT SubId,obsid_pour FROM %s" % (outlet_pt_info,)
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-9999)
    outletinfo = outletinfo.loc[outletinfo["obsid_pour"] > 0]

    sqlstat = "SELECT  %s,%s,%s,%s,%s FROM %s" % (
        obs_attributes[0],
        obs_attributes[1],
        obs_attributes[2],
        obs_attributes[3],
        obs_attributes[0] + "n",
        snapped_obs_points,
    )
    gaugeinfo = pd.read_sql_query(sqlstat, con)
    gaugeinfo = gaugeinfo.fillna(-9999)

    for i in range(0, len(outletinfo)):
        catid = outletinfo["SubId"].values[i]
        obsid = outletinfo["obsid_pour"].values[i]
        catrow = catinfo["SubId"] == catid
        subgaugeinfo = gaugeinfo.loc[gaugeinfo[obs_attributes[0] + "n"] == obsid]

        catinfo.loc[catrow, "Has_POI"] = obsid
        if obsid > 79000:
            continue

        if len(subgaugeinfo) > 0:
            catinfo.loc[catrow, "Has_POI"] = subgaugeinfo[obs_attributes[0]].values[0]
            catinfo.loc[catrow, "DA_Obs"] = subgaugeinfo[obs_attributes[2]].values[0]
            catinfo.loc[catrow, "Obs_NM"] = subgaugeinfo[obs_attributes[1]].values[0]
            catinfo.loc[catrow, "SRC_obs"] = subgaugeinfo[obs_attributes[3]].values[0]

    PERMANENT.close()
    return catinfo
