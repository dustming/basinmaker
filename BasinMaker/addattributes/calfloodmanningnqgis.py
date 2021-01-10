from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3
import pandas as pd


def calculate_flood_plain_manning_n(
    grassdb,
    grass_location,
    qgis_prefix_path,
    catinfo,
    path_landuse="#",
    path_landuse_info="#",
    riv_seg="#",
):

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

    # read lanning's n and landuse type table
    landuse_and_n_table = pd.read_csv(path_landuse_info, sep=",")

    write_grass_reclass_rule_from_table(
        landuse_and_n_table.values, os.path.join(grassdb, "landuse_manning_rules.csv")
    )

    # viturally  landuse dataset
    grass_raster_r_external(grass, input=path_landuse, output="landuse_in")
    # clip raster with mask in grass env
    grass_raster_r_clip(grass, input="landuse_in", output="landuse")
    # reclass landuse to manning's coefficient value *1000
    grass_raster_r_reclass(
        grass,
        input="landuse",
        output="landuse_Manning1",
        rules=os.path.join(grassdb, "landuse_manning_rules.csv"),
    )
    # calcuate real manning's coefficient for each landuse grid
    grass_raster_r_mapcalc(
        grass, expression="landuse_Manning = float(landuse_Manning1)/1000"
    )

    ### add averaged manning coefficent along the river network into river attribut table
    grass.run_command(
        "v.rast.stats",
        map="nstr_nfinalcat_F",
        raster="landuse_Manning",
        column_prefix="mn",
        method=["average"],
    )

    ### read length and maximum and minimum dem along channel
    sqlstat = "SELECT Gridcode,mn_average FROM nstr_nfinalcat_F"
    rivleninfo = pd.read_sql_query(sqlstat, con)
    rivleninfo = rivleninfo.fillna(-9999)

    for i in range(0, len(rivleninfo)):
        catid = rivleninfo["Gridcode"].values[i]
        catrow = catinfo["SubId"] == catid
        catinfo.loc[catrow, "FloodP_n"] = rivleninfo["mn_average"].values[i]

    PERMANENT.close()
    return catinfo
