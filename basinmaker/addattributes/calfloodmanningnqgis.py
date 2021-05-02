from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import sqlite3
import pandas as pd
from basinmaker.preprocessing.preprocessrasterqgis import preprocess_raster


def calculate_flood_plain_manning_n(
    grassdb,
    grass_location,
    qgis_prefix_path,
    catinfo,
    input_geo_names,
    path_landuse="#",
    path_landuse_info="#",
):

    mask = input_geo_names["mask"]

    preprocess_raster(
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        mask=mask,
        raster_path=path_landuse,
        raster_name="landuse",
    )

    cat_riv_info = input_geo_names["cat_riv_info"]
    outlet_pt_info = input_geo_names["outlet_pt_info"]

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

    grass_raster_r_external(
        grass, input=os.path.join(grassdb, "landuse_proj" + ".tif"), output="landuse_in"
    )
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
        map=cat_riv_info,
        raster="landuse_Manning",
        column_prefix="mn",
        method=["average"],
    )

    grass.run_command(
        "v.what.rast", map=outlet_pt_info, raster= "landuse_Manning", column="mn_average"
    )
        
    # grass.run_command(
    #     "v.rast.stats",
    #     map=cat_ply_info,
    #     raster="landuse_Manning",
    #     column_prefix="mn",
    #     method=["average"],
    # )    

    ### read length and maximum and minimum dem along channel
    sqlstat = "SELECT Gridcode,mn_average FROM %s" % (cat_riv_info)
    rivleninfo = pd.read_sql_query(sqlstat, con)
    rivleninfo = rivleninfo.fillna(-9999)
    
    ### read length and maximum and minimum dem along channel
    sqlstat = "SELECT SubId,mn_average FROM %s" % (outlet_pt_info)
    outletpoint = pd.read_sql_query(sqlstat, con)
    outletpoint = outletpoint.fillna(-9999)

    for i in range(0, len(rivleninfo)):
        catid = rivleninfo["Gridcode"].values[i]
        catrow = catinfo["SubId"] == catid
        floodn = rivleninfo["mn_average"].values[i]
        
        if floodn < 0:
            floodn = outletpoint.loc[outletpoint['SubId'] ==catid, "mn_average"].values[0]
                        
        catinfo.loc[catrow, "FloodP_n"] = floodn

    PERMANENT.close()
    return catinfo
