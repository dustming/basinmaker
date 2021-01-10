import pandas as pd
import numpy as np
from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import os


def join_pandas_table_to_vector_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    vector_name,
    pd_table,
    column_types,
    columns_names,
):

    str = ",".join(column_types)
    WriteStringToFile(str, os.path.join(grassdb, "catinfo_riv.csvt"), "w")

    pd_table.to_csv(os.path.join(grassdb, "catinfo_riv.csv"), index=None, header=True)

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

    ### add catchment info to all river segment
    grass.run_command(
        "db.in.ogr",
        input=os.path.join(grassdb, "catinfo_riv.csv"),
        output="result_riv",
        overwrite=True,
    )
    grass.run_command(
        "v.db.join",
        map=vector_name,
        column="Gridcode",
        other_table="result_riv",
        other_column="SubId",
        overwrite=True,
    )
    grass.run_command(
        "v.out.ogr",
        input=vector_name,
        output=os.path.join(grassdb, vector_name + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )

    grass.run_command(
        "v.out.ogr",
        input="Final_OL_v",
        output=os.path.join(grassdb, "Final_OL_v.shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )
