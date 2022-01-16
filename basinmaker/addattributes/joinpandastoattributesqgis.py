import pandas as pd
import numpy as np
from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import os
import sqlite3


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
    pd_table.loc[pd_table['Has_POI'] <= 0,'Has_POI'] = 0
    pd_table.loc[pd_table['Has_POI'] > 0,'Has_POI'] = 1
    pd_table.loc[pd_table['Lake_Cat'] <= 0,'Lake_Cat'] = 0
    pd_table.loc[pd_table['LakeArea'] <= 0,'LakeArea'] = 0
    pd_table.loc[pd_table['HyLakeId'] <= 0,'HyLakeId'] = 0
    pd_table.loc[pd_table['Laketype'] <= 0,'Laketype'] = 0
    pd_table.loc[pd_table['LakeVol'] <= 0,'LakeVol'] = 0
    pd_table.loc[pd_table['LakeDepth'] <= 0,'LakeDepth'] = 0
    pd_table.loc[pd_table['LakeArea'] <= 0,'LakeArea'] = 0
    pd_table.loc[pd_table['DA_Obs'] <= 0,'DA_Obs'] = 0
    pd_table.loc[pd_table['DA_error'] <= 0,'DA_error'] = 0
    
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

    con = sqlite3.connect(
        os.path.join(grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db")
    )
    
#    pd_table.to_sql('catinfo_pd',con,if_exists='replace',index_label='ix')
    
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
        "v.db.update",
        map=vector_name,
        column="DA_error",
        qcol="DrainArea/1000/1000/DA_Obs",
        where="Has_POI > 0",
    )

    grass.run_command(
        "v.out.ogr",
        input=vector_name,
        output=os.path.join(grassdb, vector_name + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )
