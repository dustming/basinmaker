from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from utilities import *
import sqlite3
from addlakeandobs.definelaketypeqgis import generate_stats_list_from_grass_raster
from addlakeandobs.defineroutinginfoqgis import generate_routing_info_of_catchments


def add_lake_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    sl_connected_lake,
    sl_non_connected_lake,
    catchments,
    path_lake_ply,
    catinfo,
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

    grass.run_command(
        "r.stats.zonal",
        base=catchments,
        cover=sl_connected_lake,
        method="max",
        output="Connect_Lake_Cat_w_Lake_ID",
        overwrite=True,
    )

    grass.run_command(
        "v.what.rast",
        map="Final_OL_v",
        raster=sl_non_connected_lake,
        column="ncl",
    )

    grass.run_command(
        "v.what.rast",
        map="Final_OL_v",
        raster="Connect_Lake_Cat_w_Lake_ID",
        column="cl",
    )

    ### read catchment
    sqlstat = "SELECT SubId,ncl,cl FROM Final_OL_v"
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-9999)
    lakeinfo = Dbf_To_Dataframe(path_lake_ply)

    for i in range(0, len(outletinfo)):
        catid = outletinfo["SubId"].values[i]

        CL_LakeId_Outlet = outletinfo["cl"].values[i]
        ### catchment can be an connect lake catchment only when outlet of this catchment is
        ### in the lake catchment
        if CL_LakeId_Outlet > 0:
            CL_LakeId = CL_LakeId_Outlet
        else:
            CL_LakeId = -9999

        NCL_LakeId = outletinfo["ncl"].values[i]
        #        print(CL_LakeId,NCL_LakeId)
        ### add lake info
        if CL_LakeId > 0 and NCL_LakeId < 0:
            catinfo.loc[i, "IsLake"] = 1
            slakeinfo = lakeinfo.loc[lakeinfo["Hylak_id"] == CL_LakeId]
            catinfo.loc[i, "HyLakeId"] = CL_LakeId
            catinfo.loc[i, "LakeVol"] = slakeinfo.iloc[0]["Vol_total"]
            catinfo.loc[i, "LakeArea"] = slakeinfo.iloc[0]["Lake_area"]
            catinfo.loc[i, "LakeDepth"] = slakeinfo.iloc[0]["Depth_avg"]
            catinfo.loc[i, "Laketype"] = slakeinfo.iloc[0]["Lake_type"]
        if NCL_LakeId > 0 and CL_LakeId < 0:
            catinfo.loc[i, "IsLake"] = 2
            slakeinfo = lakeinfo.loc[lakeinfo["Hylak_id"] == NCL_LakeId]
            catinfo.loc[i, "HyLakeId"] = NCL_LakeId
            catinfo.loc[i, "LakeVol"] = slakeinfo.iloc[0]["Vol_total"]
            catinfo.loc[i, "LakeArea"] = slakeinfo.iloc[0]["Lake_area"]
            catinfo.loc[i, "LakeDepth"] = slakeinfo.iloc[0]["Depth_avg"]
            catinfo.loc[i, "Laketype"] = slakeinfo.iloc[0]["Lake_type"]

    PERMANENT.close()
    return catinfo
