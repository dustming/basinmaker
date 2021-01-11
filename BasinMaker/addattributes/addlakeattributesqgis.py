from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3


def add_lake_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    catinfo,
):
    catchments = input_geo_names["catchment_without_merging_lakes"]
    sl_connected_lake = input_geo_names["sl_connected_lake"]
    sl_non_connected_lake = input_geo_names["sl_nonconnect_lake"]
    outlet_pt_info = input_geo_names["outlet_pt_info"]
    all_lakes = input_geo_names["all_lakes"]

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
        map=outlet_pt_info,
        raster=sl_non_connected_lake,
        column="ncl",
    )

    grass.run_command(
        "v.what.rast",
        map=outlet_pt_info,
        raster="Connect_Lake_Cat_w_Lake_ID",
        column="cl",
    )

    ### read catchment
    sqlstat = "SELECT SubId,ncl,cl FROM %s" % (outlet_pt_info)
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-9999)
    lakeinfo = Dbf_To_Dataframe(os.path.join(grassdb,all_lakes+'.shp'))

    for i in range(0, len(outletinfo)):
        catid = outletinfo["SubId"].values[i]
        catrow = catinfo["SubId"] == catid

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
            catinfo.loc[catrow, "IsLake"] = 1
            slakeinfo = lakeinfo.loc[lakeinfo["Hylak_id"] == CL_LakeId]
            catinfo.loc[catrow, "HyLakeId"] = CL_LakeId
            catinfo.loc[catrow, "LakeVol"] = slakeinfo.iloc[0]["Vol_total"]
            catinfo.loc[catrow, "LakeArea"] = slakeinfo.iloc[0]["Lake_area"]
            catinfo.loc[catrow, "LakeDepth"] = slakeinfo.iloc[0]["Depth_avg"]
            catinfo.loc[catrow, "Laketype"] = slakeinfo.iloc[0]["Lake_type"]
        if NCL_LakeId > 0 and CL_LakeId < 0:
            catinfo.loc[catrow, "IsLake"] = 2
            slakeinfo = lakeinfo.loc[lakeinfo["Hylak_id"] == NCL_LakeId]
            catinfo.loc[catrow, "HyLakeId"] = NCL_LakeId
            catinfo.loc[catrow, "LakeVol"] = slakeinfo.iloc[0]["Vol_total"]
            catinfo.loc[catrow, "LakeArea"] = slakeinfo.iloc[0]["Lake_area"]
            catinfo.loc[catrow, "LakeDepth"] = slakeinfo.iloc[0]["Depth_avg"]
            catinfo.loc[catrow, "Laketype"] = slakeinfo.iloc[0]["Lake_type"]

    PERMANENT.close()
    return catinfo
