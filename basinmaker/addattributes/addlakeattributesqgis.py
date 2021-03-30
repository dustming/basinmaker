from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import sqlite3


def add_lake_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    lake_attributes,
    catinfo,
):
    catchments = input_geo_names["catchment_without_merging_lakes"]
    sl_connected_lake = input_geo_names["sl_connected_lake"]
    sl_non_connected_lake = input_geo_names["sl_nonconnect_lake"]
    outlet_pt_info = input_geo_names["outlet_pt_info"]
    all_lakes = input_geo_names["all_lakes"]
    lake_outflow_pourpoints = input_geo_names["lake_outflow_pourpoints"]

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

    # keep large lake ids
    grass.run_command(
        "r.stats.zonal",
        base=catchments,
        cover=sl_connected_lake,
        method="min",
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

    grass.run_command(
        "v.db.addcolumn", map=lake_outflow_pourpoints, columns="lakeid int"
    )

    grass.run_command(
        "v.db.update", map=lake_outflow_pourpoints, column="lakeid", qcol="cat"
    )

    grass.run_command(
        "v.what.rast",
        map=lake_outflow_pourpoints,
        raster=catchments,
        column="SubId",
    )

    ### read catchment
    sqlstat = "SELECT SubId,lakeid FROM %s" % (lake_outflow_pourpoints)
    lakeoutinfo = pd.read_sql_query(sqlstat, con)
    lakeoutinfo = lakeoutinfo.fillna(-9999)

    ### read catchment
    sqlstat = "SELECT SubId,DowSubId,ncl,cl FROM %s" % (outlet_pt_info)
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-9999)
    lakeinfo = Dbf_To_Dataframe(os.path.join(grassdb, all_lakes + ".shp"))

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
            catinfo.loc[catrow, "Lake_Cat"] = 1
            catinfo.loc[catrow, "HyLakeId"] = CL_LakeId
            slakeinfo = lakeinfo.loc[lakeinfo[lake_attributes[0]] == CL_LakeId]
            if len(slakeinfo) > 0:
                catinfo.loc[catrow, "LakeVol"] = slakeinfo.iloc[0][lake_attributes[3]]
                catinfo.loc[catrow, "LakeArea"] = slakeinfo.iloc[0][lake_attributes[2]]*1000*1000
                catinfo.loc[catrow, "LakeDepth"] = slakeinfo.iloc[0][lake_attributes[4]]
                catinfo.loc[catrow, "Laketype"] = slakeinfo.iloc[0][lake_attributes[1]]
        if NCL_LakeId > 0 and CL_LakeId < 0:
            catinfo.loc[catrow, "Lake_Cat"] = 2
            catinfo.loc[catrow, "HyLakeId"] = NCL_LakeId
            slakeinfo = lakeinfo.loc[lakeinfo[lake_attributes[0]] == NCL_LakeId]
            if len(slakeinfo) > 0:
                catinfo.loc[catrow, "LakeVol"] = slakeinfo.iloc[0][lake_attributes[3]]
                catinfo.loc[catrow, "LakeArea"] = slakeinfo.iloc[0][lake_attributes[2]]*1000*1000
                catinfo.loc[catrow, "LakeDepth"] = slakeinfo.iloc[0][lake_attributes[4]]
                catinfo.loc[catrow, "Laketype"] = slakeinfo.iloc[0][lake_attributes[1]]
        if NCL_LakeId > 0 and CL_LakeId > 0:
            if catinfo["RivLength"].values[catrow] <= 0:
                catinfo.loc[catrow, "Lake_Cat"] = 2
                catinfo.loc[catrow, "HyLakeId"] = NCL_LakeId
                slakeinfo = lakeinfo.loc[lakeinfo[lake_attributes[0]] == NCL_LakeId]
                if len(slakeinfo) > 0:
                    catinfo.loc[catrow, "LakeVol"] = slakeinfo.iloc[0][lake_attributes[3]]
                    catinfo.loc[catrow, "LakeArea"] = slakeinfo.iloc[0][lake_attributes[2]]*1000*1000
                    catinfo.loc[catrow, "LakeDepth"] = slakeinfo.iloc[0][lake_attributes[4]]
                    catinfo.loc[catrow, "Laketype"] = slakeinfo.iloc[0][lake_attributes[1]]
            else:
                catinfo.loc[catrow, "Lake_Cat"] = 1
                catinfo.loc[catrow, "HyLakeId"] = CL_LakeId
                slakeinfo = lakeinfo.loc[lakeinfo[lake_attributes[0]] == CL_LakeId]
                if len(slakeinfo) > 0:
                    catinfo.loc[catrow, "LakeVol"] = slakeinfo.iloc[0][lake_attributes[3]]
                    catinfo.loc[catrow, "LakeArea"] = slakeinfo.iloc[0][lake_attributes[2]]*1000*1000
                    catinfo.loc[catrow, "LakeDepth"] = slakeinfo.iloc[0][lake_attributes[4]]
                    catinfo.loc[catrow, "Laketype"] = slakeinfo.iloc[0][lake_attributes[1]]

        #  #check if down subid of the lake cl lake catchments
        # if CL_LakeId > 0:
        #     downsubid = catinfo["DowSubId"].values[catrow][0]
        #     subid_inlakes = outletinfo.loc[outletinfo["cl"] == CL_LakeId]["SubId"].values
        #     # if downsubid is one of the lake subbasin continue, if not check if it is
        #     # the lake outlet
        #     if downsubid in subid_inlakes and downsubid != catid:
        #         continue
        #     i_lakeinfo = lakeoutinfo.loc[lakeoutinfo['lakeid'] == CL_LakeId]
        #     if len(i_lakeinfo) <=0:
        #         continue
        #     lake_outlet_subid = i_lakeinfo['SubId'].values[0]
        #
        #     # if catid not the lake outlet subid
        #     # and this cat is not flowinto the lake catchment
        #     # force it to flow into the lake outlet catchment
        #     if catid != lake_outlet_subid:
        #         catinfo.loc[catrow, "DowSubId"] = lake_outlet_subid
        #         print(CL_LakeId,catid,downsubid,lake_outlet_subid)

    PERMANENT.close()
    return catinfo
