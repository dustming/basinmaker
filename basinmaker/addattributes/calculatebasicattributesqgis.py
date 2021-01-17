from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3


def calculate_basic_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    projection,
    catinfo,
    cat_ply_info,
    cat_riv_info,
    outlet_pt_info,
):

    catchments = input_geo_names["catchment_without_merging_lakes"]
    river_r = input_geo_names["river_without_merging_lakes"]
    river_v = input_geo_names["str_v"]
    dem = input_geo_names["dem"]
    fdr = input_geo_names["nfdr_grass"]
    acc = input_geo_names["acc"]
    cat_use_default_acc = input_geo_names["cat_use_default_acc"]

    import grass.script as grass
    import grass.script.setup as gsetup
    from grass.pygrass.modules import Module
    from grass.pygrass.modules.shortcuts import general as g
    from grass.pygrass.modules.shortcuts import raster as r
    from grass.script import array as garray
    from grass.script import core as gcore
    from grass_session import Session

    QgsApplication.setPrefixPath(qgis_prefix_path, True)
    Qgs = QgsApplication([], False)
    Qgs.initQgis()
    from processing.core.Processing import Processing
    from processing.tools import dataobjects
    from qgis import processing

    feedback = QgsProcessingFeedback()
    Processing.initialize()
    QgsApplication.processingRegistry().addProvider(QgsNativeAlgorithms())
    context = dataobjects.createContext()
    context.setInvalidGeometryCheck(QgsFeatureRequest.GeometryNoCheck)

    os.environ.update(
        dict(GRASS_COMPRESS_NULLS="1", GRASS_COMPRESSOR="ZSTD", GRASS_VERBOSE="1")
    )
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")

    con = sqlite3.connect(
        os.path.join(grassdb, grass_location, "PERMANENT", "sqlite", "sqlite.db")
    )

    # create a catchment vector and overlay with str v
    exp = "%s = int(%s)" % (catchments, catchments)
    grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    grass.run_command(
        "r.to.vect",
        input=catchments,
        output=cat_ply_info + "t1",
        type="area",
        overwrite=True,
    )
    grass.run_command(
        "v.db.addcolumn", map=cat_ply_info + "t1", columns="GC_str VARCHAR(40)"
    )
    grass.run_command(
        "v.db.addcolumn", map=cat_ply_info + "t1", columns="Area_m double"
    )
    grass.run_command(
        "v.db.update", map=cat_ply_info + "t1", column="GC_str", qcol="value"
    )
    # dissolve based on gridcode
    grass.run_command(
        "v.dissolve",
        input=cat_ply_info + "t1",
        column="GC_str",
        output=cat_ply_info,
        overwrite=True,
    )

    grass.run_command("v.db.addcolumn", map=cat_ply_info, columns="Gridcode INT")
    grass.run_command("v.db.update", map=cat_ply_info, column="Gridcode", qcol="GC_str")

    ## obtain a stream vector, segmentation based on new catchment polygon
    grass.run_command(
        "v.overlay",
        ainput=river_v,
        alayer=2,
        atype="line",
        binput=cat_ply_info,
        operator="and",
        output=cat_riv_info + "t1",
        overwrite=True,
    )
    grass.run_command(
        "v.out.ogr",
        input=cat_riv_info + "t1",
        output=os.path.join(grassdb, cat_riv_info + "t1" + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )
    processing.run(
        "gdal:dissolve",
        {
            "INPUT": os.path.join(grassdb, cat_riv_info + "t1" + ".shp"),
            "FIELD": "b_GC_str",
            "OUTPUT": os.path.join(grassdb, cat_riv_info + ".shp"),
        },
    )
    grass.run_command(
        "v.import",
        input=os.path.join(grassdb, cat_riv_info + ".shp"),
        output=cat_riv_info,
        overwrite=True,
    )

    grass.run_command("v.db.addcolumn", map=cat_riv_info, columns="Gridcode INT")
    grass.run_command("v.db.addcolumn", map=cat_riv_info, columns="Length_m double")
    grass.run_command(
        "v.db.update", map=cat_riv_info, column="Gridcode", qcol="b_GC_str"
    )
    grass.run_command("v.db.dropcolumn", map=cat_riv_info, columns=["b_GC_str"])

    grass.run_command(
        "r.out.gdal",
        input=dem,
        output=os.path.join(grassdb, "dem_par.tif"),
        format="GTiff",
        overwrite=True,
    )

    PERMANENT.close()

    os.system(
        "gdalwarp "
        + '"'
        + os.path.join(grassdb, "dem_par.tif")
        + '"'
        + "    "
        + '"'
        + os.path.join(grassdb, "dem_proj.tif")
        + '"'
        + " -t_srs  "
        + '"'
        + projection
        + '"'
    )

    project = Session()
    project.open(
        gisdb=grassdb, location=grass_location + "_proj", create_opts=projection
    )

    grass.run_command(
        "r.import",
        input=os.path.join(grassdb, "dem_proj.tif"),
        output="dem_proj",
        overwrite=True,
    )
    grass.run_command("g.region", raster="dem_proj")
    grass.run_command(
        "v.proj",
        location=grass_location,
        mapset="PERMANENT",
        input=cat_riv_info,
        overwrite=True,
    )
    grass.run_command(
        "v.proj",
        location=grass_location,
        mapset="PERMANENT",
        input=cat_ply_info,
        overwrite=True,
    )

    grass.run_command(
        "v.to.db",
        map=cat_ply_info,
        option="area",
        columns="Area_m",
        units="meters",
        overwrite=True,
    )
    grass.run_command(
        "v.to.db",
        map=cat_riv_info,
        option="length",
        columns="Length_m",
        units="meters",
        overwrite=True,
    )

    grass.run_command(
        "r.slope.aspect",
        elevation="dem_proj",
        slope="slope",
        aspect="aspect",
        precision="DCELL",
        overwrite=True,
    )

    ### calcuate averaged DEM in each subbasin
    grass.run_command(
        "v.rast.stats",
        map=cat_ply_info,
        raster="dem_proj",
        column_prefix="d",
        method="average",
    )
    ### calcuate averaged slope and aspect of each subbasin
    grass.run_command(
        "v.rast.stats",
        map=cat_ply_info,
        raster="slope",
        column_prefix="s",
        method="average",
    )
    grass.run_command(
        "v.rast.stats",
        map=cat_ply_info,
        raster="aspect",
        column_prefix="a",
        method="average",
    )
    ### calcuate minimum and maximum dem along the channel
    grass.run_command(
        "v.rast.stats",
        map=cat_riv_info,
        raster="dem_proj",
        column_prefix="d",
        method=["minimum", "maximum"],
    )

    project.close()

    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
    grass.run_command("g.region", raster="dem")
    grass.run_command(
        "v.proj",
        location=grass_location + "_proj",
        mapset="PERMANENT",
        input=cat_riv_info,
        overwrite=True,
    )
    grass.run_command(
        "v.proj",
        location=grass_location + "_proj",
        mapset="PERMANENT",
        input=cat_ply_info,
        overwrite=True,
    )

    ## get routing structure

    exp = "%s = int(%s)" % (fdr, fdr)
    grass.run_command(
        "r.mapcalc",
        expression=exp,
        overwrite=True,
    )

    grass.run_command(
        "r.accumulate",
        direction=fdr,
        accumulation="acc_grass_CatOL",
        overwrite=True,
    )

    exp = "acc_grass_CatOL2 = if(isnull(%s),%s,%s)" % (
        cat_use_default_acc,
        "acc_grass_CatOL",
        acc,
    )
    grass.run_command(
        "r.mapcalc",
        expression=exp,
        overwrite=True,
    )

    routing_temp = generate_routing_info_of_catchments(
        grass,
        con,
        cat=catchments,
        acc=acc,#"acc_grass_CatOL2",
        Name="Final",
        str=river_r,
    )

    grass.run_command("g.copy", vector=("Final_OL_v", outlet_pt_info), overwrite=True)
    ### update dataframe

    sqlstat = "SELECT Gridcode, Length_m, d_minimum, d_maximum FROM %s" % (cat_riv_info)
    leninfo = pd.read_sql_query(sqlstat, con)
    ### read catchment
    sqlstat = "SELECT Gridcode, Area_m,d_average,s_average,a_average FROM %s" % (
        cat_ply_info
    )
    areainfo = pd.read_sql_query(sqlstat, con)

    ### read catchment
    sqlstat = "SELECT SubId, DowSubId,ILSubIdmax,ILSubIdmin FROM %s" % (outlet_pt_info)
    outletinfo = pd.read_sql_query(sqlstat, con)

    outletinfo = outletinfo.fillna(-9999)
    leninfo = leninfo.astype(float).fillna(-9999)
    areainfo = areainfo.astype(float).fillna(-9999)

    for i in range(0, len(outletinfo)):
        catid = outletinfo["SubId"].values[i]
        DownSubID = outletinfo["DowSubId"].values[i]
        catinfo.loc[i, "SubId"] = catid
        ### change the downsub id to -1 for watershed outlet
        if len(outletinfo.loc[outletinfo["SubId"] == DownSubID]) < 1:
            catinfo.loc[i, "DowSubId"] = -1
        elif catid == DownSubID:
            catinfo.loc[i, "DowSubId"] = -1
        else:
            catinfo.loc[i, "DowSubId"] = DownSubID

        catarea = np.unique(
            areainfo.loc[areainfo["Gridcode"] == catid]["Area_m"].values
        )  #'Area_m'
        catslope = np.unique(
            areainfo.loc[areainfo["Gridcode"] == catid]["s_average"].values
        )
        catelev = np.unique(
            areainfo.loc[areainfo["Gridcode"] == catid]["d_average"].values
        )
        cataspect = np.unique(
            areainfo.loc[areainfo["Gridcode"] == catid]["a_average"].values
        )
        if len(catarea) == 1:
            catinfo.loc[i, "BasArea"] = catarea
            catinfo.loc[i, "BasSlope"] = catslope
            catinfo.loc[i, "BasAspect"] = cataspect
            catinfo.loc[i, "MeanElev"] = catelev
        else:
            print(
                "Warning  basin area of stream  ",
                catid,
                "   need check   ",
                len(catarea),
            )
            catinfo.loc[i, "BasArea"] = -9999
            catinfo.loc[i, "BasSlope"] = -9999
            catinfo.loc[i, "BasAspect"] = -9999
            catinfo.loc[i, "MeanElev"] = -9999

        ### add river parameters
        rivlen = np.unique(
            leninfo.loc[leninfo["Gridcode"] == catid]["Length_m"].values
        )  #'Area_m'
        dmaxelev = np.unique(
            leninfo.loc[leninfo["Gridcode"] == catid]["d_maximum"].values
        )
        dminelev = np.unique(
            leninfo.loc[leninfo["Gridcode"] == catid]["d_minimum"].values
        )
        if len(rivlen) == 1:
            catinfo.loc[i, "RivLength"] = rivlen
            maxdem = dmaxelev
            mindem = dminelev
            catinfo.loc[i, "Min_DEM"] = mindem
            catinfo.loc[i, "Max_DEM"] = maxdem
            if rivlen >= 0:
                if max(0, float((maxdem - mindem)) / float(rivlen)) == 0:
                    catinfo.loc[i, "RivSlope"] = -9999
                else:
                    catinfo.loc[i, "RivSlope"] = max(
                        0, float((maxdem - mindem)) / float(rivlen)
                    )
            else:
                catinfo.loc[i, "RivSlope"] = -9999
        else:
            catinfo.loc[i, "RivLength"] = -9999
            catinfo.loc[i, "RivSlope"] = -9999
            catinfo.loc[i, "FloodP_n"] = -9999
            catinfo.loc[i, "Min_DEM"] = -9999
            catinfo.loc[i, "Max_DEM"] = -9999

    PERMANENT.close()
    return catinfo
