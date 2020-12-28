from processing_functions_raster_array import *
from processing_functions_raster_grass import *
from processing_functions_raster_qgis import *
from processing_functions_vector_grass import *
from processing_functions_vector_qgis import *
from utilities import *
import sqlite3
from addlakeandobs.definelaketypeqgis import generate_stats_list_from_grass_raster
from addlakeandobs.defineroutinginfoqgis import generate_routing_info_of_catchments


def calculate_basic_attributes(
    grassdb,
    grass_location,
    qgis_prefix_path,
    catchments,
    river_r,
    river_v,
    dem,
    fdr,
    acc,
    cat_use_default_acc,
    projection,
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
        output="Net_cat_F1",
        type="area",
        overwrite=True,
    )
    grass.run_command("v.db.addcolumn", map="Net_cat_F1", columns="GC_str VARCHAR(40)")
    grass.run_command("v.db.addcolumn", map="Net_cat_F1", columns="Area_m double")
    grass.run_command("v.db.update", map="Net_cat_F1", column="GC_str", qcol="value")
    # dissolve based on gridcode
    grass.run_command(
        "v.dissolve",
        input="Net_cat_F1",
        column="GC_str",
        output="Net_cat_F",
        overwrite=True,
    )

    grass.run_command("v.db.addcolumn", map="Net_cat_F", columns="Gridcode INT")
    grass.run_command("v.db.update", map="Net_cat_F", column="Gridcode", qcol="GC_str")

    ## obtain a stream vector, segmentation based on new catchment polygon
    grass.run_command(
        "v.overlay",
        ainput=river_v,
        alayer=2,
        atype="line",
        binput="Net_cat_F",
        operator="and",
        output="nstr_nfinalcat",
        overwrite=True,
    )
    grass.run_command(
        "v.out.ogr",
        input="nstr_nfinalcat",
        output=os.path.join(grassdb, "nstr_nfinalcat.shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )
    processing.run(
        "gdal:dissolve",
        {
            "INPUT": os.path.join(grassdb, "nstr_nfinalcat.shp"),
            "FIELD": "b_GC_str",
            "OUTPUT": os.path.join(grassdb, "nstr_nfinalcat_F.shp"),
        },
    )
    grass.run_command(
        "v.import",
        input=os.path.join(grassdb, "nstr_nfinalcat_F.shp"),
        output="nstr_nfinalcat_F",
        overwrite=True,
    )

    grass.run_command("v.db.addcolumn", map="nstr_nfinalcat_F", columns="Gridcode INT")
    grass.run_command(
        "v.db.addcolumn", map="nstr_nfinalcat_F", columns="Length_m double"
    )
    grass.run_command(
        "v.db.update", map="nstr_nfinalcat_F", column="Gridcode", qcol="b_GC_str"
    )
    grass.run_command("v.db.dropcolumn", map="nstr_nfinalcat_F", columns=["b_GC_str"])

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
        input="nstr_nfinalcat_F",
        overwrite=True,
    )
    grass.run_command(
        "v.proj",
        location=grass_location,
        mapset="PERMANENT",
        input="Net_cat_F",
        overwrite=True,
    )

    grass.run_command(
        "v.to.db",
        map="Net_cat_F",
        option="area",
        columns="Area_m",
        units="meters",
        overwrite=True,
    )
    grass.run_command(
        "v.to.db",
        map="nstr_nfinalcat_F",
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
        map="Net_cat_F",
        raster="dem_proj",
        column_prefix="d",
        method="average",
    )
    ### calcuate averaged slope and aspect of each subbasin
    grass.run_command(
        "v.rast.stats",
        map="Net_cat_F",
        raster="slope",
        column_prefix="s",
        method="average",
    )
    grass.run_command(
        "v.rast.stats",
        map="Net_cat_F",
        raster="aspect",
        column_prefix="a",
        method="average",
    )
    ### calcuate minimum and maximum dem along the channel
    grass.run_command(
        "v.rast.stats",
        map="nstr_nfinalcat_F",
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
        input="nstr_nfinalcat_F",
        overwrite=True,
    )
    grass.run_command(
        "v.proj",
        location=grass_location + "_proj",
        mapset="PERMANENT",
        input="Net_cat_F",
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
        acc="acc_grass_CatOL2",
        Name="Final",
        str=river_r,
    )

    ### update dataframe

    sqlstat = "SELECT Gridcode, Length_m, d_minimum, d_maximum FROM nstr_nfinalcat_F"
    leninfo = pd.read_sql_query(sqlstat, con)
    ### read catchment
    sqlstat = "SELECT Gridcode, Area_m,d_average,s_average,a_average FROM Net_cat_F"
    areainfo = pd.read_sql_query(sqlstat, con)

    ### read catchment
    sqlstat = "SELECT SubId, DowSubId,ILSubIdmax,ILSubIdmin FROM Final_OL_v"
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
