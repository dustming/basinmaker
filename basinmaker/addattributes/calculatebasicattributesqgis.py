from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
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
    problem_seg = input_geo_names["problem_seg"]
    nfdr_arcgis = input_geo_names["nfdr_arcgis"]

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
    
    grass.run_command(
        "v.rast.stats",
        map=cat_ply_info,
        raster="dem_proj",
        column_prefix="d",
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
    
    ### calcuate averaged DEM in each subbasin
    # grass.run_command(
    #     "v.rast.stats",
    #     map=cat_ply_info,
    #     raster="dem",
    #     column_prefix="d",
    #     method="average",
    # )
    # 
    # ### calcuate minimum and maximum dem along the channel
    # grass.run_command(
    #     "v.rast.stats",
    #     map=cat_riv_info,
    #     raster="dem",
    #     column_prefix="d",
    #     method=["minimum", "maximum"],
    # )
    
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

    ##### obtain catment id overlaied with problem seg
    prom_seg_id, cat_pro_id = generate_stats_list_from_grass_raster(
        grass, mode=2, input_a=problem_seg, input_b=catchments
    )
    cat_pro_id = np.unique(cat_pro_id)
    grass.run_command("g.copy", rast=(catchments, cat_use_default_acc), overwrite=True)
    if len(cat_pro_id) > 0:
        grass.run_command(
            "r.null", map=cat_use_default_acc, setnull=cat_pro_id, overwrite=True
        )
    else:
        cat_pro_id = []

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
        garray = garray,
    )

    routing_temp = generate_routing_info_of_catchments(
        grass,
        con,
        cat=river_r,
        acc="acc_grass_CatOL2",
        Name="Friv",
        str=river_r,
        garray = garray,
    )

    cat_array = garray.array(mapname=catchments)
    nfdr_arcgis_array = garray.array(mapname=nfdr_arcgis)
    cat_outlet_array = garray.array(mapname="Final_OL")
    ncols = int(cat_array.shape[1])
    nrows = int(cat_array.shape[0])

    grass.run_command("g.copy", vector=("Final_OL_v", outlet_pt_info), overwrite=True)
    PERMANENT.close()
    
    
    ## add coordinates to outlet point in wgs84 system
    project_wgs84 = Session()
    project_wgs84.open(
        gisdb=grassdb, location=grass_location + "_wgs84", create_opts='EPSG:4326'
    )
    grass.run_command(
        "v.proj",
        location=grass_location,
        mapset="PERMANENT",
        input=outlet_pt_info,
        overwrite=True,
    )
    grass.run_command(
        "v.to.db",
        map = outlet_pt_info,
        type = 'point',
        option = 'coor',
        columns = ['outletLng','outletLat'],
        overwrite=True,
    )  
    project_wgs84.close()

    ## import updated outlet points after add coordinates in wgs84 system
    PERMANENT = Session()
    PERMANENT.open(gisdb=grassdb, location=grass_location, create_opts="")
    grass.run_command("g.region", raster="dem")
    grass.run_command(
        "v.proj",
        location=grass_location + "_wgs84",
        mapset="PERMANENT",
        input=outlet_pt_info,
        overwrite=True,
    )   
 
    ### update dataframe
    ###
    sqlstat = "SELECT Gridcode, Length_m, d_minimum, d_maximum FROM %s" % (cat_riv_info)
    leninfo = pd.read_sql_query(sqlstat, con)
    ### read catchment
    sqlstat = "SELECT Gridcode, Area_m,d_average,s_average,a_average FROM %s" % (
        cat_ply_info
    )
    areainfo = pd.read_sql_query(sqlstat, con)

    ### read catchment
    sqlstat = "SELECT SubId, DowSubId,ILSubIdmax,ILSubIdmin,outletLat, outletLng FROM %s" % (outlet_pt_info)
    outletinfo = pd.read_sql_query(sqlstat, con)
    outletinfo = outletinfo.fillna(-9999)
    outletinfo = outletinfo.loc[outletinfo['SubId'] >= 0]
    ### read catchment
    sqlstat = "SELECT SubId, DSubId_str FROM %s" % ('Friv_OL_v')
    outlet_riv_info = pd.read_sql_query(sqlstat, con)
    outlet_riv_info = outlet_riv_info.fillna(-9999)


    leninfo = leninfo.astype(float).fillna(-9999)
    areainfo = areainfo.astype(float).fillna(-9999)

    for i in range(0, len(outletinfo)):
        catid = outletinfo["SubId"].values[i]
        DownSubID_cat = outletinfo["DowSubId"].values[i]
        outlet_lat = outletinfo["outletLat"].values[i]
        outlet_lon = outletinfo["outletLng"].values[i]
        
        catinfo.loc[i, "SubId"] = catid
        catinfo.loc[i, "outletLat"] = outlet_lat
        catinfo.loc[i, "outletLng"] = outlet_lon

        # load routing info based on str
        outlet_riv_info_i = outlet_riv_info.loc[outlet_riv_info['SubId'] == catid]

        # check if riv outlet exist
        if len(outlet_riv_info_i) < 1:
            DownSubID = DownSubID_cat
            DownSubID_riv = -9999
        else:
           # get down sub id based in river network
           DownSubID_riv = outlet_riv_info_i['DSubId_str'].values[0]

           # may need to modify if downsub id from river != down subid from cat
           if DownSubID_riv != DownSubID_cat:

               # check if DownSubID_cat drainage to DSubId_str
               Down_id_of_downriv_info = outletinfo.loc[outletinfo['SubId'] == DownSubID_riv]
               if len(Down_id_of_downriv_info) > 0:
                   # get down subid of down subid of river
                   DSubId_DSubId_str = Down_id_of_downriv_info['DowSubId'].values[0]
                   # if down subid of down sub id of river  = down sub id from cat
                   if catid == 10762:
                       print(DSubId_DSubId_str,DownSubID_cat)

                   if DSubId_DSubId_str == DownSubID_cat:
                       DownSubID = DownSubID_riv
                   else:
                       DownSubID = DownSubID_cat
               else:
                   DownSubID = DownSubID_cat
           else:
               DownSubID = DownSubID_cat

        if catid in cat_pro_id or DownSubID == catid or DownSubID == -9999:
            downsubid_array = return_subid_of_next_down_stream_grids(
                cat_array, catid, nfdr_arcgis_array, cat_outlet_array, ncols, nrows
            )
            #            print(catid,DownSubID,downsubid_array)
            if downsubid_array > 0:
                DownSubID = downsubid_array

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
        if len(rivlen) > 0:
            catinfo.loc[i, "RivLength"] = rivlen[0]
            maxdem = dmaxelev[0]
            mindem = dminelev[0]
            catinfo.loc[i, "Min_DEM"] = mindem
            catinfo.loc[i, "Max_DEM"] = maxdem
            if rivlen[0] >= 0:
                if max(0, float((maxdem - mindem)) / float(rivlen[0])) == 0:
                    slope_rch = -9999
                else:
                    slope_rch = max(
                        0, float((maxdem - mindem)) / float(rivlen[0])
                    )

                slope_rch = max(slope_rch,min_riv_slope)  
                slope_rch = min(slope_rch,max_riv_slope)
                  
                catinfo.loc[i, "RivSlope"] = slope_rch
                
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
