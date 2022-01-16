from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
from basinmaker.preprocessing.preprocessinglakeply import preprocessing_lake_polygon
import sqlite3
from basinmaker.addlakeandobs.definelaketypeqgis import (
    define_connected_and_non_connected_lake_type,
)
from basinmaker.addlakeandobs.filterlakesqgis import select_lakes_by_area_r
from basinmaker.addlakeandobs.pourpointsqgis import define_pour_points_with_lakes
from basinmaker.addlakeandobs.modifyfdr import modify_lakes_flow_direction
from basinmaker.addlakeandobs.modifystr import modify_str_r_to_add_sub_in_head_str


def add_lakes_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_lakefile_in,
    lake_attributes,
    threshold_con_lake,
    threshold_non_con_lake,
    only_included_lake_at_river_interction = False,
    remove_lake_inlets = False,
    path_sub_reg_lake_r="#",
    path_sub_reg_lake_bd_r="#",
    sl_connected_lake="sl_connected_lake",
    sl_non_connected_lake="sl_nonconnect_lake",
    sl_lakes="selected_lakes",
    sl_str_connected_lake="str_sl_connected_lake",
    nfdr_arcgis="narcgis_fdr",
    nfdr_grass="ngrass_fdr",
    cat_add_lake="cat_add_lake",
    pourpoints_with_lakes="pourpoints_with_lakes",
    cat_use_default_acc="cat_use_default_acc",
    lake_outflow_pourpoints="lake_outflow_pourpoints",
    problem_seg="problem_seg",
    max_memroy=1024 * 4,
):
    # define required input files names
    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]

    # define internal file names
    lake_inflow_pourpoints = Internal_Constant_Names["lake_inflow_pourpoints"]
    catchment_pourpoints_outside_lake = Internal_Constant_Names[
        "catchment_pourpoints_outside_lake"
    ]
    cat_add_lake_old_fdr = Internal_Constant_Names["cat_add_lake_old_fdr"]
    str_connected_lake = Internal_Constant_Names["str_connected_lake"]
    alllake = Internal_Constant_Names["all_lakes"]
    lake_boundary = Internal_Constant_Names["lake_boundary"]
    connected_lake = Internal_Constant_Names["connect_lake"]
    non_connected_lake = Internal_Constant_Names["nonconnect_lake"]
    lakes_lg_cl_thres = "lakes_lg_cl_thres"
    lakes_lg_ncl_thres = "lakes_lg_ncl_thres"
    # prepropessing lakes inputs
    if path_lakefile_in == "#":
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

        routing_info = generate_routing_info_of_catchments(
            grass, con, cat=cat_no_lake, acc=acc, str=str_r, Name="cat1",garray=garray
        )
        grass.run_command(
            "g.copy", rast=("cat1_OL", pourpoints_with_lakes), overwrite=True
        )

        lake_outflow_pourpoints = '#'
        grass.run_command(
            "g.copy", rast=(cat_no_lake, cat_use_default_acc), overwrite=True
        )
        grass.run_command(
            "g.copy", rast=(cat_no_lake, cat_add_lake), overwrite=True
        )
        grass.run_command("g.copy", rast=(fdr_grass, nfdr_grass), overwrite=True)
        grass.run_command("g.copy", rast=(fdr_arcgis, nfdr_arcgis), overwrite=True)
        return lake_outflow_pourpoints

    preprocessing_lake_polygon(
        path_lakefile_in=path_lakefile_in,
        lake_attributes=lake_attributes,
        mask=mask,
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        threshold_con_lake=threshold_con_lake,
        threshold_non_con_lake=threshold_non_con_lake,
        gis_platform="qgis",
        lake_name=alllake,
        lake_boundary_name=lake_boundary,
        lakes_lg_cl_thres=lakes_lg_cl_thres,
        lakes_lg_ncl_thres=lakes_lg_ncl_thres,
    )

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

    if path_sub_reg_lake_r != "#":
        grass.run_command(
            "r.unpack", input=path_sub_reg_lake_r, output="rg_lake", overwrite=True
        )
        exp = "%s = int(%s)" % ("rg_lake", "rg_lake")
        grass_raster_r_mapcalc(grass, exp)
    
        grass.run_command(
            "r.unpack",
            input=path_sub_reg_lake_bd_r,
            output="rg_lake_bd",
            overwrite=True,
        )
        exp = "%s = int(%s)" % ("rg_lake_bd", "rg_lake_bd")
        grass_raster_r_mapcalc(grass, exp)
    
        exp = "%s = if(isnull(%s),%s,%s)" % (alllake,"rg_lake",alllake,"rg_lake")
        grass_raster_r_mapcalc(grass, exp)
        exp = "%s = if(isnull(%s),%s,%s)" %  (lake_boundary,"rg_lake_bd",lake_boundary,"rg_lake")
        exp = "%s = int(%s)" % (alllake, alllake)
        grass_raster_r_mapcalc(grass, exp)
        exp = "%s = int(%s)" % (lake_boundary,lake_boundary)
        grass_raster_r_mapcalc(grass, exp)
    
    else:
        write_grass_and_arcgis_fdr_rules(grassdb)
        exp = "%s = int(%s)" % (alllake, alllake)
        grass_raster_r_mapcalc(grass, exp)
        exp = "%s = int(%s)" % (lake_boundary, lake_boundary)
        grass_raster_r_mapcalc(grass, exp)
    
    # Define connected and non connected lakes and
    # identify which str make certain lake have two outlet
    define_connected_and_non_connected_lake_type(
        grass,
        con,
        garray,
        str_r=str_r,
        lake=alllake,
        connected_lake=connected_lake,
        non_connected_lake=non_connected_lake,
        str_connected_lake=str_connected_lake,
    )
    
    select_lakes_by_area_r(
        grass=grass,
        con=con,
        garray=garray,
        str_r=str_r,
        threshold_con_lake=threshold_con_lake,
        threshold_non_con_lake=threshold_non_con_lake,
        only_included_lake_at_river_interction = only_included_lake_at_river_interction,
        lake_attributes=lake_attributes,
        lakes=alllake,
        connected_lake=connected_lake,
        non_connected_lake=non_connected_lake,
        str_connected_lake=str_connected_lake,
        sl_connected_lake=sl_connected_lake,
        sl_non_connected_lake=sl_non_connected_lake,
        sl_lakes=sl_lakes,
        sl_str_connected_lake=sl_str_connected_lake,
    )
    
    modify_str_r_to_add_sub_in_head_str(
        grass = grass,
        garray = garray,
        con = con,
        str_r = str_r,
        str_v = str_v,
        connected_lake = connected_lake,
        cat_no_lake = cat_no_lake,
        fdr_grass = fdr_grass,
        max_memroy = max_memroy,      
    )

    Lakes_WIth_Multi_Outlet, Remove_Str = define_pour_points_with_lakes(
        grass=grass,
        con=con,
        garray=garray,
        str_r=str_r,
        remove_lake_inlets = remove_lake_inlets,
        cat_no_lake=cat_no_lake,
        sl_lakes=sl_lakes,
        sl_connected_lake=sl_connected_lake,
        sl_str_connected_lake=sl_str_connected_lake,
        acc=acc,
        pourpoints_with_lakes=pourpoints_with_lakes,
        lake_inflow_pourpoints=lake_inflow_pourpoints,
        lake_outflow_pourpoints=lake_outflow_pourpoints,
        catchment_pourpoints_outside_lake=catchment_pourpoints_outside_lake,
    )

    grass.run_command(
        "r.stream.basins",
        direction=fdr_grass,
        points=pourpoints_with_lakes,
        basins=cat_add_lake_old_fdr,
        overwrite=True,
        memory=max_memroy,
    )

    grass.run_command(
        "v.what.rast",
        map=lake_outflow_pourpoints,
        raster=acc,
        column="lmax_acc",
    )
    ### read catchment
    sqlstat = "SELECT cat,lmax_acc FROM %s" % (lake_outflow_pourpoints,)
    lakeinfo = pd.read_sql_query(sqlstat, con)
    lakeinfo = lakeinfo.fillna(-9999)
    lakeinfo = lakeinfo.loc[lakeinfo["cat"] > 0]


    cat_withlake_array = garray.array(mapname=cat_add_lake_old_fdr)
    cat_withlake_array = cat_withlake_array.astype(int)
    fdr_arcgis_array = garray.array(mapname=fdr_arcgis)
    fdr_arcgis_array = fdr_arcgis_array.astype(int)
    str_r_array = garray.array(mapname=str_r)
    sl_lakes_array = garray.array(mapname=sl_lakes)
    sl_lakes_array = sl_lakes_array.astype(int)
    acc_array = garray.array(mapname=acc)
    ncols = int(cat_withlake_array.shape[1])
    nrows = int(cat_withlake_array.shape[0])
    lake_boundary_array = garray.array(mapname=lake_boundary)

    maximumLakegrids = 1000000000
    pec_grid_outlier = 1
    un_modify_fdr_lakeids = []
    outlakeids, chandir, ndir, bd_problem = modify_lakes_flow_direction(
        cat_withlake_array,
        sl_lakes_array,
        acc_array,
        fdr_arcgis_array,
        str_r_array,
        lakeinfo,
        nrows,
        ncols,
        lake_boundary_array,
        pec_grid_outlier,
        maximumLakegrids,
        un_modify_fdr_lakeids,
    )

    temparray = garray.array()

    temparray[:, :] = ndir[:, :]
    temparray.write(mapname=nfdr_arcgis, overwrite=True)
    grass.run_command("r.null", map=nfdr_arcgis, setnull=-9999)

    temparray[:, :] = chandir[:, :]
    temparray.write(mapname="chandir", overwrite=True)
    grass.run_command("r.null", map="chandir", setnull=-9999)

    temparray[:, :] = bd_problem[:, :]
    temparray.write(mapname="bd_problem", overwrite=True)
    grass.run_command("r.null", map="bd_problem", setnull=-9999)

    grass.run_command(
        "r.reclass",
        input=nfdr_arcgis,
        output=nfdr_grass,
        rules=os.path.join(grassdb, "Arcgis2GrassDIR.txt"),
        overwrite=True,
    )
    # cat5
    grass.run_command(
        "r.stream.basins",
        direction=nfdr_grass,
        points=pourpoints_with_lakes,
        basins=cat_add_lake,
        overwrite=True,
        memory=max_memroy,
    )

    grass.run_command("g.copy", rast=(cat_no_lake, cat_use_default_acc), overwrite=True)
    grass.run_command("g.copy", rast=(str_r, "good_seg"), overwrite=True)

    if len(Remove_Str) > 0:
        grass.run_command(
            "r.null", map=cat_use_default_acc, setnull=Remove_Str, overwrite=True
        )
        grass.run_command("r.null", map="good_seg", setnull=Remove_Str, overwrite=True)

        exp = "prob_seg_str = if(isnull(good_seg),%s,null())" % (str_r)
        grass_raster_r_mapcalc(grass, expression=exp)

        exp = "prob_seg_str_w_lid = if(isnull(prob_seg_str),null(),%s)" % (
            sl_connected_lake
        )
        grass_raster_r_mapcalc(grass, expression=exp)

        grass.run_command(
            "g.copy", rast=("prob_seg_str_w_lid", "prob_seg_lake"), overwrite=True
        )

        grass.run_command(
            "r.null",
            map="prob_seg_lake",
            setnull=Lakes_WIth_Multi_Outlet,
            overwrite=True,
        )

        exp = "%s = if(isnull(prob_seg_lake),prob_seg_str_w_lid,null())" % (problem_seg)
        grass_raster_r_mapcalc(grass, expression=exp)

    else:
        exp = "%s = if(isnull(%s),null(),null())" % (problem_seg, "good_seg")
        grass_raster_r_mapcalc(grass, expression=exp)

    grass.run_command(
        "r.to.vect",
        input=problem_seg,
        output=problem_seg,
        type="point",
        overwrite=True,
        flags="v",
    )

    print("Following lake have multi outlet ")
    print(Lakes_WIth_Multi_Outlet)
    print("following str are corrected to make one lake one outlet")
    print(Remove_Str)

    grass.run_command(
        "v.out.ogr",
        input=catchment_pourpoints_outside_lake,
        output=os.path.join(grassdb, catchment_pourpoints_outside_lake + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )
    grass.run_command(
        "v.out.ogr",
        input=pourpoints_with_lakes,
        output=os.path.join(grassdb, pourpoints_with_lakes + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )

    grass.run_command(
        "v.out.ogr",
        input=lake_outflow_pourpoints,
        output=os.path.join(grassdb, lake_outflow_pourpoints + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )
    grass.run_command(
        "v.out.ogr",
        input=lake_inflow_pourpoints,
        output=os.path.join(grassdb, lake_inflow_pourpoints + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
        quiet="Ture",
    )

    PERMANENT.close()
    return lake_outflow_pourpoints
