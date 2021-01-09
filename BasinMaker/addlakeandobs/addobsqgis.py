from func.grassgis import *
from func.qgis import *
from func.pdtable import *
from func.rarray import *
from utilities.utilities import *
import sqlite3
from preprocessing.preprocessingobs import preprocessing_obs_point
from addlakeandobs.definelaketypeqgis import generate_stats_list_from_grass_raster


def add_obs_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_obsfile_in,
    path_lakefile_in,
    obs_attributes=[],
    search_radius=100,
    lake_pourpoints="#",
    path_sub_reg_outlets_v="#",
    max_memroy=1024 * 4,
    obsname="obs",
    pourpoints_add_obs="pourpoints_add_obs",
):

    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]

    # prepropessing lakes inputs
    preprocessing_obs_point(
        mask=mask,
        path_obsin_in=path_obsfile_in,
        obs_attributes=obs_attributes,
        grassdb=grassdb,
        grass_location=grass_location,
        qgis_prefix_path=qgis_prefix_path,
        gis_platform="qgis",
        obsname=obsname + "t1",
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

    # obtain maximum current cat id
    if path_lakefile_in == "#":
        catids, temp = generate_stats_list_from_grass_raster(
            grass, mode=1, input_a=lake_pourpoints
        )
        maxcatid = max(catids)
    else:
        # find catchment without lake outlets
        grass.run_command(
            "r.stats.zonal",
            base=cat_no_lake,
            cover=acc,
            method="max",
            output="catnolake_maxacc",
            overwrite=True,
        )
        ### Find the grid that equal to the max acc, thus this is the outlet grids
        exp = "'%s' =if(%s == int('%s'),%s,null())" % (
            "catnolake_OL",
            acc,
            "catnolake_maxacc",
            cat_no_lake,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

        catids, temp = generate_stats_list_from_grass_raster(
            grass, mode=1, input_a="catnolake_OL"
        )
        maxcatid = max(catids)

    # snap obs points
    grass_raster_r_stream_snap(
        grass,
        input=obsname + "t1",
        output=obsname + "_snap",
        stream_rast=str_r,
        accumulation=acc,
        radius=search_radius,
        memory=max_memroy,
    )

    # obtain use Obs_ID as observation point raster value
    grass_raster_v_to_raster(
        grass,
        input=obsname + "_snap",
        output=obsname + "_snap",
        column="#",
        use="cat",
    )
    grass_raster_r_to_vect(
        grass,
        input=obsname + "_snap",
        output=obsname + "_snap_r2v",
        type="point",
        flags="v",
    )
    grass_raster_v_db_join(
        grass,
        map=obsname + "_snap_r2v",
        column="cat",
        other_table=obsname + "t1",
        other_column="cat",
    )
    exp = obs_attributes[0] + "n int"
    grass.run_command(
        "v.db.addcolumn", map=obsname + "_snap_r2v", columns=obs_attributes[0] + "n int"
    )
    grass.run_command(
        "v.db.update",
        map=obsname + "_snap_r2v",
        column=obs_attributes[0] + "n",
        qcol=obs_attributes[0] + " + " + str(int(maxcatid) + 1),
    )

    grass.run_command(
        "v.out.ogr",
        input=obsname + "_snap_r2v",
        output=os.path.join(grassdb, obsname + "_snap_r2v.shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )

    grass_raster_v_to_raster(
        grass,
        input=obsname + "_snap_r2v",
        output=obsname + "1",
        column=obs_attributes[0] + "n",
        use="attr",
    )

    if path_sub_reg_outlets_v != "#":
        # unpack subregion outlet point
        grass_raster_v_unpack(
            grass, input=path_sub_reg_outlets_v, output="Sub_reg_outlets_pt"
        )
        # convert it to raster
        grass_raster_v_to_raster(
            grass,
            input="Sub_reg_outlets_pt",
            output="Sub_reg_outlets",
            column="value",
            use="attr",
        )
        # added into observation raster point
        exp = "%s = if(isnull(Sub_reg_outlets),%s,Sub_reg_outlets)" % (
            obsname,
            obsname + "1",
        )
        grass_raster_r_mapcalc(
            grass,
            expression=exp,
        )

        grass_raster_setnull(
            grass,
            raster_nm=obsname,
            null_values=[-9999, 0],
            create_new_raster=False,
            new_raster_nm="#",
        )
    else:
        grass.run_command("g.copy", rast=(obsname + "1", obsname), overwrite=True)
    ####

    # remove obs point located in lake catchments
    if path_lakefile_in != "#":

        ##### obtain lake id and correspond catchment id
        lake_id, cat_id = generate_stats_list_from_grass_raster(
            grass, mode=2, input_a="pourpoints_sl_lakes", input_b=lake_pourpoints
        )

        lake_new_cat_ids = np.column_stack((lake_id, cat_id))

        # remove obs that located within the lake catchments
        obsid, cat_add_lake_id = generate_stats_list_from_grass_raster(
            grass, mode=2, input_a=obsname, input_b=lake_pourpoints
        )
        lakecat_obs = np.column_stack((cat_add_lake_id, obsid))
        obsinlake_mask = np.isin(lakecat_obs[:, 0], lake_new_cat_ids[:, 1])
        obsid_inlake = lakecat_obs[obsinlake_mask, 1]
        if len(obsid_inlake) > 0:
            grass.run_command(
                "r.null", map=obsname, setnull=obsid_inlake, overwrite=True
            )

        # combine lake and obs pourpoints

        # combine obsoutlets and outlet from cat no lake
        exp = "'%s' =if(isnull(%s),%s,%s)" % (
            pourpoints_add_obs,
            lake_pourpoints,
            obsname,
            lake_pourpoints,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    else:
        # combine obsoutlets and outlet from cat no lake
        exp = "'%s' =if(isnull(%s),null(),%s)" % (
            pourpoints_add_obs,
            "catnolake_OL",
            obsname,
            "catnolake_OL",
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    return
