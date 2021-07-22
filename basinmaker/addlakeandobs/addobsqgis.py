from basinmaker.func.grassgis import *
from basinmaker.func.qgis import *
from basinmaker.func.pdtable import *
from basinmaker.func.rarray import *
from basinmaker.utilities.utilities import *
import sqlite3
from basinmaker.preprocessing.preprocessingobs import preprocessing_obs_point


def add_obs_into_existing_watershed_delineation(
    grassdb,
    grass_location,
    qgis_prefix_path,
    input_geo_names,
    path_obsfile_in,
    obs_attributes=[],
    search_radius=100,
    path_sub_reg_outlets_v="#",
    max_memroy=1024 * 4,
    pourpoints_add_obs="pourpoints_add_obs",
    snapped_obs_points="snapped_obs_points",
):

    fdr_arcgis = input_geo_names["fdr_arcgis"]
    fdr_grass = input_geo_names["fdr_grass"]
    str_r = input_geo_names["str_r"]
    str_v = input_geo_names["str_v"]
    acc = input_geo_names["acc"]
    cat_no_lake = input_geo_names["cat_no_lake"]
    mask = input_geo_names["mask"]
    dem = input_geo_names["dem"]
    pourpoints_with_lakes = input_geo_names["pourpoints_with_lakes"]
    lake_outflow_pourpoints = input_geo_names["lake_outflow_pourpoints"]
    cat_add_lake = input_geo_names["cat_add_lake"]

    # define internal file names
    obsname = Internal_Constant_Names["obs"]

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

    catids, temp = generate_stats_list_from_grass_raster(
        grass, mode=1, input_a=pourpoints_with_lakes
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
        output=snapped_obs_points,
        column="#",
        use="cat",
    )
    grass_raster_r_to_vect(
        grass,
        input=snapped_obs_points,
        output=snapped_obs_points,
        type="point",
        flags="v",
    )
    grass_raster_v_db_join(
        grass,
        map=snapped_obs_points,
        column="cat",
        other_table=obsname + "t1",
        other_column="cat",
    )
    exp = obs_attributes[0] + "n int"
    grass.run_command(
        "v.db.addcolumn", map=snapped_obs_points, columns=obs_attributes[0] + "n int"
    )
    grass.run_command(
        "v.db.update",
        map=snapped_obs_points,
        column=obs_attributes[0] + "n",
        qcol=obs_attributes[0] + " + " + str(int(maxcatid) + 1),
    )

    grass.run_command(
        "v.out.ogr",
        input=snapped_obs_points,
        output=os.path.join(grassdb, snapped_obs_points + ".shp"),
        format="ESRI_Shapefile",
        overwrite=True,
    )

    grass_raster_v_to_raster(
        grass,
        input=snapped_obs_points,
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
            column="reg_subid",
            use="attr",
        )
        # added into observation raster point
        exp = "%s = if(isnull(int(Sub_reg_outlets)),%s,Sub_reg_outlets)" % (
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
    if lake_outflow_pourpoints != "#":

        ##### obtain lake id and correspond catchment id
        lake_id, cat_id = generate_stats_list_from_grass_raster(
            grass,
            mode=2,
            input_a=lake_outflow_pourpoints,
            input_b=pourpoints_with_lakes,
        )

        lake_new_cat_ids = np.column_stack((lake_id, cat_id))
        grass.run_command("g.copy", rast=(obsname, obsname + "2"), overwrite=True)
        # remove obs that located within the lake catchments
        obsid, cat_add_lake_id = generate_stats_list_from_grass_raster(
            grass, mode=2, input_a=obsname, input_b=cat_add_lake
        )
        lakecat_obs = np.column_stack((cat_add_lake_id, obsid))
        obsinlake_mask = np.isin(lakecat_obs[:, 0], lake_new_cat_ids[:, 1])
        obsid_inlake = lakecat_obs[obsinlake_mask, 1]
        if len(obsid_inlake) > 0:
            grass.run_command(
                "r.null", map=obsname + "2", setnull=obsid_inlake, overwrite=True
            )

        # combine lake and obs pourpoints
        # combine obsoutlets and outlet from cat no lake
        exp = "'%s' =if(isnull(int(%s)),%s,%s)" % (
            pourpoints_add_obs,
            pourpoints_with_lakes,
            obsname + "2",
            pourpoints_with_lakes,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    else:
        # combine obsoutlets and outlet from cat no lake
        exp = "'%s' =if(isnull(int(%s)),%s,%s)" % (
            pourpoints_add_obs,
            pourpoints_with_lakes,
            obsname,
            pourpoints_with_lakes,
        )
        grass.run_command("r.mapcalc", expression=exp, overwrite=True)

    return
